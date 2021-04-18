################################################################################################################################
# General setup
################################################################################################################################

# Import libraries
import sys
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import romb, simpson, dblquad

# Numerical parameters
array_factor_est = 8
tol_residual_rho = 1e-88
dblquad_epsabs = 1e-3
dblquad_epsrel = 1e-3
num_interp_quad_specialized_ignore_yield = int(1e2)
minimum_exponent = np.log(sys.float_info.min)/np.log(10)

################################################################################################################################
# Deformation application class
################################################################################################################################

class deform_network:

	############################################################################################################################
	# Initialization
	############################################################################################################################

	def __init__(self, F, deformation_type, total_time_in_seconds, single_chain_model, relaxation_function = None, J_sw = 1, \
		max_F_dot = None, max_RAM_usage_in_bytes = None, nondimensional_timestep_suggestion = 1e-2, num_grid_suggestion = 129, \
		interp_kind_2D = 'quintic', use_spatial_grid = True, enumerate_full_arrays = True, \
		ignore_yield = False, ignore_reforming = False):

		# Initialize and retain certain variables
		self.csv_initialized = False
		self.ignore_yield = ignore_yield
		self.ignore_reforming = ignore_reforming
		self.use_spatial_grid = use_spatial_grid
		self.enumerate_full_arrays = enumerate_full_arrays
		self.initialized_single_chain_model = single_chain_model

		########################################################################################################################
		# Deformation
		########################################################################################################################

		# Retain certain variables
		self.F = F
		self.J_sw = J_sw
		self.deformation_type = deformation_type
		self.total_time_in_seconds = total_time_in_seconds

		# Estimate the maximum rate of deformation if not given
		if max_F_dot is None:
			t_temp = np.linspace(0, total_time_in_seconds, int(1e5))
			self.max_F_dot = np.max(np.abs(np.diff(F(t_temp))/np.diff(t_temp)))
		else:
			self.max_F_dot = max_F_dot

		# If gamma_c specified, best method to use may change
		if single_chain_model.gamma_c is None:
			self.use_specialized = False
		else:
			self.use_specialized = True

		########################################################################################################################
		# Spatial discretization or quadrature
		########################################################################################################################

		# Spatial integration using a grid
		if use_spatial_grid is True:

			# Retain the 2D interpolation kind
			self.interp_kind_2D = interp_kind_2D

			# Create the symmetry-conscious spatial grid
			self.num_grid = self.adjust_for_romb(num_grid_suggestion)
			self.z = np.linspace(0, single_chain_model.gamma_TS, self.num_grid)
			self.r = self.z
			self.dz = self.z[1] - self.z[0]
			self.dr = self.dz
			self.Z, self.R = np.meshgrid(self.z, self.r)
			ELL = np.sqrt(self.Z*self.Z + self.R*self.R)

			# Integration element specialized for stress calculation
			self.ELEMENT_stress = self.element_stress(self.Z, self.R, single_chain_model)

			# Adjust normalization of P_A_eq on the grid
			P_A_eq_ELL_non_normalized = np.nan_to_num(single_chain_model.P_A_eq(ELL), nan = 0)
			self.P_A_eq_normalization = self.integral_grid_d_3_xi(P_A_eq_ELL_non_normalized)/single_chain_model.P_A_tot_eq

			# Total reverse reaction rate coefficient on the grid
			P_A_eq_ELL = np.nan_to_num(single_chain_model.P_A_eq(ELL, normalization = self.P_A_eq_normalization), nan = 0)
			self.k_ELL = np.nan_to_num(single_chain_model.k(ELL), nan = 0)
			if single_chain_model.gamma_c is None:
				self.K_hat = self.integral_grid_d_3_xi(self.k_ELL*P_A_eq_ELL)/single_chain_model.P_B_tot_eq
			else:
				self.K_hat = single_chain_model.K_hat

			# Maximum reverse reaction rate coefficient on the grid
			if single_chain_model.gamma_c is None:
				self.max_k_rev = np.max(self.k_ELL*P_A_eq_ELL/single_chain_model.P_B_tot_eq)
			else:
				self.max_k_rev = single_chain_model.max_k_rev

		# Spatial integration using quadrature
		else:

			# Inherit from single_chain_model since will use the same integration scheme
			self.P_A_eq_normalization = 1
			self.k_0 = single_chain_model.k_0
			self.K_hat = single_chain_model.K_hat
			self.max_k_rev = single_chain_model.max_k_rev

		########################################################################################################################
		# Relaxation function related setup
		########################################################################################################################

		if relaxation_function is None:
			self.g_timescale = np.inf
			self.g = lambda t, tau: 1 + 0*t + 0*tau
			self.d_g_d_tau = lambda t, tau: 0*t + 0*tau
			self.g_K_hat = np.exp(minimum_exponent)
		else:
			self.g_timescale = relaxation_function.timescale
			self.g = relaxation_function.g
			self.d_g_d_tau = relaxation_function.d_g_d_tau
			try: 
				self.g_K_hat = single_chain_model.P_A_tot_eq/single_chain_model.P_B_tot_eq/self.g_timescale
			except AttributeError:
				self.g_K_hat = np.exp(minimum_exponent)

		########################################################################################################################
		# Time discretization
		########################################################################################################################

		# Estimate timestep based on the smallest timescales
		timescales = np.array([1/self.max_F_dot, 1/self.K_hat, 1/self.max_k_rev, self.g_timescale, 1/self.g_K_hat])
		estimated_timestep = float(nondimensional_timestep_suggestion*np.min(timescales))

		# Enumerating full arrays requires large memory, so have to do it in chunks rather than over the full history
		if use_spatial_grid is True:

			# Memory considerations
			if max_RAM_usage_in_bytes is None:
				import psutil
				max_RAM_usage_in_bytes = psutil.virtual_memory().available
			max_array_numel = max_RAM_usage_in_bytes/8/array_factor_est
			if enumerate_full_arrays is True:
				max_num_time_chunk = np.floor(np.sqrt(max_array_numel/self.num_grid**2)).astype(int)
			else:
				max_num_time_chunk = np.floor(max_array_numel/self.num_grid**2).astype(int)

			# Chunk history and decrease timestep in order to satisfy memory requirements and Romberg integration
			self.num_chunks = 0
			self.num_time = 2*max_num_time_chunk
			while self.num_time > max_num_time_chunk:
				self.num_chunks += 1
				self.num_time = self.adjust_for_romb(total_time_in_seconds/estimated_timestep/self.num_chunks)
			self.timestep = total_time_in_seconds/self.num_chunks/(self.num_time - 1)

			# Enumerate time limits for each chunk
			t_lims_all = total_time_in_seconds/self.num_chunks*np.arange(0, self.num_chunks + 1, 1)
			self.t_lims = np.zeros((self.num_chunks, 2))
			for index_chunk in range(self.num_chunks):
				self.t_lims[index_chunk,:] = [t_lims_all[index_chunk], t_lims_all[index_chunk + 1]]

		# No memory considerations and corresponding history chunking since will not enumerate full arrays
		else:
			self.num_time = self.adjust_for_romb(total_time_in_seconds/estimated_timestep)
			self.timestep = total_time_in_seconds/(self.num_time - 1)

	############################################################################################################################
	# Function to solve for results over the applied deformation history
	############################################################################################################################

	def solve(self, display_progress = True, csv_directory = None, checkpoint_directory = None):

		# Methods using a grid for spatial integrals
		if self.use_spatial_grid is True:

			# Enumerate full arrays to mimimize computation time, chunking time history to satisfy memory requirements
			if self.enumerate_full_arrays is True:

				# Function to remove an initial point from results
				def remove_initial_point(results):
					results_out = list(results)
					for index in range(len(results)):
						results_out[index] = results[index][1:]
					return tuple(results_out)
				
				# Allocate results
				t = []
				F = []
				P_A_tot = []
				total_rate_break = []
				total_rate_reform = []
				beta_sigma_over_n = []

				# Loop over all chunks in the history
				for index_chunk in range(self.num_chunks):

					# Display progress if opted
					if display_progress is True:
						print('  On chunk', index_chunk + 1, 'of', self.num_chunks, end = '\r')

					# Equilibrium initial distribution for first chunk, or initial distribution from end of previous chunk
					if index_chunk == 0:
						results_chunk, P_A_end = self.compute_results_grid_chunky(self.t_lims[index_chunk,:])
					else:
						results_chunk, P_A_end = self.compute_results_grid_chunky(self.t_lims[index_chunk,:], P_A_0 = P_A_end)

					# Remove repeated points in time among chunks
					if index_chunk > 0 or self.J_sw != 1 or self.initialized_single_chain_model.gamma_c is not None:
						results_chunk = remove_initial_point(results_chunk)

					# Collect results
					t = np.append(t, results_chunk[0])
					F = np.append(F, results_chunk[1])
					P_A_tot = np.append(P_A_tot, results_chunk[2])
					total_rate_break = np.append(total_rate_break, results_chunk[3])
					total_rate_reform = np.append(total_rate_reform, results_chunk[4])
					beta_sigma_over_n = np.append(beta_sigma_over_n, results_chunk[5])
					results = t, F, P_A_tot, total_rate_break, total_rate_reform, beta_sigma_over_n

					# Create a checkpoint .csv if opted after each chunk
					if checkpoint_directory is not None:
						checkpoint(checkpoint_directory).create(t[-1], P_A_end)

			# Refrain from enumerating full arrays but requires no history chunking, takes longer
			else:
				results, P_A_end = self.compute_results_grid([0, self.total_time_in_seconds])

		# Methods using quadrature for spatial integrals
		else:

			# Efficient method in the special case for when k is constant and the single-chain model is infinitely-extensible
			if self.use_specialized is True and self.ignore_yield is True:
				results, P_A_end = self.compute_results_quad_specialized_ignore_yield([0, self.total_time_in_seconds])

			# Method for the general case
			else:
				results, P_A_end = self.compute_results_quad([0, self.total_time_in_seconds])
		
		# Append results to the .csv if opted
		if csv_directory is not None:
			if self.csv_initialized is False: 
				results_csv_initialized = results_csv(csv_directory)
				self.csv_initialized = True
			results_csv_initialized.append(0, results)

		# Create a checkpoint .csv if opted
		if checkpoint_directory is not None:
			checkpoint(checkpoint_directory).create(t[-1], P_A_end)

		# Return the results
		return results

	############################################################################################################################
	# Function for computing results using spatial quadrature
	############################################################################################################################

	def compute_results_quad(self, t_span, P_A_0 = None):

		# Enumerate the time and relative deformation components
		t, F_zz_rel, F_rr_rel = self.enumerate_t_and_F_rel(t_span)
		total_time = t_span[-1] - t_span[0]
		num_time = len(t)

		# Function for the relatively-deformed coordinates
		def ell_rel(z, r, index_t, index_tau):
			return np.sqrt(z*z*F_zz_rel[index_t, index_tau]**2 + r*r*F_rr_rel[index_t, index_tau]**2)

		########################################################################################################################
		# Time-dependent quantities for spatial integrals
		########################################################################################################################

		# Function for the relatively-deformed initial distribution
		if P_A_0 is None:
			def P_A_0_rel_t(z, r, index_t):
				return self.initialized_single_chain_model.P_A_eq(ell_rel(z, r, index_t, 0)/self.J_sw**(1/3), \
					normalization = self.P_A_eq_normalization*self.J_sw)
		else:
			sys.exit('Error: Beginning from a nonequilibrium initial distribution not yet implemented.')

		# Function for the relatively-deformed equilibrium distribution
		def P_A_eq_rel(z, r, index_t, index_tau):
			return self.initialized_single_chain_model.P_A_eq(ell_rel(z, r, index_t, index_tau), \
				normalization = self.P_A_eq_normalization)

		# Function for the relatively-deformed reaction rate coefficient function
		def k_rel(z, r, index_t, index_tau):
			return self.initialized_single_chain_model.k(ell_rel(z, r, index_t, index_tau))

		# Function for the reaction propagator
		def Xi(z, r, index_t, index_tau):
			if index_t == index_tau:
				return 1 + 0*(z + r)
			else:
				integrand = np.zeros((index_t + 1 - index_tau))*(z + r)
				for index_s in range(index_t + 1 - index_tau):
					integrand[:,:,index_s] = k_rel(z, r, index_t, index_s)
				return np.exp(-self.integral_ds(integrand))

		# Function for the homogeneous solution for P_A
		def P_A_h(z, r, index_t):
			return np.nan_to_num(P_A_0_rel_t(z, r, index_t)*Xi(z, r, index_t, 0)*self.g(t[index_t], 0), nan = 0)

		# Function for the integrand of K
		def integrand_K(z, r, index_t, index_tau):
			return np.nan_to_num(P_A_eq_rel(z, r, index_t, index_tau)*Xi(z, r, index_t, index_tau)*( \
				k_rel(z, r, index_t, index_tau)*self.g(t[index_t], t[index_tau]) \
					+ self.d_g_d_tau(t[index_t], t[index_tau])), nan = 0)

		########################################################################################################################
		# Solve the integral equation
		########################################################################################################################

		# Amount of initially-intact chains that have been broken
		P_B_tot_h = np.zeros(num_time)
		for index_t in range(num_time):
			P_B_tot_h[index_t] = 1 - self.integral_quad_d_3_xi(lambda z,r: P_A_h(z, r, index_t))

		# Integral equation only defined when P_B_tot_eq is nonzero
		if self.initialized_single_chain_model.P_B_tot_eq > 0 and self.ignore_reforming is False:

			# Kernel K(t,tau) and right-hand side b(t)
			b = self.initialized_single_chain_model.P_B_tot_eq*P_B_tot_h
			K = np.zeros((num_time, num_time))
			for index_t in range(num_time):
				for index_tau in range(index_t + 1):
					K[index_t, index_tau] = self.integral_quad_d_3_xi( \
						lambda z,r: integrand_K(z, r, index_t, index_tau))/self.initialized_single_chain_model.P_B_tot_eq

			# Successive approximations to retrieve rho(t)
			rho = self.solve_Volterra(K, b, total_time)

			# Total probability of broken chains
			P_B_tot = self.initialized_single_chain_model.P_B_tot_eq*rho

		# Integral equation undefined when P_B_tot_eq = 0; also for rate-independent irreversible breaking
		else:
			P_B_tot = P_B_tot_h

		########################################################################################################################
		# Compute and return the results
		########################################################################################################################

		# Total probability of intact chains
		P_A_tot = 1 - P_B_tot

		# Distribution at the end of the partition
		P_A_end = np.nan

		# Total breaking and reforming rates
		total_rate_reform = self.K_hat*(1 - P_A_tot)
		total_rate_break = np.gradient(P_A_tot)/np.gradient(t) - total_rate_reform

		# Nondimensional stress corresponding to applied the deformation
		beta_sigma_h_over_n = np.zeros(num_time)
		beta_sigma_p_over_n = np.zeros(num_time)
		for index_t in range(num_time):
			beta_sigma_h_over_n[index_t] = self.integral_quad_d_3_xi(lambda z,r: P_A_h(z, r, index_t), element = 'stress')
			if self.initialized_single_chain_model.P_B_tot_eq > 0 and self.ignore_reforming is False:
				integrand_beta_sigma_p_over_n = np.zeros(num_time)
				for index_tau in range(index_t + 1):
					integrand_beta_sigma_p_over_n[index_tau] = rho[index_tau]*self.integral_quad_d_3_xi( \
						lambda z,r: integrand_K(z, r, index_t, index_tau), element = 'stress')
				beta_sigma_p_over_n[index_t] = self.integral_d_tau(integrand_beta_sigma_p_over_n)
		beta_sigma_over_n = beta_sigma_h_over_n + beta_sigma_p_over_n

		# Return results
		results = t, self.F(t), P_A_tot, total_rate_break, total_rate_reform, beta_sigma_over_n
		return results, P_A_end

	############################################################################################################################
	# Function for computing results on a spatial grid
	############################################################################################################################

	def compute_results_grid(self, t_span, P_A_0 = None):

		# Enumerate the time and relative deformation components
		t, F_zz_rel, F_rr_rel = self.enumerate_t_and_F_rel(t_span)
		total_time = t_span[-1] - t_span[0]
		num_time = len(t)

		# Function for the relatively-deformed coordinates
		def ell_rel(index_t, index_tau):
			return np.sqrt(self.Z*self.Z*F_zz_rel[index_t, index_tau]**2 + self.R*self.R*F_rr_rel[index_t, index_tau]**2)

		########################################################################################################################
		# Time-dependent quantities for spatial integrals
		########################################################################################################################

		# Function for the relatively-deformed initial distribution
		if P_A_0 is None:
			def P_A_0_rel_t(index_t):
				return self.initialized_single_chain_model.P_A_eq(ell_rel(index_t, 0)/self.J_sw**(1/3), \
					normalization = self.P_A_eq_normalization*self.J_sw)
		else:
			sys.exit('Error: Beginning from a nonequilibrium initial distribution not yet implemented.')

		# Function for the relatively-deformed equilibrium distribution
		def P_A_eq_rel(index_t, index_tau):
			return self.initialized_single_chain_model.P_A_eq(ell_rel(index_t, index_tau), \
				normalization = self.P_A_eq_normalization)

		# Function for the relatively-deformed reaction rate coefficient function
		def k_rel(index_t, index_tau):
			return self.initialized_single_chain_model.k(ell_rel(index_t, index_tau))

		# Function for the reaction propagator
		def Xi(index_t, index_tau):
			if index_t == index_tau:
				return np.ones((self.num_grid, self.num_grid))
			else:
				integrand = np.zeros((self.num_grid, self.num_grid, index_t + 1 - index_tau))
				for index_s in range(index_t + 1 - index_tau):
					integrand[:,:,index_s] = k_rel(index_t, index_s)
				return np.exp(-self.integral_ds(integrand))

		# Function for the homogeneous solution for P_A
		def P_A_h(index_t):
			out = P_A_0_rel_t(index_t)*Xi(index_t, 0)*self.g(t[index_t], 0)
			out[np.isnan(out) + np.isinf(out)] = 0
			return out

		# Function for the integrand of K
		def integrand_K(index_t, index_tau):
			out = P_A_eq_rel(index_t, index_tau)*Xi(index_t, index_tau)*( \
				k_rel(index_t, index_tau)*self.g(t[index_t], t[index_tau]) + self.d_g_d_tau(t[index_t], t[index_tau]))
			out[np.isnan(out) + np.isinf(out)] = 0
			return out

		########################################################################################################################
		# Solve the integral equation
		########################################################################################################################

		# Amount of initially-intact chains that have been broken
		P_B_tot_h = np.zeros(num_time)
		for index_t in range(num_time):
			P_B_tot_h[index_t] = 1 - self.integral_grid_d_3_xi(P_A_h(index_t))

		# Integral equation only defined when P_B_tot_eq is nonzero
		if self.initialized_single_chain_model.P_B_tot_eq > 0 and self.ignore_reforming is False:

			# Kernel K(t,tau) and right-hand side b(t)
			b = self.initialized_single_chain_model.P_B_tot_eq*P_B_tot_h
			K = np.zeros((num_time, num_time))
			for index_t in range(num_time):
				for index_tau in range(index_t + 1):
					K[index_t, index_tau] = \
						self.integral_grid_d_3_xi(integrand_K(index_t,index_tau))/self.initialized_single_chain_model.P_B_tot_eq

			# Successive approximations to retrieve rho(t)
			rho = self.solve_Volterra(K, b, total_time)

			# Total probability of broken chains
			P_B_tot = self.initialized_single_chain_model.P_B_tot_eq*rho

		# Integral equation undefined when P_B_tot_eq = 0; also for rate-independent irreversible breaking
		else:
			P_B_tot = P_B_tot_h

		########################################################################################################################
		# Compute and return the results
		########################################################################################################################

		# Total probability of intact chains
		P_A_tot = 1 - P_B_tot

		# Distribution at the end of the partition
		P_A_end = np.nan

		# Total breaking and reforming rates
		total_rate_reform = self.K_hat*(1 - P_A_tot)
		total_rate_break = np.gradient(P_A_tot)/np.gradient(t) - total_rate_reform

		# Nondimensional stress corresponding to applied the deformation
		beta_sigma_h_over_n = np.zeros(num_time)
		beta_sigma_p_over_n = np.zeros(num_time)
		for index_t in range(num_time):
			beta_sigma_h_over_n[index_t] = self.integral_grid_d_3_xi(P_A_h(index_t), element = 'stress')
			if self.initialized_single_chain_model.P_B_tot_eq > 0 and self.ignore_reforming is False:
				integrand_beta_sigma_p_over_n = np.zeros(num_time)
				for index_tau in range(index_t + 1):
					integrand_beta_sigma_p_over_n[index_tau] = \
						rho[index_tau]*self.integral_grid_d_3_xi(integrand_K(index_t, index_tau), element = 'stress')
				beta_sigma_p_over_n[index_t] = self.integral_d_tau(integrand_beta_sigma_p_over_n)
		beta_sigma_over_n = beta_sigma_h_over_n + beta_sigma_p_over_n

		# Return results
		results = t, self.F(t), P_A_tot, total_rate_break, total_rate_reform, beta_sigma_over_n
		return results, P_A_end

	############################################################################################################################
	# Function for computing results on a spatial grid; enumerates full arrays and uses vectorized operations
	############################################################################################################################

	def compute_results_grid_chunky(self, t_span, P_A_0 = None):

		# Enumerate the time and relative deformation components
		t, F_zz_rel, F_rr_rel = self.enumerate_t_and_F_rel(t_span)
		total_time = t_span[-1] - t_span[0]
		num_time = len(t)
		Delta_t = t - t[0]

		# Enumerate the relatively-deformed coordinates
		Z_rel = np.tensordot(self.Z, F_zz_rel, axes = 0)
		R_rel = np.tensordot(self.R, F_rr_rel, axes = 0)
		ELL_rel = np.sqrt(Z_rel*Z_rel + R_rel*R_rel)
		z_rel_t = np.tensordot(self.z, F_zz_rel[:,0], axes = 0)
		r_rel_t = np.tensordot(self.r, F_rr_rel[:,0], axes = 0)

		# Cleanup
		del F_zz_rel, F_rr_rel, Z_rel, R_rel

		########################################################################################################################
		# Time-dependent quantities for spatial integrals
		########################################################################################################################

		# Enumerate the relatively-deformed equilibrium distribution
		P_A_eq_rel = self.initialized_single_chain_model.P_A_eq(ELL_rel, normalization = self.P_A_eq_normalization)
		
		# Enumerate the relatively-deformed initial distribution
		if P_A_0 is None:
			P_A_0_rel_t = self.initialized_single_chain_model.P_A_eq(ELL_rel[:,:,:,0]/self.J_sw**(1/3), \
				normalization = self.P_A_eq_normalization*self.J_sw)
		else:
			P_A_0_rel_t = np.zeros((self.num_grid, self.num_grid, num_time))
			for index_t in range(num_time):
				P_A_0_rel_t[:,:,index_t] = self.interp_fun_2D(z_rel_t[:,index_t], r_rel_t[:,index_t], P_A_0)

		# Cleanup
		del z_rel_t, r_rel_t

		# Enumerate the relatively-deformed reaction rate coefficient function
		k_rel = self.initialized_single_chain_model.k(ELL_rel)

		# Cleanup
		del ELL_rel

		# Enumerate the reaction propagator
		Xi = np.zeros((self.num_grid, self.num_grid, num_time, num_time))
		for index_t in range(num_time):
			for index_tau in range(index_t + 1):
				Xi[:,:,index_t,index_tau] = \
					np.exp(-self.integral_ds(k_rel[:,:,index_t,index_tau:index_t + 1]))

		# Enumerate the relaxation function and its derivative
		g = np.zeros((num_time, num_time))
		d_g_d_tau = np.zeros((num_time, num_time))
		for index_t in range(num_time):
			for index_tau in range(index_t + 1):
				g[index_t, index_tau] = self.g(Delta_t[index_t], Delta_t[index_tau])
				d_g_d_tau[index_t, index_tau] = self.d_g_d_tau(Delta_t[index_t], Delta_t[index_tau])

		# Homogeneous solution for P_A
		P_A_h = Xi[:,:,:,0]*P_A_0_rel_t

		# Cleanup
		del P_A_0_rel_t

		# Integrand of K
		if self.initialized_single_chain_model.gamma_c is not None:
			k_rel[np.isinf(k_rel)] = 0
		integrand_K = P_A_eq_rel*Xi*(k_rel*g[None,None,:,:] + d_g_d_tau[None,None,:,:])

		# Cleanup
		del Xi, k_rel, P_A_eq_rel

		#
		P_A_h[np.isnan(P_A_h) + np.isinf(P_A_h)] = 0
		integrand_K[np.isnan(integrand_K) + np.isinf(integrand_K)] = 0

		########################################################################################################################
		# Solve the integral equation
		########################################################################################################################

		# Amount of initially-intact chains that have been broken
		P_B_tot_h = 1 - self.integral_grid_d_3_xi(P_A_h)

		# Integral equation only defined when P_B_tot_eq is nonzero
		if self.initialized_single_chain_model.P_B_tot_eq > 0 and self.ignore_reforming is False:

			# Kernel K(t,tau) and right-hand side b(t)
			K = self.integral_grid_d_3_xi(integrand_K)/self.initialized_single_chain_model.P_B_tot_eq
			b = P_B_tot_h/self.initialized_single_chain_model.P_B_tot_eq

			# Successive approximations to retrieve rho(t)
			rho = self.solve_Volterra(K, b, total_time)

			# Total probability of broken chains
			P_B_tot = self.initialized_single_chain_model.P_B_tot_eq*rho

		# Integral equation undefined when P_B_tot_eq = 0; also for rate-independent irreversible breaking
		else:
			P_B_tot = P_B_tot_h
			rho = np.zeros(num_time)

		# Integrand of particular solution for P_A
		integrand_P_A_p = rho*integrand_K

		# Cleanup
		del integrand_K

		########################################################################################################################
		# Compute and return the results
		########################################################################################################################

		# Total probability of intact chains
		P_A_tot = 1 - P_B_tot

		# Distribution of intact chains
		P_A = P_A_h + self.integral_d_tau(integrand_P_A_p)

		# Cleanup
		del P_A_h, integrand_P_A_p

		# Distribution at the end of the partition
		P_A_end = P_A[:,:,-1]

		# Total breaking and reforming rates
		total_rate_break = -self.integral_grid_d_3_xi(P_A*self.k_ELL[:,:,None])
		total_rate_reform = self.K_hat*(1 - P_A_tot)

		# Nondimensional stress corresponding to the applied deformation
		beta_sigma_over_n = self.integral_grid_d_3_xi(P_A, element = 'stress')

		# Return results
		results = t, self.F(t), P_A_tot, total_rate_break, total_rate_reform, beta_sigma_over_n
		return results, P_A_end

	############################################################################################################################
	# Function for (one with k0 AND ignores yield)
	############################################################################################################################

	def compute_results_quad_specialized_ignore_yield(self, t_span):

		# Enumerate the time and relative deformation components
		t, F_zz_rel, F_rr_rel = self.enumerate_t_and_F_rel(t_span)

		# Function for the relatively-deformed coordinates
		def ell_rel(z, r, index_t, index_tau):
			return np.sqrt(z*z*F_zz_rel[index_t, index_tau]**2 + r*r*F_rr_rel[index_t, index_tau]**2)

		########################################################################################################################
		# Time-dependent quantities for spatial integrals
		########################################################################################################################

		# Function for the relatively-deformed equilibrium distribution
		def P_A_eq_rel(z, r, index_t, index_tau):
			return np.nan_to_num(self.initialized_single_chain_model.P_A_eq(ell_rel(z, r, index_t, index_tau)/self.J_sw**(1/3),\
				normalization = self.P_A_eq_normalization*self.J_sw), nan = 0)

		# Function for the reaction propagator
		def Xi(t, tau):
			return np.exp(-self.k_0*(t - tau))

		# Function for relative time derivative of the reaction propagator
		def d_Xi_d_tau(t, tau):
			return Xi(t, tau)*self.k_0

		########################################################################################################################
		# Compute stress only at unique deformations to interpolate from
		########################################################################################################################

		# Choose deformation for stress response based on deformation mode
		if self.deformation_type == 'uniaxial':
			F_use = F_zz_rel

		# Store a limited number of unique deformations over the history
		F_unique = np.unique(F_use)
		indices = np.unique(np.round(np.linspace(0, len(F_unique) - 1, num_interp_quad_specialized_ignore_yield)).astype(int))
		F_store = F_unique[indices]

		# Compute the stress at these deformations to interpolate from
		beta_sigma_h0_over_n = np.zeros(len(F_store))
		for index in range(len(F_store)):
			index_t, index_tau = np.argwhere(F_store[index] == F_use)[0]
			beta_sigma_h0_over_n[index] = \
				self.integral_quad_d_3_xi(lambda z,r: P_A_eq_rel(z, r, index_t, index_tau), element = 'stress')

		########################################################################################################################
		# Compute and return the results
		########################################################################################################################

		# Function to interpolate from computed stress response
		interp_sigma_fun = interp1d(F_store, beta_sigma_h0_over_n, kind = 'cubic', bounds_error = True)

		# Homogeneous solution for the nondimensional stress
		beta_sigma_h_over_n = self.g(t, 0)*Xi(t, 0)*interp_sigma_fun(F_zz_rel[:,0])

		# Particular solution for the nondimensional stress
		beta_sigma_p_over_n = np.zeros(self.num_time)
		for index_t in range(self.num_time):
			integrand = interp_sigma_fun(F_zz_rel[index_t,:index_t + 1])*( \
				Xi(t[index_t], t[:index_t + 1])*self.d_g_d_tau(t[index_t], t[:index_t + 1]) \
				+ d_Xi_d_tau(t[index_t], t[:index_t + 1])*self.g(t[index_t], t[:index_t + 1]))
			beta_sigma_p_over_n[index_t] = self.integral_d_tau(np.append(integrand, np.zeros(self.num_time - 1 - index_t)))

		# Return results
		beta_sigma_over_n = beta_sigma_h_over_n + beta_sigma_p_over_n
		P_A_tot = self.initialized_single_chain_model.P_A_tot_eq*np.ones(self.num_time)
		total_rate_reform = self.K_hat*(1 - P_A_tot)
		total_rate_break = -total_rate_reform
		results = t, self.F(t), P_A_tot, beta_sigma_h_over_n, beta_sigma_p_over_n, beta_sigma_over_n
		return results, np.nan

	############################################################################################################################
	# Function to adjust discretization for Romberg integration 
	############################################################################################################################

	def adjust_for_romb(self, num_discretization, decrease = False):
		if ((np.log(num_discretization - 1)/np.log(2)).is_integer()):
			return int(round(num_discretization))
		else:
			n = 0
			dos_check = 3
			while dos_check >= 2:
				n += 1
				dos_check = (num_discretization - 1)**(1/n)
			if decrease is True and dos_check < 2:
				return int(1 + 2**(n - 1))
			else:
				return int(1 + 2**n)

	############################################################################################################################
	# Function for integration over the spatial grid
	############################################################################################################################
	
	def integral_grid_d_3_xi(self, FUN, element = None):
		if element is None:
			element = self.R
		elif element == 'stress':
			element = self.ELEMENT_stress
		if FUN.ndim == 2:
			return 4*np.pi*romb(romb(FUN*element, dx = self.dr, axis = 0), dx = self.dz, axis = 0)
		elif FUN.ndim == 3:
			return 4*np.pi*romb(romb(FUN*element[:,:,None], dx = self.dr, axis = 0), dx = self.dz, axis = 0)
		elif FUN.ndim == 4:
			return 4*np.pi*romb(romb(FUN*element[:,:,None,None], dx = self.dr, axis = 0), dx = self.dz, axis = 0)

	############################################################################################################################
	# Function for integration over continuous space
	############################################################################################################################
	
	def integral_quad_d_3_xi(self, fun, element = None):
		if element is None:
			def integrand(z, r):
				return fun(z, r)*r
		elif element == 'stress':
			def integrand(z, r):
				return fun(z, r)*self.element_stress(z, r, self.initialized_single_chain_model)
		if self.initialized_single_chain_model.gamma_c is None:
			lim = self.initialized_single_chain_model.gamma_TS
		else:
			lim = self.initialized_single_chain_model.gamma_c
		return 4*np.pi*dblquad(integrand, 0, lim, lambda r: 0, lambda r: np.sqrt(lim**2 - r**2), \
			epsabs = dblquad_epsabs, epsrel = dblquad_epsrel)[0]

	############################################################################################################################
	# Function for integration element specialized for stress calculation
	############################################################################################################################

	def element_stress(self, z, r, single_chain_model):
		ell = np.sqrt(z*z + r*r)
		eta = np.nan_to_num(single_chain_model.eta(ell), nan = 0)
		if isinstance(z, np.ndarray):
			eta_over_ell = np.zeros(z.shape)
			eta_over_ell[ell != 0] = eta[ell != 0]/ell[ell != 0]
		else:
			if ell == 0:
				eta_over_ell = 0
			else:
				eta_over_ell = eta/ell
		C = (single_chain_model.N_b + single_chain_model.varsigma*single_chain_model.N_b_H)/self.J_sw
		if self.deformation_type == 'uniaxial':
			return C*eta_over_ell*r*(z*z - r*r/2)
		elif self.deformation_type == 'equibiaxial' or deformation_type == 'simple_shear':
			return C*eta_over_ell*r*(r*r/2 - z*z)

	############################################################################################################################
	# Function to interpolate from a stored 2D function on the spatial grid
	############################################################################################################################
	
	def interp_fun_2D(self, z_query, r_query, FUN):
		return interp2d(self.z, self.r, FUN, kind = self.interp_kind_2D)(z_query, r_query)

	############################################################################################################################
	# Functions for integration in time
	############################################################################################################################
	
	def integral_ds(self, FUN):
		return simpson(FUN, dx = self.timestep, axis = -1, even = 'last')

	def integral_d_tau(self, FUN):
		return romb(FUN, dx = self.timestep, axis = -1)

	############################################################################################################################
	# Function to enumerate the time and relative deformation
	############################################################################################################################

	def enumerate_t_and_F_rel(self, t_span):

		# Enumerate the time
		t = np.linspace(t_span[0], t_span[-1], self.num_time)

		# Relative deformation gradient components
		if self.deformation_type == 'uniaxial':
			F_zz_rel = np.tensordot(1/self.F(t), self.F(t), axes = 0)
			F_rr_rel = np.tensordot(np.sqrt(self.F(t)), 1/np.sqrt(self.F(t)), axes = 0)
		elif self.deformation_type == 'equibiaxial':
			F_zz_rel = np.tensordot(self.F(t), 1/self.F(t), axes = 0)
			F_rr_rel = np.tensordot(1/self.F(t), self.F(t), axes = 0)
		elif self.deformation_type == 'simple_shear':
			pass

		return t, F_zz_rel, F_rr_rel

	############################################################################################################################
	# Function to solve the Volterra integral equation
	############################################################################################################################

	def solve_Volterra(self, K, b, total_time):
		M = 0
		rho = b
		residual_bound_rho = 1
		while residual_bound_rho > tol_residual_rho:
			M += 1
			rho = b - self.integral_d_tau(K*rho)
			residual_bound_rho = \
				(self.K_hat*total_time)**(M + 1)/self.initialized_single_chain_model.P_B_tot_eq/np.math.factorial(M + 1)
		return rho

################################################################################################################################
# Checkpoint creation class
################################################################################################################################

class checkpoint:

	# Initialization also clears any previous checkpoint
	def __init__(self, checkpoint_directory):
		self.checkpoint_directory_and_file = checkpoint_directory + 'checkpoint.csv'
		open(self.checkpoint_directory_and_file, 'w').close()

	# Function to create checkpoints
	def create(self, t_end, P_A_end):
		f = open(self.checkpoint_directory_and_file, 'a')
		f.write("%.8e" % t_end)
		for i in range(len(P_A_end[:,0])):
			f.write("\n")
			for j in range(len(P_A_end[0,:])):
				f.write("%.8e\t" % P_A_end[i,j])
		f.close()

	# Function to read checkpoints
	def read(self, existing_checkpoint_directory_and_file):
		pass

################################################################################################################################
# Results writing class
################################################################################################################################

class results_csv:

	# Initialization also clears any previous results
	def __init__(self, csv_directory):
		self.csv_directory_and_file = csv_directory + 'results.csv'
		open(self.csv_directory_and_file, "w").close()

	# Function to append results
	def append(self, index_chunk, results):
		f = open(self.csv_directory_and_file, "a")
		for index_t in range(len(results[0])):
			for index_results in range(len(results)):
				f.write("%.8e\t" % results[index_results][index_t])
			f.write("\n")
		f.close()