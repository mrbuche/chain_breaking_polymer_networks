################################################################################################################################
# General setup
################################################################################################################################

# Import libraries
import sys
import copy
import numpy as np
from scipy.integrate import romb
from scipy.optimize import root_scalar

################################################################################################################################
# Deformation application class
################################################################################################################################

class deform_network_no_symmetry:

	############################################################################################################################
	# Initialization
	############################################################################################################################

	def __init__(self, single_chain_model, J_sw = 1, num_grid_suggestion = 129):

		# Inherit single_chain_model
		self.scm = single_chain_model

		# Volumetric swelling ratio
		self.J_sw = J_sw

		# Spatial grid octant (F diagonal, i.e. no rotation)
		self.num_grid = self.adjust_for_romb(num_grid_suggestion)
		self.x = np.linspace(0, self.scm.gamma_TS, self.num_grid)
		self.y = self.x
		self.z = self.x
		self.dx = self.x[1] - self.x[0]
		self.dy = self.dx
		self.dz = self.dx
		self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z)
		self.ELL = np.sqrt(self.X*self.X + self.Y*self.Y + self.Z*self.Z)

		# Integration element specialized for stress calculation
		self.ELEMENT_stress = self.element_stress(self.X, self.Y, self.Z, self.scm)

		# Adjust normalization of P_A_eq on the grid
		P_A_eq_ELL_non_normalized = np.nan_to_num(self.scm.P_A_eq(self.ELL), nan = 0)
		self.P_A_eq_normalization = self.integral_grid_d_3_xi(P_A_eq_ELL_non_normalized)/self.scm.P_A_tot_eq

	############################################################################################################################
	# Function to compute the results given a deformation and yield function
	############################################################################################################################

	def compute_results(self, index_t, F, update_yield_function = False):

		# Enumerate the relatively-deformed initial distribution
		F_inv = np.linalg.inv(F)
		P_A_0_rel_t = self.scm.P_A_eq(self.ELL_rel(F_inv)/self.J_sw**(1/3))/self.P_A_eq_normalization/self.J_sw
		P_A_0_rel_t[np.isnan(P_A_0_rel_t)] = 0

		# Enumerate the yield function
		ELL_rel = np.zeros((self.num_grid, self.num_grid, self.num_grid, index_t + 1))
		ELL_rel[:, :, :, -1] = self.ELL
		for index_index_t in range(index_t):
			ELL_rel[:, :, :, index_index_t] = self.ELL_rel(self.F[:, :, index_index_t]*F_inv)
		yield_function = np.trapz(ELL_rel > self.scm.gamma_c, x = self.t[:index_t + 1], axis = -1) == 0

		# Homogeneous solution for P_A
		P_A = P_A_0_rel_t*yield_function

		# Total probability that a chain is intact
		P_A_tot = self.integral_grid_d_3_xi(P_A)

		# Nondimensional stress corresponding to the applied deformation
		beta_sigma_over_n_plus_p = np.diag([ \
				self.integral_grid_d_3_xi(P_A, element = 'stress', component = [1, 1]), \
				self.integral_grid_d_3_xi(P_A, element = 'stress', component = [2, 2]), \
				self.integral_grid_d_3_xi(P_A, element = 'stress', component = [3, 3])  \
			])

		# Return results
		return P_A_tot, beta_sigma_over_n_plus_p

	############################################################################################################################
	# Function to solve for results over the applied deformation history
	############################################################################################################################

	def solve(self, total_time_in_seconds, F_applied, timestep_in_seconds = None, n_over_beta = None):
		
		# Enumerate time discretization
		if timestep_in_seconds is None:
			self.timestep_in_seconds = total_time_in_seconds/25
		else:
			self.timestep_in_seconds = timestep_in_seconds
		num_time = 1 + np.ceil(total_time_in_seconds/self.timestep_in_seconds).astype(int)
		self.t = np.linspace(0, total_time_in_seconds, num_time)

		# Verify maximum array size is acceptable (memory) before starting
		A = np.zeros((self.num_grid, self.num_grid, self.num_grid, num_time))
		del A

		# Allocate results
		self.F = np.zeros((3, 3, num_time))
		self.P_A_tot = np.zeros(num_time)
		self.p = np.zeros(num_time)
		self.beta_sigma_over_n = np.zeros((3, 3, num_time))
		self.beta_sigma_nominal_over_n = np.zeros((3, 3, num_time))

		# Loop through time
		for index_t in range(num_time):

			# When appling F_11(t)
			if F_applied.component(self.t[index_t]) == 11:
				F_11 = F_applied(self.t[index_t])
				def fun(F_33):
					F_22 = 1/(F_11*F_33)
					beta_sigma_over_n_plus_p = self.compute_results(index_t, np.diag([F_11, F_22, F_33]))[-1]
					return beta_sigma_over_n_plus_p[1, 1] - beta_sigma_over_n_plus_p[2, 2]

			# When applying F_22(t)
			elif F_applied.component(self.t[index_t]) == 22:
				F_22 = F_applied(self.t[index_t])
				def fun(F_33):
					F_11 = 1/(F_22*F_33)
					beta_sigma_over_n_plus_p = self.compute_results(index_t, np.diag([F_11, F_22, F_33]))[-1]
					return beta_sigma_over_n_plus_p[0, 0] - beta_sigma_over_n_plus_p[2, 2]

			# Solve for F_33(t) using incompressibility and the traction-free boundary conditions
			if index_t == 0:
				F_33 = 1
			F_33 = root_scalar(fun, x0 = F_33*0.95, x1 = F_33/0.95).root

			# Compute F_22(t) or F_11(t) and thus retrieve F(t)
			if F_applied.component(self.t[index_t]) == 11:
				F_22 = 1/(F_11*F_33)
			elif F_applied.component(self.t[index_t]) == 22:
				F_11 = 1/(F_22*F_33)
			F = np.diag([F_11, F_22, F_33])

			# Compute and store the results, update yield function
			P_A_tot, beta_sigma_over_n_plus_p = self.compute_results(index_t, F, update_yield_function = True)
			self.F[:,:, index_t] = F
			self.P_A_tot[index_t] = P_A_tot
			self.p[index_t] = beta_sigma_over_n_plus_p[2, 2]
			self.beta_sigma_over_n[:, :, index_t] = beta_sigma_over_n_plus_p - self.p[index_t]*np.diag([1, 1, 1])
			self.beta_sigma_nominal_over_n[:, :, index_t] = np.linalg.det(self.F[:,:, index_t]) \
				*np.linalg.inv(self.F[:,:, index_t]).dot(self.beta_sigma_over_n[:, :, index_t])

		# Compute and store other results
		if n_over_beta is not None:
			self.n_over_beta = n_over_beta
			self.sigma = n_over_beta*self.beta_sigma_over_n
			self.sigma_nominal = n_over_beta*self.beta_sigma_nominal_over_n

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
	
	def integral_grid_d_3_xi(self, FUN, element = 1, component = None):
		if element == 'stress':
			element = self.ELEMENT_stress[:, :, :, component[0] - 1, component[1] - 1]
		if FUN.ndim == 3:
			return 8*romb(romb(romb(FUN*element, \
				dx = self.dx, axis = 0), dx = self.dy, axis = 0), dx = self.dz, axis = 0)
		elif FUN.ndim == 4:
			return 8*romb(romb(romb(FUN*element[:,:,:,None], \
				dx = self.dx, axis = 0), dx = self.dy, axis = 0), dx = self.dz, axis = 0)
		elif FUN.ndim == 5:
			return 8*romb(romb(romb(FUN*element[:,:,:,None,None], \
				dx = self.dx, axis = 0), dx = self.dy, axis = 0), dx = self.dz, axis = 0)

	############################################################################################################################
	# Function for integration element specialized for stress calculation
	############################################################################################################################

	def element_stress(self, x, y, z, single_chain_model):
		ell = np.sqrt(x*x + y*y + z*z)
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
		element_out = np.zeros((ell.shape[0], ell.shape[1], ell.shape[2], 3, 3))
		element_out[:, :, :, 0, 0] = C*eta_over_ell*x*x
		element_out[:, :, :, 0, 1] = C*eta_over_ell*x*y
		element_out[:, :, :, 0, 2] = C*eta_over_ell*x*z
		element_out[:, :, :, 1, 0] = C*eta_over_ell*y*x
		element_out[:, :, :, 1, 1] = C*eta_over_ell*y*y
		element_out[:, :, :, 1, 2] = C*eta_over_ell*y*z
		element_out[:, :, :, 2, 0] = C*eta_over_ell*z*x
		element_out[:, :, :, 2, 1] = C*eta_over_ell*z*y
		element_out[:, :, :, 2, 2] = C*eta_over_ell*z*z
		return element_out

	############################################################################################################################
	# Function to return relatively-deformed coordinates
	############################################################################################################################

	def ELL_rel(self, F_rel):
		X_rel = F_rel[0, 0]*self.X + F_rel[0, 1]*self.Y + F_rel[0, 2]*self.Z
		Y_rel = F_rel[1, 0]*self.X + F_rel[1, 1]*self.Y + F_rel[1, 2]*self.Z
		Z_rel = F_rel[2, 0]*self.X + F_rel[2, 1]*self.Y + F_rel[2, 2]*self.Z
		return np.sqrt(X_rel*X_rel + Y_rel*Y_rel + Z_rel*Z_rel)