################################################################################################################################
# General setup
################################################################################################################################

# Import libraries
import sys
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

# Interpolation parameters
num_interp = int(3e3)
interp_kind_1D = 'cubic'

# Numerical tolerance parameters
cutoff_for_log_over_sinh = 3e1
cutoff_stretch_for_harmonic_eta_EFJC = 3
minimum_exponent = np.log(sys.float_info.min)/np.log(10)
maximum_exponent = np.log(sys.float_info.max)/np.log(10)
eta_small = 10**minimum_exponent

# Function to invert a function
def inv_fun_1D(x_query, fun, bounds = None):

	# Change method depending on whether bounds are involved
	if bounds is None:
		return minimize_scalar(lambda x: np.abs(fun(x) - x_query)).x
	else:
		return minimize_scalar(lambda x: np.abs(fun(x) - x_query), bounds = bounds, method = 'bounded').x

# Function to create interpolation function from stored function
def interp_fun_1D(x_store, y_store):
	return interp1d(x_store, y_store, kind = interp_kind_1D, bounds_error = False, fill_value = np.nan)

# Function to avoid overflow when computing ln(x/sinh(x))
def log_over_sinh(x):

	# Determine when argument is sufficiently large
	where_x_large = np.nan_to_num(x, nan = -1) > cutoff_for_log_over_sinh
	log_of_x_over_sinh_x = np.zeros(x.size)

	# Use asymptotic relation valid for sufficiently large arguments
	if where_x_large.any():
		log_of_x_over_sinh_x[where_x_large] = np.log(2*x[where_x_large]) - x[where_x_large]

	# Compute analytically otherwise, and zero where argument is zero
	where_x_zero = x == 0
	where_compute = ~(where_x_large + where_x_zero)
	if where_compute.any():
		log_of_x_over_sinh_x[where_compute] = np.log(x[where_compute]/np.sinh(x[where_compute]))
	return log_of_x_over_sinh_x

# Hyperbolic cotangent function
def coth_safe(eta):
	eta = np.where(eta == 0, eta_small, eta)
	return 1/np.tanh(eta)

# Langevin function
def Langevin(eta):
	eta = np.where(eta == 0, eta_small, eta)
	return 1/np.tanh(eta) - 1/eta

################################################################################################################################
# Ideal chain model
################################################################################################################################

class ideal:

	# For more information, see:
	#	Statistical mechanical constitutive theory of polymer networks: 
	#		The inextricable links between distribution, behavior, and ensemble
	# 	Michael R. Buche and Meredith N. Silberstein
	# 	Physical Review E, 2021, 102, 012501
	# 	doi.org/10.1103/PhysRevE.102.012501

	# Class initialization
	def __init__(self, **kwargs):

		# Default parameter values
		N_b = 88
		k_0 = np.exp(minimum_exponent)
		gamma_c = np.inf

		# Retrieve specified parameters
		for key, value in kwargs.items():
			if key == 'N_b':
				N_b = value
			elif key == 'k_0':
				if value > 0:
					k_0 = value
			elif key == 'gamma_c':
				gamma_c = value

		# Retain for certain purposes
		self.N_b = N_b
		self.k_0 = k_0
		self.gamma_c = gamma_c

		# Model-specific modifications
		self.P_A_tot_eq = 1
		self.gamma_TS = np.inf
		self.k = lambda gamma_in: k_0 + 0*gamma_in
		self.K_hat = k_0
		self.max_k_rev = k_0
		self.N_b_H = 0
		self.varsigma = 1

	# Nondimensional mechamical response of the chain
	def eta(self, gamma_in):
		return 3*gamma_in

	# Nondimensional equilibrium distribution function
	def P_A_eq(self, gamma_in, normalization = 1):
		return (gamma_in <= self.gamma_c)*(3*self.N_b/2/np.pi)**(3/2)*np.exp(-3/2*self.N_b*gamma_in**2)/normalization

	# Nondimensional equilibrium radial distribution function
	def g_A_eq(self, gamma_in, normalization = 1):
		return 4*np.pi*gamma_in**2*self.P_A_eq(gamma_in, normalization)

################################################################################################################################
# Extensible freely-joined chain model
################################################################################################################################

class EFJC:

	# For more information, see:
	# 	Analytical results of the extensible freely jointed chain model
	# 	Alessandro Fiasconaro and Fernando Falo
	#	Physica A 2019, 532, 121929
	# 	doi.org/10.1016/j.physa.2019.121929
	# See also:
	#	Statistical mechanical constitutive theory of polymer networks: 
	#		The inextricable links between distribution, behavior, and ensemble
	# 	Michael R. Buche and Meredith N. Silberstein
	# 	Physical Review E, 2021, 102, 012501
	# 	doi.org/10.1103/PhysRevE.102.012501

	# Class initialization
	def __init__(self, **kwargs):

		# Default parameter values
		N_b = None
		kappa = None
		k_0 = np.exp(minimum_exponent)
		gamma_c = np.inf

		# Retrieve specified parameters
		for key, value in kwargs.items():
			if key == 'N_b':
				N_b = value
			elif key == 'k_0':
				if value > 0:
					k_0 = value
			elif key == 'kappa':
				kappa = value
			elif key == 'gamma_c':
				gamma_c = value

		# Check parameter specifications
		if N_b is None:
			sys.exit('Error: Need to specify N_b.')
		elif kappa is None:
			sys.exit('Error: Need to specify kappa.')

		# Retain for certain purposes
		self.N_b = N_b
		self.k_0 = k_0
		self.kappa = kappa
		self.gamma_c = gamma_c

		# Model-specific modifications
		self.P_A_tot_eq = 1
		self.gamma_TS = np.inf
		self.k = lambda gamma_in: k_0 + 0*gamma_in
		self.K_hat = k_0
		self.max_k_rev = k_0
		self.N_b_H = 0
		self.varsigma = 1

		# Nondimensional mechanical response of the chain
		def gamma_fun(eta):
			coth = coth_safe(eta)
			L = Langevin(eta)
			return L + eta/kappa*(1 + (1 - L*coth)/(1 + eta/kappa*coth))

		# Compute and store the inverted nondimensional mechanical response to interpolate from
		self.gamma_store = np.linspace(0, cutoff_stretch_for_harmonic_eta_EFJC, num_interp)
		self.eta_store = np.zeros(self.gamma_store.size)
		for i in range(1, len(self.gamma_store)):
			self.eta_store[i] = inv_fun_1D(self.gamma_store[i], gamma_fun)

		# Function to interpolate from the inverted nondimensional mechamical response of the chain
		self.eta_interp_fun = interp_fun_1D(self.gamma_store, self.eta_store)
		def eta_fun(gamma_in):
			if isinstance(gamma_in, np.ndarray):
				eta_out = np.zeros(gamma_in.shape)
				harmonic_region = gamma_in > cutoff_stretch_for_harmonic_eta_EFJC
				eta_out[harmonic_region] = kappa*(gamma_in[harmonic_region] - 1)
				eta_out[~harmonic_region] = self.eta_interp_fun(gamma_in[~harmonic_region])
			else:
				if gamma_in > cutoff_stretch_for_harmonic_eta_EFJC:
					eta_out = kappa*(gamma_in - 1)
				else:
					eta_out = self.eta_interp_fun(gamma_in)
			return eta_out

		# Nondimensional equilibrium distribution function
		def P_A_eq_fun(gamma_in, normalization = 1):

			# Compute mechanical response
			eta = np.array(eta_fun(gamma_in))
			eta[eta == 0] = eta_small

			# Compute nondimensional Helmholtz free energy per link
			coth = coth_safe(eta)
			L = Langevin(eta)
			vartheta = eta*L + log_over_sinh(eta) - np.log(1 + eta/kappa*coth) \
				+ eta**2/kappa/2*(1/2 + (1 - L*coth)/(1 + eta/kappa*coth))

			# Compute P_A_eq below the yield surface
			return (gamma_in <= self.gamma_c)*np.exp(-N_b*vartheta)/normalization

		# Nondimensional equilibrium radial distribution function
		def g_A_eq_fun(gamma_in, normalization = 1):
			return 4*np.pi*gamma_in**2*P_A_eq_fun(gamma_in, normalization)

		# Normalize the equilibrium distribution
		P_A_eq_normalization = quad(g_A_eq_fun, 0, np.inf, full_output = 1)[0]/self.P_A_tot_eq

		# Compute and store the function to interpolate from
		self.P_A_eq_store = P_A_eq_fun(self.gamma_store, normalization = P_A_eq_normalization)
		self.P_A_eq_interp_fun = interp_fun_1D(self.gamma_store, self.P_A_eq_store)

	# Function to interpolate from the inverted nondimensional mechamical response of the chain
	def eta(self, gamma_in):
		if isinstance(gamma_in, np.ndarray):
			eta_out = np.zeros(gamma_in.shape)
			harmonic_region = gamma_in > cutoff_stretch_for_harmonic_eta_EFJC
			eta_out[harmonic_region] = self.kappa*(gamma_in[harmonic_region] - 1)
			eta_out[~harmonic_region] = self.eta_interp_fun(gamma_in[~harmonic_region])
		else:
			if gamma_in > cutoff_stretch_for_harmonic_eta_EFJC:
				eta_out = self.kappa*(gamma_in - 1)
			else:
				eta_out = self.eta_interp_fun(gamma_in)
		return eta_out

	# Function to interpolate from the stored nondimensional equilibrium distribution function
	def P_A_eq(self, gamma_in, normalization = 1):
		return (gamma_in <= self.gamma_c)*self.P_A_eq_interp_fun(gamma_in)/normalization

	# Function for the nondimensional equilibrium radial distribution function
	def g_A_eq(self, gamma_in, normalization = 1):
		return 4*np.pi*gamma_in**2*self.P_A_eq(gamma_in, normalization)

################################################################################################################################
# Morse potential-supplemented freely-joined chain model
################################################################################################################################

class Morse_FJC:

	# For more information, see:
	# 	(title)
	# 	Michael R. Buche and Meredith N. Silberstein
	#	(journal)
	# 	doi.org/

	# Class initialization
	def __init__(self, **kwargs):

		# Default parameter values
		N_b = None
		N_b_H = 0
		k_0 = np.exp(minimum_exponent)
		omega = None
		grumbo = None
		beta_u_b = None
		kappa = None
		kappa_H = None
		varsigma = 1
		beta_Delta_Psi_0 = 0
		gamma_c = None

		# Retrieve specified parameters
		for key, value in kwargs.items():
			if key == 'N_b':
				N_b = value
			elif key == 'N_b_H':
				N_b_H = value
			elif key == 'k_0':
				if value > 0:
					k_0 = value
			elif key == 'omega':
				omega = value
			elif key == 'kappa':
				kappa = value
			elif key == 'grumbo':
				grumbo = value
			elif key == 'beta_u_b':
				beta_u_b = value
			elif key == 'kappa_H':
				kappa_H = value
			elif key == 'varsigma':
				varsigma = value
			elif key == 'beta_Delta_Psi_0':
				beta_Delta_Psi_0 = value
			elif key == 'gamma_c':
				gamma_c = value

		# Check parameter specifications
		if N_b is None:
			sys.exit('Error: Need to specify N_b in Morse_FJC().')
		elif kappa is None:
			sys.exit('Error: Need to specify kappa in Morse_FJC().')
		elif N_b_H == 0:
			kappa_H = 1 # kappa_H just has to be nonzero if N_b_H = 0
		elif kappa_H is None:
			sys.exit('Error: Specify nonzero kappa_H when specifying nonzero N_b_H in Morse_FJC().')
		elif N_b == 0:
			sys.exit('Error: For N_b = 0, use the class EFJC() instead of Morse_FJC().')
		if beta_u_b is not None:
			if grumbo is not None:
				sys.exit('Error: Either gamma or beta_u_b need to be specified in Morse_FJC(), but not both.')
			else:
				grumbo = np.sqrt(kappa/2/beta_u_b)
		elif beta_u_b is None:
			if grumbo is None:
				sys.exit('Error: Need to specify beta_u_b in Morse_FJC().')

		# Retain for certain purposes
		self.N_b = N_b
		self.N_b_H = N_b_H
		self.varsigma = varsigma
		self.k_0 = k_0
		self.gamma_c = gamma_c

		# Automatic physical parameters
		eta_max = kappa/4/grumbo
		self.P_A_tot_eq = 1/(1 + N_b*np.exp(-beta_Delta_Psi_0))
		self.P_B_tot_eq = 1 - self.P_A_tot_eq

		# Check if Python thinks P_A_tot_eq or P_B_tot_eq = 0
		if self.P_B_tot_eq**2 == 0:
			sys.exit('Error: P_B_tot_eq is smaller than precision allows; need to decrease beta_Delta_Psi_0 in Morse_FJC().')
		elif self.P_A_tot_eq**2 == 0:
			sys.exit('Error: P_A_tot_eq is smaller than precision allows; need to increase beta_Delta_Psi_0 in Morse_FJC().')

		# Nondimensional incremental mechanical response of the Morse potential
		def Delta_lambda(eta):
			return np.log(2/(1 + np.sqrt(1 - eta/eta_max)))/grumbo

		# Nondimensional Morse potential
		def beta_u(eta):
			return kappa/2/grumbo**2*(1 - np.exp(-grumbo*Delta_lambda(eta) ))**2

		# Nondimensional mechanical response of the chain
		def gamma_fun(eta):
			r = varsigma*N_b_H/N_b
			eta_c = varsigma*eta
			return 1/(1 + r)*(Langevin(eta) + Delta_lambda(eta) + r*(Langevin(eta_c) + eta_c/kappa_H))

		# Limiting stretch for an intact chain
		self.gamma_TS = gamma_fun(eta_max)

		# Ensure that gamma_c is below the limiting stretch for an intact chain
		if gamma_c is not None:
			if gamma_c > self.gamma_TS:
				sys.exit('Error: Cannot have gamma_c > gamma_TS = ' + str(self.gamma_TS))

		# Initial nondimensional free energy barrier to transition state
		beta_Delta_Psi_TS_0 = kappa/8/grumbo**2 - 2*np.log(self.gamma_TS)

		# Ensure either attempt frequency or k(0) have been specified
		if omega is None:
			if k_0 is None:
				sys.exit('Error: Either omega or k_0 need to be specified in Morse_FJC().')
			else:
				omega = 2*np.pi/N_b*np.exp(beta_Delta_Psi_TS_0 + np.log(k_0))
		elif k_0 is None:
			k_0 = N_b/2/np.pi*np.exp(-beta_Delta_Psi_TS_0 + np.log(omega))
		else:
			sys.exit('Error: Either omega or k_0 need to be specified in Morse_FJC(), but not both.')

		# Compute and store the inverted nondimensional mechanical response to interpolate from
		self.gamma_store = np.linspace(0, self.gamma_TS, num_interp)
		self.eta_store = np.zeros(self.gamma_store.size)
		for i in range(len(self.gamma_store)):
			self.eta_store[i] = inv_fun_1D(self.gamma_store[i], gamma_fun, bounds = (0, eta_max))

		# Function to interpolate from the inverted nondimensional mechamical response of the chain
		self.eta = interp_fun_1D(self.gamma_store, self.eta_store)

		# Nondimensional equilibrium distribution function
		def P_A_eq_fun(gamma_in, normalization = 1):

			# Compute nondimensional Helmholtz free energy per link
			eta = self.eta(gamma_in)
			r_N = N_b_H/N_b
			eta_c = varsigma*eta
			vartheta = 1/(1 + r_N)*(eta/np.tanh(eta) - 1 + log_over_sinh(eta) + beta_u(eta) \
				+ r_N*(eta_c/np.tanh(eta_c) - 1 + log_over_sinh(eta_c) + eta_c**2/2/kappa_H))

			# Avoid overflow before computing and returning P_A_eq
			exponent = -np.log(normalization) - (N_b + N_b_H)*vartheta
			exponent[exponent > maximum_exponent] = maximum_exponent
			exponent[exponent < minimum_exponent] = minimum_exponent
			return np.exp(exponent)

		# Nondimensional equilibrium radial distribution function
		def g_A_eq_fun(gamma_in, normalization = 1):
			return 4*np.pi*gamma_in**2*P_A_eq_fun(gamma_in, normalization)

		# Normalize the equilibrium distribution
		P_A_eq_normalization = quad(g_A_eq_fun, 0, self.gamma_TS, full_output = 1)[0]/self.P_A_tot_eq

		# Compute and store the function to interpolate from
		self.P_A_eq_store = P_A_eq_fun(self.gamma_store, normalization = P_A_eq_normalization)
		self.P_A_eq_interp_fun = interp_fun_1D(self.gamma_store, self.P_A_eq_store)

		# Reaction rate coefficient function
		if self.gamma_c is None:
			def k_fun(gamma_in):

				# Compute nondimensional Helmholtz free energy barrier
				eta = self.eta(gamma_in)
				adj_beta_Delta_Psi_TS = log_over_sinh(self.gamma_TS*eta) - log_over_sinh(eta) - beta_u(eta) \
				    + self.gamma_TS*eta/np.tanh(self.gamma_TS*eta) - eta/np.tanh(eta)

				# Avoid overflow before computing and returning k
				exponent = np.log(k_0) - N_b*adj_beta_Delta_Psi_TS
				exponent[exponent > maximum_exponent] = maximum_exponent
				exponent[exponent < minimum_exponent] = minimum_exponent
				return np.exp(exponent)

			# Compute and store the function to interpolate from
			self.log_k_store = np.log(k_fun(self.gamma_store))
			self.log_k_interp_fun = interp_fun_1D(self.gamma_store, self.log_k_store)

		# Total reverse reaction rate coefficient
		if self.gamma_c is None:
			integrand_K_hat = lambda xi: k_fun(xi)*g_A_eq_fun(xi, normalization = P_A_eq_normalization)
			self.K_hat = quad(integrand_K_hat, 0, self.gamma_TS, full_output = 1)[0]/self.P_B_tot_eq
		else:
			self.K_hat = k_0*self.P_A_tot_eq/self.P_B_tot_eq

		# Maximum reverse reaction rate coefficient
		if self.gamma_c is None:
			fun = lambda xi: -k_fun(xi)*P_A_eq_fun(xi, normalization = P_A_eq_normalization)/self.P_B_tot_eq
			self.max_k_rev = -fun(minimize_scalar(fun, bounds = (0, self.gamma_TS), method = 'bounded').x)[0]
		else:
			self.max_k_rev = k_0*P_A_eq_fun(0)[0]/self.P_B_tot_eq

	# Function to interpolate from the stored nondimensional equilibrium distribution function
	def P_A_eq(self, gamma_in, normalization = 1):
		return self.P_A_eq_interp_fun(gamma_in)/normalization

	# Function for the nondimensional equilibrium radial distribution function
	def g_A_eq(self, gamma_in, normalization = 1):
		return 4*np.pi*gamma_in**2*self.P_A_eq(gamma_in, normalization)

	# Function to interpolate from the stored reaction rate coefficient function
	def k(self, gamma_in):
		if self.gamma_c is None:
			return np.exp(self.log_k_interp_fun(gamma_in))
		else:
			if isinstance(gamma_in, np.ndarray):
				k_out = self.k_0*np.ones(gamma_in.shape)
				k_out[gamma_in > self.gamma_c] = np.inf
			else:
				if gamma_in > self.gamma_c:
					k_out = np.inf
				else:
					k_out = self.k_0
			return k_out
