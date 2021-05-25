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
	log_of_x_over_sinh_x = np.zeros(x.shape)

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
		self.P_B_tot_eq = 0
		self.K_hat = k_0
		self.max_k_rev = k_0
		self.N_b_H = 0
		self.varsigma = 1
		self.gamma_TS = self.gamma_c
		self.k = lambda gamma_in: k_0 + 0*gamma_in

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
				eta_out[harmonic_region] = self.kappa*(gamma_in[harmonic_region] - 1)
				eta_out[~harmonic_region] = self.eta_interp_fun(gamma_in[~harmonic_region])
			else:
				if gamma_in > cutoff_stretch_for_harmonic_eta_EFJC:
					eta_out = self.kappa*(gamma_in - 1)
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

		# Retain each single-chain function
		self.eta = eta_fun
		self.P_A_eq = lambda gamma_in: P_A_eq_fun(gamma_in, normalization = P_A_eq_normalization)
		self.g_A_eq = lambda gamma_in: g_A_eq_fun(gamma_in, normalization = P_A_eq_normalization)