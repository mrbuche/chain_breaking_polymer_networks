################################################################################################################################
# General setup
################################################################################################################################

# Import libraries
import sys
import numpy as np
from scipy.integrate import quad

################################################################################################################################
# Relaxation function from Long et al. 2014
################################################################################################################################

class Long_et_al_2014:

	# For more information, see:
	# 	Time Dependent Behavior of a Dual Cross-Link Self-Healing Gel: Theory and Experiments
	# 	Rong Long, Koichi Mayumi, Costantino Creton, Tetsuharu Narita, and Chung-Yuen Hui
	# 	Macromolecules 2014, 47, 7243−7250
	# 	doi.org/10.1021/ma501290h
	# See also:
	#	Mechanics of a Dual Cross-Link Gel with Dynamic Bonds: Steady State Kinetics and Large Deformation Effects
	# 	Jingyi Guo, Rong Long, Koichi Mayumi, and Chung-Yuen Hui
	# 	Macromolecules 2016, 49, 3497−3507
	#	doi.org/10.1021/acs.macromol.6b00421
	# See also:
	# 	(title)
	# 	Michael R. Buche and Meredith N. Silberstein
	#	(journal)
	# 	doi.org/

	def __init__(self, **kwargs):

		# Default parameter values
		self.alpha = None
		self.t_R = None
		self.x_p = 0

		# Retrieve specified parameters
		for key, value in kwargs.items():
			if key == 'alpha':
				self.alpha = value
			elif key == 't_R':
				self.t_R = value
			elif key == 'x_p':
				self.x_p = value

		# Check parameter specifications
		if self.alpha is None:
			sys.exit("Error: Need to specify alpha for Long_et_al_2014().")
		elif self.alpha < 1:
			sys.exit("Error: Need to specify alpha > 1 for Long_et_al_2014().")
		elif self.t_R is None:
			sys.exit("Error: Need to specify t_R for Long_et_al_2014().")

		# Smallest timescale associated with this relaxation function
		self.timescale = self.t_R

	# Relaxation function
	def g(self, t, tau):
		return self.x_p + (1 - self.x_p)*(1 + (self.alpha - 1)*(t - tau)/self.t_R)**(1/(1 - self.alpha))

	# Relative time derivative of the relaxation function
	def d_g_d_tau(self, t, tau):
		return (1 - self.x_p)/self.t_R*(1 + (self.alpha - 1)*(t - tau)/self.t_R)**(self.alpha/(1 - self.alpha))

	# Storage function
	def g_p(self, omega):
		g_p_fun = lambda omega: \
			self.x_p + omega*quad(lambda s: (self.g(s, 0) - self.x_p)*np.sin(omega*s), 0, np.inf, full_output = 1)[0]
		return np.vectorize(g_p_fun)(omega)

	# Loss function
	def g_pp(self, omega):
		g_pp_fun = lambda omega: \
			omega*quad(lambda s: (self.g(s, 0) - self.x_p)*np.cos(omega*s), 0, np.inf, full_output = 1)[0]
		return np.vectorize(g_pp_fun)(omega)

################################################################################################################################
# Relaxation function for the Rouse model
################################################################################################################################

class Rouse:

	# For more information, see:
	# 	A Theory of the Linear Viscoelastic Properties of Dilute Solutions of Coiling Polymers 
	#	Prince E Rouse
	#	Journal of Chemical Physics, 1953, 21, 1272–1280
	#	doi.org/10.1063/1.1699180
	# See also:
	#	Extensions of the Rouse Theory of Viscoelastic Properties to Undiluted Linear Polymers
	#	John D. Ferry, Robert F. Landel, and Malcolm L. Williams
	#	Journal of Applied Physics, 1955, 26, 359
	#	doi.org/10.1063/1.1721997
	# See also:
	# 	Bridging experiments and theory: 
	#		isolating the effects of metal–ligand interactions on viscoelasticity of reversible polymer networks
	# 	Xinyue Zhang, Yuval Vidavsky, Sinai Aharonovich, Steven J. Yang, Michael R. Buche, 
	#		Charles E. Diesendruck and Meredith N. Silberstein
	#	Soft Matter, 2020, 16, 8591-8601
	#	doi.org/10.1039/D0SM01115K

	def __init__(self, **kwargs):

		# Default parameter values
		self.N_b = None
		self.t_0 = None

		# Retrieve specified parameters
		for key, value in kwargs.items():
			if key == 'N_b':
				self.N_b = value
			elif key == 't_0':
				self.t_0 = value

		# Check parameter specifications
		if self.N_b is None:
			sys.exit("Error: Need to specify N_b for Rouse().")
		elif self.t_0 is None:
			sys.exit("Error: Need to specify t_0 for Rouse().")

		# Smallest timescale associated with this relaxation function
		self.timescale = self.t_0

	# Relaxation function
	def g(self, t, tau):
		g_out = np.zeros(omega.shape)
		for p in range(1, int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			g_out = g_out + np.exp(-(t - tau)/t_p)/self.N_b
		return g_out

	# Relative time derivative of the relaxation function
	def d_g_d_tau(self, t, tau):
		d_g_d_tau_out = np.zeros(omega.shape)
		for p in range(1, int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			d_g_d_tau_out = d_g_d_tau_out + np.exp(-(t - tau)/t_p)/self.N_b/t_p
		return d_g_d_tau_out

	# Storage function
	def g_p(self, omega):
		g_p_out = np.zeros(omega.shape)
		for p in range(1, int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			g_p_out = g_p_out + ((omega*t_p)**2/(1 + (omega*t_p)**2))/self.N_b
		return g_p_out

	# Loss function
	def g_pp(self, omega):
		g_pp_out = np.zeros(omega.shape)
		for p in range(1, int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			g_pp_out = g_pp_out + (omega*t_p/(1 + (omega*t_p)**2))/self.N_b
		return g_pp_out

################################################################################################################################
# 
################################################################################################################################

class sticky_Rouse:

	# For more information, see:
	# 	Ionomer dynamics and the sticky Rouse model 
	#	Quan Chen, Gregory J. Tudryn, and Ralph H. Colby
	#	Journal of Rheology 2013, 57, 1441
	#	doi.org/10.1122/1.4818868
	# See also:
	# 	Bridging experiments and theory: 
	#		isolating the effects of metal–ligand interactions on viscoelasticity of reversible polymer networks
	# 	Xinyue Zhang, Yuval Vidavsky, Sinai Aharonovich, Steven J. Yang, Michael R. Buche, 
	#		Charles E. Diesendruck and Meredith N. Silberstein
	#	Soft Matter, 2020, 16, 8591-8601
	#	doi.org/10.1039/D0SM01115K

	def __init__(self, **kwargs):

		# Default parameter values
		self.N_b = None
		self.N_x = None
		self.t_0 = None
		self.t_x = None
		self.beta_E_A = None
		self.G_x_over_G_0 = 0

		# Retrieve specified parameters
		for key, value in kwargs.items():
			if key == 'N_b':
				self.N_b = value
			elif key == 'N_x':
				self.N_x = value
			elif key == 't_0':
				self.t_0 = value
			elif key == 't_x':
				self.t_x = value
			elif key == 'beta_E_a':
				self.beta_E_a = value
				self.t_x = self.t_0*np.exp(self.beta_E_a)
			elif key == 'G_x_over_G_0':
				self.G_x_over_G_0 = value

		# Check parameter specifications
		if self.N_b is None:
			sys.exit("Error: Need to specify N_b for sticky_Rouse().")
		elif self.N_x is None:
			sys.exit("Error: Need to specify N_x for sticky_Rouse().")
		elif self.t_0 is None:
			sys.exit("Error: Need to specify t_0 for sticky_Rouse().")
		elif self.t_x is None:
			sys.exit("Error: Need to specify t_x for sticky_Rouse().")

		# Smallest timescale associated with this relaxation function
		self.timescale = np.min([self.t_0, self.t_x])

	# Relaxation function
	def g(self, t, tau):
		g_out = np.zeros(omega.shape)
		for p in range(int(self.N_x + 1), int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			g_out = g_out + np.exp(-(t - tau)/t_p)/self.N_b
		for p in range(1, int(self.N_x + 1)):
			t_p_x = self.t_x*(self.N_x/p)**2
			g_out = g_out + (1 + self.G_x_over_G_0)*np.exp(-(t - tau)/t_p_x)/self.N_b
		return g_out

	# Relative time derivative of the relaxation function
	def d_g_d_tau(self, t, tau):
		d_g_d_tau_out = np.zeros(omega.shape)
		for p in range(int(self.N_x + 1), int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			d_g_d_tau_out = d_g_d_tau_out + np.exp(-(t - tau)/t_p)/self.N_b/t_p
		for p in range(1, int(self.N_x + 1)):
			t_p_x = self.t_x*(self.N_x/p)**2
			d_g_d_tau_out = d_g_d_tau_out + (1 + self.G_x_over_G_0)*np.exp(-(t - tau)/t_p_x)/self.N_b/t_p_x
		return d_g_d_tau_out

	# Storage function
	def g_p(self, omega):
		g_p_out = np.zeros(omega.shape)
		for p in range(int(self.N_x + 1), int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			g_p_out = g_p_out + ((omega*t_p)**2/(1 + (omega*t_p)**2))/self.N_b
		for p in range(1, int(self.N_x + 1)):
			t_p_x = self.t_x*(self.N_x/p)**2
			g_p_out = g_p_out + (1 + self.G_x_over_G_0)*((omega*t_p_x)**2/(1 + (omega*t_p_x)**2))/self.N_b
		return g_p_out

	# Loss function
	def g_pp(self, omega):
		g_pp_out = np.zeros(omega.shape)
		for p in range(int(self.N_x + 1), int(self.N_b + 1)):
			t_p = self.t_0*(self.N_b/p)**2
			g_pp_out = g_pp_out + (omega*t_p/(1 + (omega*t_p)**2))/self.N_b
		for p in range(1, int(self.N_x + 1)):
			t_p_x = self.t_x*(self.N_x/p)**2
			g_pp_out = g_pp_out + (1 + self.G_x_over_G_0)*(omega*t_p_x/(1 + (omega*t_p_x)**2))/self.N_b
		return g_pp_out
