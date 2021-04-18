################################################################################################################################
# General setup
################################################################################################################################

# Import libraries
import numpy as np
import matplotlib.pyplot as plt

################################################################################################################################
# Plotting class
################################################################################################################################
class plotter:

	def __init__(self, plot_directory = './'):
		self.plot_directory = plot_directory

	# For general plotting
	def save_current_figure(self, xlabel, ylabel, name):
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		plt.tight_layout()
		plt.show()
		fig = plt.gcf()
		plt.savefig(self.plot_directory + name)
		plt.close()

	############################################################################################################################
	# Function for plotting deformation
	############################################################################################################################

	def plot_deformation(self, F, total_time_in_seconds, csv_directory = None):

		# Enumerate time
		t_temp = np.linspace(0, total_time_in_seconds, int(1e5))

		# Plot the deformation
		fig = plt.figure()
		plt.plot(t_temp, F(t_temp))
		self.save_current_figure('$t$ [seconds]', '$F(t)$', 'F(t).png')

		# Save values to .csv if opted
		if csv_directory is not None:
			self.save_csv(csv_directory + 'deformation.csv', t_temp, F(t_temp))

	############################################################################################################################
	# Function for plotting single-chain functions
	############################################################################################################################

	def plot_single_chain(self, single_chain_model, J_sw = None, csv_directory = None):

		# Enumerate stretch
		if np.isinf(single_chain_model.gamma_TS):
			gamma_plot = np.linspace(0, 1.5, 1000)
		else:
			gamma_plot = np.linspace(0, single_chain_model.gamma_TS, 1000)

		# Plot the nondimensional mechanical response
		fig = plt.figure()
		plt.plot(gamma_plot, single_chain_model.eta(gamma_plot))
		self.save_current_figure('$\gamma$', '$\eta(\gamma)$', 'eta.png')

		# Plot the nondimensional equilibrium distribution
		fig = plt.figure()
		plt.plot(gamma_plot, single_chain_model.P_A_eq(gamma_plot), label = 'original')
		if J_sw is not None:
			plt.plot(gamma_plot, single_chain_model.P_A_eq(gamma_plot/J_sw**(1/3))/J_sw, label = 'swollen')
			plt.legend()
		self.save_current_figure('$\gamma$', r'$\mathscr{P}_A^{eq}(\gamma)$', 'P_A_eq.png')

		# Plot the nondimensional equilibrium radial distribution function
		fig = plt.figure()
		plt.plot(gamma_plot, single_chain_model.g_A_eq(gamma_plot), label = 'original')
		if J_sw is not None:
			plt.plot(gamma_plot, single_chain_model.g_A_eq(gamma_plot/J_sw**(1/3))/J_sw**(1/3), label = 'swollen')
			plt.legend()
		self.save_current_figure('$\gamma$', r'$\mathscr{g}_A^{eq}(\gamma)$', 'g_A_eq.png')

		# Plot the reaction rate coefficient function
		if np.all(np.isclose(single_chain_model.k(gamma_plot[:-2]), single_chain_model.k_0)) == False:
			fig = plt.figure()
			plt.plot(gamma_plot, single_chain_model.k(gamma_plot))
			plt.yscale('log')
			self.save_current_figure('$\gamma$', '$k(\gamma)$ [1/seconds]', 'k.png')

		# Save values to .csv if opted
		if csv_directory is not None:
			y_tuple = single_chain_model.eta(gamma_plot), single_chain_model.P_A_eq(gamma_plot), \
				single_chain_model.g_A_eq(gamma_plot), single_chain_model.k(gamma_plot)
			self.save_csv(csv_directory + 'single_chain.csv', gamma_plot, y_tuple)

	############################################################################################################################
	# Function to save plot values to .csv
	############################################################################################################################

	def save_csv(self, csv_file, x, y_tuple):

		# Initialization also clears any previous results
		open(csv_file, "w").close()

		# Save values
		f = open(csv_file, "a")
		for index in range(len(x)):
			f.write("%.8e\t" % x[index])
			if isinstance(y_tuple, tuple):
				for y in y_tuple:
					f.write("%.8e\t" % y[index])
			else:
				f.write("%.8e\t" % y_tuple[index])
			f.write("\n")
		f.close()

	############################################################################################################################
	# Function for plotting results
	############################################################################################################################

	def plot_results(self, deformation_object, results, use_nominal = False, n_over_beta = None, \
		data_F_stress = None, F_stress_1 = None, F_stress_2 = None, F_stress_3 = None):

		# Retrieve results
		t = results[0]
		F = results[1]
		beta_sigma_over_n = results[5]

		# Simpler plotting if given None as the deformation_object
		if deformation_object is not None:

			# Plot the component of the stress solution if ignoring chain breaking
			if deformation_object.ignore_yield is True:

				# Plot the homogeneous and particular solutions for the stress
				beta_sigma_h_over_n = results[3]
				beta_sigma_p_over_n = results[4]
				fig = plt.figure()
				plt.plot(F, beta_sigma_h_over_n, label = 'homogeneous')
				plt.plot(F, beta_sigma_p_over_n, label = 'particular')
				plt.plot(F, beta_sigma_over_n, label = 'total')
				plt.legend()
				self.save_current_figure('$F(t)$', r'$\beta\sigma(t)/n$', 'sigma_hpt(F).png')

				fig = plt.figure()
				plt.plot(t, beta_sigma_h_over_n, label = 'homogeneous')
				plt.plot(t, beta_sigma_p_over_n, label = 'particular')
				plt.plot(t, beta_sigma_over_n, label = 'total')
				plt.legend()
				self.save_current_figure('$t$ [seconds]', r'$\beta\sigma(t)/n$', 'sigma_hpt(t).png')

			# Typical plotting
			else:

				# Plot total probability of intact chains
				P_A_tot = results[2]
				fig = plt.figure()
				plt.plot(F, deformation_object.initialized_single_chain_model.P_A_tot_eq + 0*t, linestyle = 'dashed')
				plt.plot(F, P_A_tot)
				self.save_current_figure('$F(t)$', '$P_\mathrm{A}^\mathrm{tot}(t)$', 'P_A_tot(F).png')

				fig = plt.figure()
				plt.plot(t, deformation_object.initialized_single_chain_model.P_A_tot_eq + 0*t, linestyle = 'dashed')
				plt.plot(t, P_A_tot)
				self.save_current_figure('$t$ [seconds]', '$P_\mathrm{A}^\mathrm{tot}(t)$', 'P_A_tot(t).png')

				# Plot total rate of breaking and reforming chains
				total_rate_break = results[3]
				total_rate_reform = results[4]
				fig = plt.figure()
				plt.plot(F, total_rate_break, label = 'breaking rate')
				plt.plot(F, total_rate_reform, label = 'reforming rate')
				plt.plot(F, total_rate_reform + total_rate_break, label = 'net rate')
				plt.legend()
				self.save_current_figure('$F(t)$', r'$\frac{d}{dt}\,P_\mathrm{A}^\mathrm{tot}(t)$', 'd_P_A_tot_dt(F).png')

				fig = plt.figure()
				plt.plot(t, total_rate_break, label = 'breaking rate')
				plt.plot(t, total_rate_reform, label = 'reforming rate')
				plt.plot(t, total_rate_reform + total_rate_break, label = 'net rate')
				plt.legend()
				self.save_current_figure('$t$ [seconds]', \
					r'$\frac{d}{dt}\,P_\mathrm{A}^\mathrm{tot}(t)$', 'd_P_A_tot_dt(t).png')

		# Plot the stress
		fig = plt.figure()
		if n_over_beta is None:
			n_over_beta = 1
			y_label = r'$\beta\sigma(t)/n$'
		else:
			y_label = r'$\sigma(t)$'

		if use_nominal is True:
			y_label = 'nominal ' + y_label
			stress = n_over_beta*beta_sigma_over_n/F
		else:
			stress = n_over_beta*beta_sigma_over_n

		# Plot the stress data
		if data_F_stress is not None:
			plt.plot(data_F_stress[0], data_F_stress[1], 'o')

		# Plot any extra stress
		if F_stress_1 is not None:
			if use_nominal is True:
				sigma = F_stress_1[1]/F_stress_1[0]
			else:
				sigma = F_stress_1[1]
			plt.plot(F_stress_1[0], sigma, '--')
		if F_stress_2 is not None:
			if use_nominal is True:
				sigma = F_stress_2[1]/F_stress_2[0]
			else:
				sigma = F_stress_2[1]
			plt.plot(F_stress_2[0], sigma, '--')
		if F_stress_3 is not None:
			if use_nominal is True:
				sigma = F_stress_3[1]/F_stress_3[0]
			else:
				sigma = F_stress_3[1]
			plt.plot(F_stress_3[0], sigma, '--')

		# Plot the stress as a function of deformation
		plt.plot(F, stress)
		self.save_current_figure('$F(t)$', y_label, 'sigma(F).png')

		# Plot the stress as a function of time
		if deformation_object is not None:
			fig = plt.figure()
			plt.plot(t, stress)
			self.save_current_figure('$t$ [seconds]', y_label, 'sigma(t).png')
