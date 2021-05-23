# Example 3

	# Import the library
	from chain_breaking_polymer_networks import *

	# Create the single-chain models
	x_p = 5.85e-2
	single_chain_model_1 = EFJC(N_b = 50, kappa = 40, k_0 = 0.12)
	single_chain_model_2 = EFJC(N_b = 50, kappa = 40)

	# Define the deformation
	strain_rate = 0.01
	maximum_strain = 1.55
	total_time_in_seconds = maximum_strain/strain_rate
	def F(t):
		return 1 + strain_rate*t

	# Apply the deformation and solve
	results_1 = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model_1, \
		max_F_dot = strain_rate, use_spatial_grid = False, ignore_yield = True, \
		nondimensional_timestep_suggestion = 1e-1).solve()
	results_2 = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model_2, \
		max_F_dot = strain_rate, use_spatial_grid = False, ignore_yield = True, \
		nondimensional_timestep_suggestion = 1e-1).solve()

	# Compute the total stress and place into a results tuple
	beta_sigma_over_n_2 = np.interp(results_1[1], results_2[1], results_2[5])
	beta_sigma_over_n_tot = x_p*beta_sigma_over_n_2 + (1 - x_p)*results_1[5]
	results = results_1[0], results_1[1], results_1[2], results_1[3], results_1[4], beta_sigma_over_n_tot

	# Plot the results
	plotter().plot_results(None, results, n_over_beta = 37.78, use_nominal = True)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_3/sigma(F).png?raw=true)
