# Example 2

	# Import the library
	from chain_breaking_polymer_networks import *

	# Create the single-chain model
	single_chain_model = ideal()

	# Create the relaxation function
	relaxation_function = Long_et_al_2014(alpha = 2.6, t_R = 0.6, x_p = 0.1)

	# Define the deformation
	strain_rate = 0.03
	maximum_strain = 1
	total_time_in_seconds = 57
	def F(t):
		return 1 + strain_rate*t*np.heaviside(maximum_strain - strain_rate*t, 0.5) \
			+ (2*maximum_strain - strain_rate*t)*np.heaviside(strain_rate*t - maximum_strain, 0.5)

	# Apply the deformation and solve
	network_model = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, \
		relaxation_function = relaxation_function, use_spatial_grid = False, ignore_yield = True)
	results = network_model.solve()

	# Plot the results
	plotter().plot_results(network_model, results, n_over_beta = 24.15, use_nominal = True)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_2/sigma_hpt(F).png?raw=true)
