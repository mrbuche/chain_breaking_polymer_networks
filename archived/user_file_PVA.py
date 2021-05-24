# Import the package
from chain_breaking_polymer_networks import *

# Create the single-chain model
single_chain_model = EFJC(N_b = 50, kappa = 40)

# Create the relaxation function
relaxation_function = Long_et_al_2014(alpha = 3.23, t_R = 1.99, x_p = 4.68e-2)

# Plot the single-chain model
output_directory = 'output_PVA/'
plotter_object = plotter(plot_directory = output_directory)
plotter_object.plot_single_chain(single_chain_model)

# Loop over all deformation rates
order = [4, 0, 1, 2, 3, 5]
lambda_dot_all   = [0.90, 0.10, 0.03, 0.01, 3e-3, 1e-4]
Delta_lambda_all = [0.65, 1.55, 2.66, 2.90, 3.90, 4.56]
for index in order:

	# Define the deformation
	lambda_dot = lambda_dot_all[index]
	Delta_lambda = Delta_lambda_all[index]
	total_time_in_seconds = Delta_lambda/lambda_dot
	def F(t):
		return 1 + lambda_dot*t

	# Take larger timesteps for very slow deformation rates
	if lambda_dot == 0.003:
		nondim_dt = 1e0
	elif lambda_dot == 0.0001:
		nondim_dt = 5e1
	else:
		nondim_dt = 1e-2

	# Output file names
	output_directory_index = output_directory + str(lambda_dot) + '_'
	plotter_object = plotter(plot_directory = output_directory_index)

	# Apply the deformation and solve
	deformation_object = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, \
		relaxation_function = relaxation_function, max_F_dot = lambda_dot, use_spatial_grid = False, ignore_yield = True, \
		nondimensional_timestep_suggestion = nondim_dt)
	results = deformation_object.solve(csv_directory = output_directory_index)

	# Plot the results with the experimental data
	filename = 'data/mayumi_3a_0p' + str(lambda_dot)[2:] + '.csv'
	data_F_stress = np.genfromtxt(filename, usecols = (0)), np.genfromtxt(filename, usecols = (1))
	plotter_object.plot_results(None, results, n_over_beta = 37.78, use_nominal = True, data_F_stress = data_F_stress)