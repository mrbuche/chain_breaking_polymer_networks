# Import the package
from chain_breaking_polymer_networks import *

# Create the single-chain model
single_chain_model = ideal()

# Create the relaxation function
relaxation_function = Long_et_al_2014(alpha = 2.6, t_R = 0.6, x_p = 0.1)

# Plot the single-chain model
output_directory = 'output_Long_relax/'
plotter_object = plotter(plot_directory = output_directory)
plotter_object.plot_single_chain(single_chain_model)

# Define the deformation
Delta_lambda = 1
lambda_dot = 0.03
total_time_in_seconds = 2*Delta_lambda/lambda_dot
def F(t):
	return 1 + lambda_dot*t*np.heaviside(Delta_lambda/lambda_dot - t, 0.5) \
		+ (2*Delta_lambda - lambda_dot*t)*np.heaviside(t - Delta_lambda/lambda_dot, 0.5)

# Plot the deformation
plotter_object.plot_deformation(F, total_time_in_seconds)

# Apply the deformation and solve
deformation_object = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, \
	relaxation_function = relaxation_function, max_F_dot = lambda_dot, use_spatial_grid = False, ignore_yield = True)
results = deformation_object.solve(csv_directory = output_directory)

# Plot the results with the experimental data
data_F_stress = np.genfromtxt('data/long2014_2e.csv', usecols = (0)), np.genfromtxt('data/long2014_2e.csv', usecols = (1))
plotter_object.plot_results(deformation_object, results, n_over_beta = 24.15, use_nominal = True, data_F_stress = data_F_stress)