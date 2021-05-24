# Import the package
from chain_breaking_polymer_networks import *

# Create the single-chain model
single_chain_model = Morse_FJC(N_b = 1, N_b_H = 9, kappa = 200, kappa_H = 400, \
	beta_u_b = 100, k_0 = 2e-4, beta_Delta_Psi_0 = 8.55)

# Plot the single-chain model
output_directory = 'output_Cu/1.5_'
plotter_object = plotter(plot_directory = output_directory)
plotter_object.plot_single_chain(single_chain_model)

# Define the deformation
Delta = 1.5
lambda_dot = 8.4/60
total_time_in_seconds = 2*Delta/lambda_dot
def F(t):
	return 1 + lambda_dot*t*np.heaviside(total_time_in_seconds/2 - t, 0.5) \
		+ (2*Delta - lambda_dot*t)*np.heaviside(t - total_time_in_seconds/2, 0.5)

# Plot the deformation
plotter_object.plot_deformation(F, total_time_in_seconds)

# Apply the deformation and solve
deformation_object = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, max_F_dot = lambda_dot, \
	nondimensional_timestep_suggestion = 1e-1)
results = deformation_object.solve(csv_directory = output_directory)

# Plot the results with the experimental data
filename = 'data/lin_5b_all.csv'
data_F_stress = np.genfromtxt(filename, usecols = (0)), np.genfromtxt(filename, usecols = (1))
plotter_object.plot_results(deformation_object, results, \
	n_over_beta = 0.96/2, use_nominal = True, data_F_stress = data_F_stress)