# Import the package
from chain_breaking_polymer_networks import *

# Create the single-chain model
single_chain_model = Morse_FJC(N_b = 1, N_b_H = 8, kappa = 200, kappa_H = 500, \
	beta_u_b = 100, k_0 = 1e-2, varsigma = 1, beta_Delta_Psi_0 = 5)

# Plot the single-chain model
output_directory = 'baseline_'
plotter_object = plotter(plot_directory = output_directory)
plotter_object.plot_single_chain(single_chain_model)

# Define the deformation
Delta = 8
total_time_in_seconds = 2*Delta
def F(t):
	return 1 + t*np.heaviside(Delta - t, 0.5) + (2*Delta - t)*np.heaviside(t - Delta, 0.5)

# Plot the deformation
plotter_object.plot_deformation(F, total_time_in_seconds)

# Apply the deformation and solve
deformation_object = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, num_grid_suggestion = 1025)
results = deformation_object.solve(csv_directory = output_directory)

# Plot the results
plotter_object.plot_results(deformation_object, results)