# Import the library
from chain_breaking_polymer_networks import *

# Create the single-chain model
single_chain_model = Morse_FJC(N_b = 1, N_b_H = 8, kappa = 2e2, kappa_H = 5e2, beta_u_b = 1e2, k_0 = 1e-2, beta_Delta_Psi_0 = 2)

# Plot the single-chain model
plotter_object = plotter()
plotter_object.plot_single_chain(single_chain_model)

# Define the deformation
strain_rate = 1
maximum_strain = 8
total_time_in_seconds = 2*maximum_strain/strain_rate
def F(t):
	return 1 + strain_rate*t*np.heaviside(maximum_strain - strain_rate*t, 0.5) \
		+ (2*maximum_strain - strain_rate*t)*np.heaviside(strain_rate*t - maximum_strain, 0.5)

# Plot the deformation
plotter_object.plot_deformation(F, total_time_in_seconds)

# Apply the deformation and solve
network_model = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, num_grid_suggestion = 513)
results = network_model.solve(csv_directory = './')

# Plot the results
plotter_object.plot_results(network_model, results)
