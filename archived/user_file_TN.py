# Import the package
from chain_breaking_polymer_networks import *

# Create the single-chain model
single_chain_model = Morse_FJC(N_b = 1, N_b_H = 38, beta_u_b = 61.57, kappa = 9000, 
	kappa_H = 6000, k_0 = 1e-88, beta_Delta_Psi_0 = 5, varsigma = 4)

# Volumetric swelling ratio
J_sw = 15.625

# Plot the single-chain model
output_directory = 'output_TN/'
plotter_object = plotter(plot_directory = output_directory)
plotter_object.plot_single_chain(single_chain_model, J_sw = J_sw, csv_directory = output_directory)

# Define the deformation
lambda_dot = 0.025
total_time_in_seconds = 15/lambda_dot
def F(t):
	return 1 + lambda_dot*t*np.heaviside(0.5 - lambda_dot*t, 0.5) \
		+ (1 - lambda_dot*t)*np.heaviside(lambda_dot*t - 0.5, 0.5)*np.heaviside(1 - lambda_dot*t, 0.5) \
		+ (lambda_dot*t - 1)*np.heaviside(lambda_dot*t - 1, 0.5)*np.heaviside(2 - lambda_dot*t, 0.5) \
		+ (3 - lambda_dot*t)*np.heaviside(lambda_dot*t - 2, 0.5)*np.heaviside(3 - lambda_dot*t, 0.5) \
		+ (lambda_dot*t - 3)*np.heaviside(lambda_dot*t - 3, 0.5)*np.heaviside(4.5 - lambda_dot*t, 0.5) \
		+ (6 - lambda_dot*t)*np.heaviside(lambda_dot*t - 4.5, 0.5)*np.heaviside(6 - lambda_dot*t, 0.5) \
		+ (lambda_dot*t - 6)*np.heaviside(lambda_dot*t - 6, 0.5)*np.heaviside(8 - lambda_dot*t, 0.5) \
		+ (10 - lambda_dot*t)*np.heaviside(lambda_dot*t - 8, 0.5)*np.heaviside(10 - lambda_dot*t, 0.5) \
		+ (lambda_dot*t - 10)*np.heaviside(lambda_dot*t - 10, 0.5)*np.heaviside(12.5 - lambda_dot*t, 0.5) \
		+ (15 - lambda_dot*t)*np.heaviside(lambda_dot*t - 12.5, 0.5)*np.heaviside(15 - lambda_dot*t, 0.5)

# Plot the deformation
plotter_object.plot_deformation(F, total_time_in_seconds, csv_directory = output_directory)

# Apply the deformation and solve
deformation_object = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, J_sw = J_sw)
results = deformation_object.solve(csv_directory = output_directory)

# Combine with a Neo-Hookean model for the filler network
E = 1.5
n_over_beta_1 = 0.2
n_over_beta_23 = E/3 - n_over_beta_1
F = results[1]
sigma_1 = n_over_beta_1/single_chain_model.P_A_tot_eq*results[5]
F_stress_23 = F, n_over_beta_23*(F**2 - 1/F)
F_stress_tot = F, F_stress_23[1] + sigma_1

# Plot the results with the experimental data
filename = 'data/ducrot_fig4c_stress.csv'
data_F_stress = np.genfromtxt(filename, usecols = (0)), np.genfromtxt(filename, usecols = (1))
plotter_object.plot_results(deformation_object, results, n_over_beta = n_over_beta_1, \
	data_F_stress = data_F_stress, F_stress_1 = F_stress_23, F_stress_2 = F_stress_tot)