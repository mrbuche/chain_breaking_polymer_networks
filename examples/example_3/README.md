# Example 3

This is an example containing several features of this package.

[user_file.py](user_file.py) contains the full code. Begin by importing the package:

	from chain_breaking_polymer_networks import *
	
## Single-chain models

Create two single-chain models, both using the extensible freely-jointed chain (EFJC) model with 50 links and a nondimensional stiffness of 40 (read more [here](https://dx.doi.org/10.1016/j.physa.2019.121929), [here](https://dx.doi.org/10.1103/PhysRevE.102.012501), and [here](https://arxiv.org/abs/2104.08866)). The first model has a constant reaction rate coefficient function of 0.12/s (when specifying the intitial reaction rate coefficient for the ideal or EFJC models, the transient network model is employed). The second model has permanent chains; the fraction of permanent chains is defined here as 5.85%:

	x_p = 5.85e-2
	single_chain_model_1 = EFJC(N_b = 50, kappa = 40, k_0 = 0.12)
	single_chain_model_2 = EFJC(N_b = 50, kappa = 40)

## Deformation

The stretch for monotonic uniaxial tension is applied at a strain rate of 0.01/s to a maximum stretch of 2.55:

	strain_rate = 0.01
	maximum_strain = 1.55
	total_time_in_seconds = maximum_strain/strain_rate
	def F(t): return 1 + strain_rate*t
	
## Network and results

Here we create a network using each single-chain model and solve them separately (parallel configuration) for the relevant results in a single line. We also suggest decreasing the nondimensional timestep to 0.1 (default is 0.01), specify the maximum strain rate (helps with determining an appropriate minimum timestep; default is to estimate it from the given deformation), choose to use quadrature over a grid, and ignore any limiting intact chain extension (EFJC model is infinitely extensible):

	results_1 = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model_1, max_F_dot = strain_rate, use_spatial_grid = False, ignore_yield = True, nondimensional_timestep_suggestion = 1e-1).solve()
	results_2 = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model_2, max_F_dot = strain_rate, use_spatial_grid = False, ignore_yield = True, nondimensional_timestep_suggestion = 1e-1).solve()
	
Although we suggested to utilize the same timestep for each model, generally (due to difference in parameters) the results will not occur at the same times. Therefore before adding the stresses from each network (parallel configuration) and creating a new `results` tuple, we interpolate:

	beta_sigma_over_n_2 = np.interp(results_1[1], results_2[1], results_2[5])
	beta_sigma_over_n_tot = x_p*beta_sigma_over_n_2 + (1 - x_p)*results_1[5]
	results = results_1[0], results_1[1], results_1[2], results_1[3], results_1[4], beta_sigma_over_n_tot

These results are then plotted, where using `None` in place of `network_model` simplifies the plot output; we also specify the shear modulus and choose to show the nominal (engineering) stress:

	plotter().plot_results(None, results, n_over_beta = 37.78, use_nominal = True)

![sigma(F)](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_3/sigma(F).png)
