# chain_breaking_polymer_networks

[![DOI:10.5281/zenodo.4699349](https://mbuche.github.io/web/badges/zenodo.4699349.svg)](https://doi.org/10.5281/zenodo.4699349) &nbsp; [![GitHub release](https://mbuche.github.io/web/badges/releasev1.0.1.svg)](https://github.com/mbuche/chain_breaking_polymer_networks/releases/) &nbsp; [![PyPI pyversions](https://mbuche.github.io/web/badges/python3.7.svg)](https://pypi.org/project/chain-breaking-polymer-networks/) &nbsp; [![GitHub license](https://mbuche.github.io/web/badges/licenseMIT.svg)](https://github.com/mbuche/chain_breaking_polymer_networks/blob/master/LICENSE)

This is the Python package corresponding to "Chain breaking in the statistical mechanical constitutive theory of polymer networks" by Michael R. Buche and Meredith N. Silberstein, 2021.

[![arXiv:2104.08866](https://mbuche.github.io/web/badges/badgearXiv210408866.svg)](https://arxiv.org/abs/2104.08866) 

# Installation

The package is on [PyPI](https://pypi.org/project/chain-breaking-polymer-networks/) and can be installed using `pip`:

	pip install chain_breaking_polymer_networks

It was written for Python 3, and uses some typical packages: `numpy`, `scipy`, and `matplotlib`.

# Basic usage

The package is best imported using:

	from chain_breaking_polymer_networks import *

## single_chain

The package contains a `single_chain` module of classes corresponding to different single-chain models. For example, 

	single_chain_model = ideal(N_b = 88)
	
creates an ideal chain model with 88 links. The single-chain models are nearly fully nondimensional (model parameters, inputs and outputs for functions) and they use keyword arguments for the model parameters, though not all parameters are optional (see source code, or examples). The single-chain models contain four main functions:
* the nondimensional single-chain mechanical response, `eta(gamma)`
* the nondimensional equilibrium probability density distribution of intact chain extensions, `P_A_eq(gamma)`
* the nondimensional equilibrium radial distribution function, `g_A_eq(gamma)`
* the net forward reaction rate coefficient function, `k(gamma)` (is the only function with units, 1/seconds)

Example: for the ideal chain model above,

	single_chain_model.P_A_eq(0.88)

returns the probability density that a chain is both intact and at a nondimensional end-to-end length of 0.88.

## network

The `network` module contains a few classes, most notably the `deform_network` class. Given an applied deformation, such as

	def F(t): return 1 + t
	
the total testing time in seconds, the deformation mode (currently supports 'uniaxial' and 'equibiaxial' stress; the deformation is of course incompressible), and a, the `deform_network` class is used to create a network model from the single-chain model. Here we apply uniaxial stress for 3 seconds:

	network_model = deform_network(F, 'uniaxial', 3, single_chain_model, ignore_yield = True, use_spatial_grid = False)
	
Since this initialization also prepares the solution method, many optional keyword arguments are available. In this example we ignore the breaking of chains via meeting a yield surface at some critical extension (the ideal chain model is infinitely extensible) using `ignore_yield = True`, and we choose to utilize quadrature for spatial integrations rather than a pre-specified spatial grid using `use_spatial_grid = False`; the converse in either case is the default. The [Examples](#examples) that follow illustrate the critical aspects of creating the network model, and more information can be found in the helpful comments in the `network.py` file, as well as in the Appendix of our [paper](https://arxiv.org/abs/2104.08866).

The results (stress, total probability that a chain is intact, etc.) are solved for over the specified testing time using

	results = network_model.solve(csv_directory = './')
	
where the optional keyword argument here is used to write the results to a .csv file. The default is None (no .csv file is written).

## relaxation_function

The `relaxation_function` module contains several different classes corresponding to different relaxation functions, `.g(t, tau)`, as well as their derivatives, `.d_g_d_tau(t, tau)`, and their corresponding loss and storage functions, `.g_p(t, tau)` and `g_pp(t, tau)`. For example, 

	sticky_Rouse(N_b = 50, N_x = 5, t_0 = 4e-4, beta_E_A = 10)
	
returns an object that is the sticky Rouse model (read more [here](https://dx.doi.org/10.1122/1.4818868) or [here](https://dx.doi.org/10.1039/D0SM01115K)) with 50 Kuhn monomers, 5 cross-links per chain, a Kuhn monomer relaxation time of 0.0001 seconds, and a nondimensional cross-link dissociation energy of 10. The relaxation function is an optional keyword argument when creating the network model (the default is None); see our [paper](https://arxiv.org/abs/2104.08866) for more details.

## plotting

The `plotting` module alows object for plotting to be quickly created, i.e.

	plotter_object = plotter(plot_directory = './')

where the default for the optional keyword argument is shown (the current directory). The plotter() class is the only instance requiring matplotlib in this package. The single-chain model functions can be plotted using

	plotter_object.plot_single_chain(single_chain_model)
	
and all the results from deforming the network using

	plotter_object.plot_results(network_model, results)
	
If only the nondimensional stress-stretch response of the network is desired, one can use

	plotter_object.plot_results(None, results)

# Examples

## Example 1

	# Import the library
	from chain_breaking_polymer_networks import *

	# Create the single-chain model
	single_chain_model = Morse_FJC(N_b = 1, N_b_H = 8, kappa = 2e2, kappa_H = 5e2, \
		beta_u_b = 1e2, k_0 = 1e-2, beta_Delta_Psi_0 = 2)

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

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/eta.png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/P_A_eq.png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/g_A_eq.png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/k.png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/sigma(F).png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/sigma(t).png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/P_A_tot(t).png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/d_P_A_tot_dt(t).png?raw=true)

## Example 2

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

## Example 3

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

## Example 4

	# Import the library
	from chain_breaking_polymer_networks import *

	# Create the single-chain model
	single_chain_model = Morse_FJC(N_b = 1, N_b_H = 38, beta_u_b = 61.57, kappa = 9e3, 
		kappa_H = 6e3, k_0 = 0, beta_Delta_Psi_0 = 5, varsigma = 4)

	# Volumetric swelling ratio
	J_sw = 15.625

	# Plot the single-chain model
	plotter_object = plotter()
	plotter_object.plot_single_chain(single_chain_model, J_sw = J_sw)

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
	plotter_object.plot_deformation(F, total_time_in_seconds)

	# Apply the deformation and solve
	deformation_object = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, J_sw = J_sw)
	results = deformation_object.solve()

	# Combine the results with a Neo-Hookean model for the filler network
	modulus_TN = 1.5
	n_over_beta_1 = 0.2
	n_over_beta_23 = modulus_TN/3 - n_over_beta_1
	F = results[1]
	sigma_1 = n_over_beta_1/single_chain_model.P_A_tot_eq*results[5]
	F_stress_23 = F, n_over_beta_23*(F**2 - 1/F)
	F_stress_tot = F, F_stress_23[1] + sigma_1

	# Plot the results with the experimental data
	plotter_object.plot_results(deformation_object, results, n_over_beta = n_over_beta_1, \
		F_stress_1 = F_stress_23, F_stress_2 = F_stress_tot)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_4/F(t).png?raw=true)

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_4/sigma(F).png?raw=true)
