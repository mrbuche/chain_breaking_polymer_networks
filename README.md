# chain_breaking_polymer_networks

[![DOI:10.5281/zenodo.4699349](https://mbuche.github.io/web/badges/zenodo.4699349.svg)](https://doi.org/10.5281/zenodo.4699349) &nbsp; [![GitHub release](https://mbuche.github.io/web/badges/releasev1.0.1.svg)](https://github.com/mbuche/chain_breaking_polymer_networks/releases/) &nbsp; [![PyPI pyversions](https://mbuche.github.io/web/badges/python3.7.svg)](https://pypi.org/project/chain-breaking-polymer-networks/) &nbsp; [![GitHub license](https://mbuche.github.io/web/badges/licenseMIT.svg)](https://github.com/mbuche/chain_breaking_polymer_networks/blob/master/LICENSE)

This is the Python package corresponding to: Buche, Michael R., and Meredith N. Silberstein. Chain breaking in the statistical mechanical constitutive theory of polymer networks. Journal of the Mechanics and Physics of Solids (2021).

[![DOI:10.1016/j.jmps.2021.104593](https://mbuche.github.io/web/badges/j.jmps.2021.104593.svg)](https://doi.org/10.1016/j.jmps.2021.104593) 

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
	
the total testing time in seconds, the deformation mode (currently supports 'uniaxial' and 'equibiaxial' stress; the deformation is of course incompressible), and a single-chain model, the `deform_network` class is used to create a network model from the single-chain model. Here we apply uniaxial stress for 3 seconds:

	network_model = deform_network(F, 'uniaxial', 3, single_chain_model, ignore_yield = True, use_spatial_grid = False)
	
Since this initialization also prepares the solution method, many optional keyword arguments are available. In this example we ignore the breaking of chains via meeting a yield surface at some critical extension (the ideal chain model is infinitely extensible) using `ignore_yield = True`, and we choose to utilize quadrature for spatial integrations rather than a pre-specified spatial grid using `use_spatial_grid = False`; the converse in either case is the default. The [Examples](#examples) that follow illustrate the critical aspects of creating the network model, and more information can be found in the helpful comments in the `network.py` file, as well as in the Appendix of our [paper](https://arxiv.org/abs/2104.08866).

The results (stress, total probability that a chain is intact, etc.) are solved for over the specified testing time using

	results = network_model.solve(csv_directory = './')
	
where the optional keyword argument here is used to write the results to a .csv file. The default is None (no .csv file is written).

## relaxation_function

The `relaxation_function` module contains several different classes corresponding to different relaxation functions, `.g(t, tau)`, as well as their derivatives, `.d_g_d_tau(t, tau)`, and their corresponding loss and storage functions, `.g_p(t, tau)` and `.g_pp(t, tau)`. For example, 

	sticky_Rouse(N_b = 50, N_x = 5, t_0 = 4e-4, beta_E_A = 10)
	
returns an object that is the sticky Rouse model with 50 Kuhn monomers, 5 cross-links per chain, a Kuhn monomer relaxation time of 0.0004 seconds, and a nondimensional cross-link dissociation energy of 10 (read more [here](https://dx.doi.org/10.1122/1.4818868) or [here](https://dx.doi.org/10.1039/D0SM01115K)). The relaxation function is an optional keyword argument when creating the network model (the default is None); see our [paper](https://arxiv.org/abs/2104.08866) for more details.

## plotting

The `plotting` module allows object for plotting to be quickly created, i.e.

	plotter_object = plotter(plot_directory = './')

where the default for the optional keyword argument is shown (the current directory). The `plotter()` class is the only instance requiring `matplotlib` in this package. The single-chain model functions can be plotted using

	plotter_object.plot_single_chain(single_chain_model)
	
and all the results from deforming the network using

	plotter_object.plot_results(network_model, results)
	
If only the nondimensional stress-stretch response of the network is desired, one can use

	plotter_object.plot_results(None, results)

# Examples

* [Example 1](/examples/example_1/#readme)
* [Example 2](/examples/example_2/#readme)
* [Example 3](/examples/example_3/#readme)
* [Example 4](/examples/example_4/#readme)
