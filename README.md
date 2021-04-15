# chain_breaking_polymer_networks

This is the Python package corresponding to "Chain breaking in the statistical mechanical constitutive theory of polymer networks" by Michael R. Buche and Meredith N. Silberstein, 2021.

It was written for Python 3, and uses some typical packages: numpy, scipy, and matplotlib.

# Installation

The package is on [PyPI](https://pypi.org/project/chain-breaking-polymer-networks/) and may be installed using `pip`:

	pip install chain_breaking_polymer_networks

## Extended installation

On a Linux machine, start by updating and then installing Python:

	sudo apt-get update
	sudo apt-get install python
	
Next, install the pip package manger,

	sudo apt-get install python-pip
	
and use `pip` to install the packages that this packages will use:

	pip install numpy, scipy, matplotlib
	
Finally, install this package:
	
	pip install chain_breaking_polymer_networks

# Basic usage

The package is best imported using:

	from chain_breaking_polymer_networks import *

## single_chain

The package contains a library of 'single_chain' classes corresponding to different single-chain models. For example, 

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

The `network` file contains a few classes, most notably the `deform_network` class. Given an applied deformation, such as

	def F(t): return 1 + t
	
the total testing time in seconds, the deformation mode (currently supports 'uniaxial' and 'equibiaxial' stress; the deformation is of course incompressible), and a, the `deform_network` class is used to create a network model from the single-chain model. Here we apply uniaxial stress for 3 seconds:

	network_model = deform_network(F, 'uniaxial', 3, single_chain_model, ignore_yield = True, use_spatial_grid = False)
	
Since this initialization also prepares the solution method, many optional keyword arguments are available. Here we ignore the breaking of chains via meeting a yield surface at some critical extension (the ideal chain model is infinitely extensible), and we choose to utilize quadrature for spatial integrations rather than a pre-specified spatial grid; the converse in either case is the default. The examples that follow illustrate the critical aspects of creating the network model, and more information can be found in the helpful comments in the `network.py` file, as well as in the Appendix of our paper.

The results (stress, total probability that a chain is intact, etc.) are solved for over the specified testing time using

	results = network_model.solve(csv_directory = './')
	
where the optional keyword argument here is used to write the results to a .csv file. The default is None (no .csv file is written).

## plotting

The object used for quick plotting is created using

	plotter_object = plotter(plot_directory = './')

where the default for the optional keyword argument is shown (the current directory). The plotter() class is the only instance requiring matplotlib in this package. The single-chain model functions can be plotted using

	plotter_object.plot_single_chain(single_chain_model)
	
and all the results from deforming the network using

	plotter_object.plot_results(network_model, results)
	
If only the nondimensional stress-stretch response of the network is desired, one can use

	plotter_object.plot_results(None, results)

# Examples

## Example 1

![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_0/sigma(F).png?raw=true)
