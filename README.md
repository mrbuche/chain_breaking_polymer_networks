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
* the net forward reaction rate coefficient function, `k(gamma)`

Example: for the ideal chain model above,

	single_chain_model.P_A_eq(0.88)

returns the probability density that a chain is both intact and at a nondimensional end-to-end length of 0.88.

## network

some suitable deformation, as well as the total time in seconds (or use a deformation rate of unity and the time is nondimensional) and the deformation mode (BC) (currently supports uniaxial and equibiaxial stresses)

	network_model = deform_network(F, 'uniaxial', 3, single_chain_model)
	
optional keyword arguments later
this is where all the magic happens -- complicated solvers, memory shit, more on methods later too

	results = network_model.solve(csv_directory = './')
	
optional keyword argument csv_directory to write a .csv file with the results

## plotting

The object used for quick plotting is created using

	plotter_object = plotter(plot_directory = './')

where the default for the optional keyword argument is that shown (the current directory). The plotter() class is the only instance requiring matplotlib in this package. The single-chain model functions can be plotted using

	plotter_object.plot_single_chain(single_chain_model)
	
and all the results from deforming the network using

	plotter_object.plot_results(network_model, results)
	
If only the nondimensional stress-stretch response of the network is desired, one can use

	plotter_object.plot_results(None, results)
	
which would then output:
	
![alt text](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_0/sigma(F).png?raw=true)

# Examples

## Example 1
