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
	
creates an ideal chain model with 88 links. The single-chain models are nearly fully nondimensional (model parameters, inputs and outputs for functions) and they use keyword arguments for the model parameters, though not all parameters are optional (see source code, or examples). The single-chain models contain three main functions:
* the nondimensional single-chain mechanical response, `eta(gamma)`
* the nondimensional equilibrium distribution of intact chains, `P_A_eq(gamma)`
* the net forward reaction rate coefficient function, `k(gamma)`

## network

## plotting

## relaxation_function

# Examples

## Example 1
