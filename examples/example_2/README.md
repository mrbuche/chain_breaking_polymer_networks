# Example 2

This is an example of using the `relaxation_function` module of this package.

[user_file.py](user_file.py) contains the full code. Begin by importing the package:

	from chain_breaking_polymer_networks import *
	
## Single-chain model

Here we utilize the ideal single-chain model:

	single_chain_model = ideal()
	
## Relaxation function

The relaxation function from Long et al. 2014 (read more [here](https://dx.doi.org/10.1021/ma501290h), [here](https://dx.doi.org/10.1021/acs.macromol.6b00421), and [here](https://arxiv.org/abs/2104.08866)) is created with an exponent parameter 2.6, a characteristic bond breaking time of 2.6 seconds, and a 10% fraction of permanent chains:

	relaxation_function = Long_et_al_2014(alpha = 2.6, t_R = 0.6, x_p = 0.1)

## Deformation

One cycle of uniaxial tension will be applied with a strain rate of 0.03/s and a maximum stretch of 2:

	strain_rate = 0.03
	maximum_strain = 1
	total_time_in_seconds = 57
	def F(t): return 1 + strain_rate*t*np.heaviside(maximum_strain - strain_rate*t, 0.5) + (2*maximum_strain - strain_rate*t)*np.heaviside(strain_rate*t - maximum_strain, 0.5)

## Network model

The network model is created using the above specified single-chain model and relaxation function, while also choosing to use quadrature rather than the spatial grid for completing spatial integrals, and ignoring chain breaking via reaching any intact chain limiting extension (the ideal chain model is infinitely extensible):

	network_model = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, relaxation_function = relaxation_function, use_spatial_grid = False, ignore_yield = True)
		
## Results
	
The results are solved for; no .csv file is written since no directory is specified:
	
	results = network_model.solve()

In this case (`ignore_yield = True`) the homogeneous (from initially-intact chains) and particular (from reformed chains) solutions for the stress are returned in the `results` tuple since there are no breaking/reforming rates. When plotting the stress, these solutions are automatically shown. We also specify the shear modulus and choose to plot the nominal (engineering) stress:
	
	plotter().plot_results(network_model, results, n_over_beta = 24.15, use_nominal = True)

![sigma(F)](https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_2/sigma_hpt(F).png)
