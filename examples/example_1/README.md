# Example 1

This is simplest self-complete example of using this package. [user_file.py](user_file.py) contains the full code.


As always, begin by importing the package:

	from chain_breaking_polymer_networks import *
	
## Single-chain model

In executing the line:

	single_chain_model = Morse_FJC(N_b = 1, N_b_H = 8, kappa = 2e2, kappa_H = 5e2, beta_u_b = 1e2, k_0 = 1e-2, beta_Delta_Psi_0 = 2)

we create the Morse-FJC single-chain model consisting of:

* 1 breakable link with nondimensional stiffness 200 and energy 100, 
* 8 unbreakable links with nondimensional stiffness 500,
* an initial reaction rate coefficient of 0.01/s, and
* a nondimensional free energy change of 2 when breaking.

For more details on this single-chain model, see our [paper](https://arxiv.org/abs/2104.08866). 

The plotting object is then created and used to plot (saving the image in the local directory) each relevant single-chain function via:

	plotter_object = plotter()
	plotter_object.plot_single_chain(single_chain_model)
	
<table>
	<tr>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/eta.png" alt="eta(gamma)" style="width: 250px;"/>
		</td>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/P_A_eq.png" alt="P_A_eq(gamma)" style="width: 250px;"/>
		</td>
	</tr>
	<tr>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/g_A_eq.png" alt="g_A_eq(gamma)" style="width: 250px;"/>
		</td>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/k.png" alt="k(gamma)" style="width: 250px;"/>
		</td>
	</tr>
</table>
	
## Deformation



	strain_rate = 1
	maximum_strain = 8
	total_time_in_seconds = 2*maximum_strain/strain_rate
	def F(t): return 1 + strain_rate*t*np.heaviside(maximum_strain - strain_rate*t, 0.5) + (2*maximum_strain - strain_rate*t)*np.heaviside(strain_rate*t - maximum_strain, 0.5)

Using the same plotting object, 

	plotter_object.plot_deformation(F, total_time_in_seconds)
	
<table>
	<tr>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/F(t).png" alt="F(t)" style="width: 250px;"/>
		</td>
		<td></td>
	</tr>
</table>
	
## Network model, results

	# Apply the deformation and solve
	network_model = deform_network(F, 'uniaxial', total_time_in_seconds, single_chain_model, num_grid_suggestion = 513)
	results = network_model.solve(csv_directory = './')

	# Plot the results
	plotter_object.plot_results(network_model, results)
	
<table border="0">
	<tr>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/sigma(F).png" alt="sigma(F)" style="width: 250px;"/>
		</td>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/sigma(t).png" alt="sigma(t)" style="width: 250px;"/>
		</td>
	</tr>
	<tr>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/P_A_tot(t).png" alt="P_A_tot(t)" style="width: 250px;"/>
		</td>
		<td> 
			<img src="https://github.com/mbuche/chain_breaking_polymer_networks/blob/main/examples/example_1/d_P_A_tot_dt(t).png" alt="d_P_A_tot_dt(t)" style="width: 250px;"/>
		</td>
	</tr>
</table>
