# no_symmetry_twofold

This is a partially-finished side-development that will be applied to the source modules (`single_chain` and `network`) in the future. 
It does not rely on symmetry to either (1) compute spatial integrals, or (2) determine the unknown components of the deformation gradient. 
Instead, it employs a solver to iteratively compute the stress given the past deformation history and vary the unknown deformation gradient components until the traction boundary conditions are satisfied.
It is currently only applicable when using the EFJC model and rate-independent irreversible breaking, in addition to the deformation gradient existing in principal coordinates.
In the future it will be fully generalized, although certain cases will quickly become computationally prohibitive, so the original implementation is still highly recommended.
