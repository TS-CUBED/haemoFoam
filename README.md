# haemoFoam
Haemodynamics simulation framework based on OpenFOAM. Includes particle migration model and advanced haemorheology models, as well as windkessel boundary conditions.

This is the public front for the haemoFoam code. There is no public code here at the moment. If you are interested in using this code for your own work, please contact me with a short description of what you plan to do, and I will give you access to the private development repository. I just want to keep an overview over who is using this.

The current development is on OpenFOAM 7, but the haemodynamics/haemorheology model will work on foam-extend. If you need it on OpenFoam v1912, I can port it, but will need a few days to do that (mostly library paths, I assume).

Reference:
Schenkel, T. and Halliday, I., 2020. Continuum Scale Non Newtonian Particle Transport Model for Haemorheology--Implementation and Validation. arXiv preprint arXiv:2004.12380. https://arxiv.org/abs/2004.12380
