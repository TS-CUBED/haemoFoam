# haemoFoam

Haemodynamics simulation framework based on OpenFOAM. Includes particle migration model and advanced haemorheology models, as well as windkessel boundary conditions.

The current development is on OpenFOAM v2206.

## If you want to use haemoFoam

HaemoFoam is published under the GPL (as is OpenFoam). You are, therefore, free to use it. I would ask you for three things, though:

1. Please reference my work if you use haemoFoam in your research and publications.
2. Please contribute back to haemoFoam and don't just fork it. I would hope we can develop this together in a mutually useful direction.
3. Please get in touch and let me know what you want to do. I an happy to collaborate.

## Reference:

Schenkel, T.; Halliday, I. Continuum Scale Non Newtonian Particle Transport Model for Hæmorheology. Mathematics 2021, 9, 2100. https://doi.org/10.3390/math9172100 

This is the preprint of the 2021 Mathematics paper:

Schenkel, T. and Halliday, I., 2020. Continuum Scale Non Newtonian Particle Transport Model for Haemorheology--Implementation and Validation. arXiv preprint arXiv:2004.12380. https://arxiv.org/abs/2004.12380

## Installation

Initialise OpenFOAM (e.g. `. /usr/lib/openfoam/openfoam2306/etc/bashrc`) and then run

`./Allwmake`

from the main `haemoFoam-0.n.n` directory. This will install all the solvers and models.

## Documentation

Documentation is not complete yet. See the `docs` folder. There are also a [few videos I recorded for my students.](https://www.youtube.com/playlist?list=PLWHQIdms-YHSVrFf5qchNdjX-lFuDj4kK)

## Acknowledgements:

The windkessel model is based on work by: Andris Piebalgs, Boram Gu, Emily Manchester, Imperial College London, who kindly let me have their [original code](https://github.com/KeepFloyding/OpenFOAM-phys-flow).




