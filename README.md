# turbinesFoam

[![DOI](https://zenodo.org/badge/4234/turbinesFoam/turbinesFoam.svg)](https://zenodo.org/badge/latestdoi/4234/turbinesFoam/turbinesFoam)
![OpenFOAM v2412](https://img.shields.io/badge/OpenFOAM-v2412-brightgreen.svg)
![OpenFOAM v2406](https://img.shields.io/badge/OpenFOAM-v2406-brightgreen.svg)
![OpenFOAM v2312](https://img.shields.io/badge/OpenFOAM-v2312-brightgreen.svg)
![OpenFOAM v2306](https://img.shields.io/badge/OpenFOAM-v2306-brightgreen.svg)
![OpenFOAM v2212](https://img.shields.io/badge/OpenFOAM-v2212-brightgreen.svg)

turbinesFoam is a library for simulating wind and marine hydrokinetic turbines
in OpenFOAM using the actuator line method.

[![](https://cloud.githubusercontent.com/assets/4604869/10141523/f2e3ad9a-65da-11e5-971c-b736abd30c3b.png)](https://www.youtube.com/watch?v=THZvV4R1vow)

Be sure to check out the
[development snapshot videos on YouTube](https://www.youtube.com/playlist?list=PLOlLyh5gytG8n8D3V1lDeZ3e9fJf9ux-e).

## Installation

### Docker

Spin up an interactive shell with:

```sh
docker run --rm -it -v $PWD:/work ghcr.io/turbinesfoam/turbinesfoam
```

### Compile from source

```sh
cd $WM_PROJECT_USER_DIR
git clone https://github.com/turbinesFoam/turbinesFoam.git
cd turbinesFoam
./Allwmake
```

## Usage

See the tutorials located in the `tutorials` directory.

## Contributing

Pull requests are very welcome!
See the [issue tracker](https://github.com/petebachant/turbinesFoam/issues)
for more details.

## Features

`fvOptions` classes for adding actuator lines and turbines constructed from
actuator lines to any compatible solver or turbulence model, e.g.,
`simpleFoam`, `pimpleFoam`, `interFoam`, etc.

## Publications

Bachant, P., Goude, A., and Wosnik, M. (2016) [_Actuator line modeling of vertical-axis turbines_](https://arxiv.org/abs/1605.01449). arXiv preprint 1605.01449.

## How to cite

The latest release of turbinesFoam can be cited via DOI thanks to Zenodo: [![DOI](https://zenodo.org/badge/4234/turbinesFoam/turbinesFoam.svg)](https://zenodo.org/badge/latestdoi/4234/turbinesFoam/turbinesFoam)

## Acknowledgements

This work was funded through a National Science Foundation CAREER award,
principal investigator Martin Wosnik ([NSF CBET
1150797](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1150797), Energy for
Sustainability, original program manager Geoffrey A. Prentice, current program
manager Gregory L. Rorrer).

OpenFOAM is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.

Interpolation, Gaussian projection, and vector rotation functions adapted from
NREL's [SOWFA](https://github.com/NREL/SOWFA).
