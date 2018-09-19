# BoseHubbardDiagonalize

A Julia exact diagonalization code for the [Bose-Hubbard model](https://en.wikipedia.org/wiki/Bose%E2%80%93Hubbard_model) with a focus on entanglement entropy.

The basis is enumerated using the method of [Szabados et al., 2011](http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf).
See also [Zhang et al., 2011](http://arxiv.org/pdf/1102.4006v1.pdf).

Tested with Julia 1.0.


## Installation

```
pkg> add https://github.com/MelkoCollective/BH_diagonalize.git
```

In order to run the driver scripts, you will also need to
```
pkg> add ArgParse
pkg> add Arpack
pkg> add JeszenszkiBasis
pkg> add LsqFit
pkg> add ProgressMeter
pkg> add Qutilities
```


## Examples

To run the following examples, you should set the project (e.g. using `--project` or `JULIA_PROJECT`) to a Julia project that has the prerequisites installed.

* `julia BH_main.jl --help`
* `julia BH_main.jl --out output.dat --ee 2 4 4`
* `julia BH_main.jl --out output.dat --site-max 4 --ee 3 6 6`

- - - -

* `julia tools/eop_peak_finder.jl --help`
* `julia tools/eop_peak_finder.jl --out output.dat --u-min 6 --u-max 10 4 4 2`

- - - -

* `julia tools/qubit_entanglement.jl --help`
* `julia tools/qubit_entanglement.jl --out output.dat 6 6 1 2 5 6`
