# BoseHubbardDiagonalize

A Julia exact diagonalization code for the [Bose-Hubbard model](https://en.wikipedia.org/wiki/Bose%E2%80%93Hubbard_model) with a focus on entanglement entropy.

The basis is enumerated using the method of [Szabados et al., 2011](http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf).
See also [Zhang et al., 2011](http://arxiv.org/pdf/1102.4006v1.pdf).

Tested with Julia 1.3.


## Installation

```
pkg> add https://github.com/MelkoCollective/BH_diagonalize.git
```

In order to run the driver scripts in `bin/`, you will also need to
```
pkg> add ArgParse
pkg> add Arpack
pkg> add JeszenszkiBasis
pkg> add LsqFit
pkg> add ProgressMeter
pkg> add Qutilities
```

### Application project

If you're working with a clone of this repository, you can use the basic application project in `bin/`, which already has the necessary dependencies.
From the repository root, run
```
julia --project=bin
```
and then
```
pkg> dev .
```
to create `bin/Manifest.toml` with a development version of `BoseHubbardDiagonalize`.


## Examples

To run the following examples, you should set the project (e.g. using `--project` or `JULIA_PROJECT`) to a Julia project that has the prerequisites installed.

* `julia bin/diagonalize.jl --help`
* `julia bin/diagonalize.jl --out output.dat --ee 2 4 4`
* `julia bin/diagonalize.jl --out output.dat --site-max 4 --ee 3 6 6`

- - - -

* `julia bin/eop_peak_finder.jl --help`
* `julia bin/eop_peak_finder.jl --out output.dat --u-min 6 --u-max 10 4 4 2`

- - - -

* `julia bin/qubit_entanglement.jl --help`
* `julia bin/qubit_entanglement.jl --out output.dat 6 6 1 2 5 6`
