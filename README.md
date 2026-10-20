# SemiAnalyticWECs.jl: A Julia-based solver for surface/submerged piezoelectric wave energy converters

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://zjwegert.github.io/SemiAnalyticWECs.jl/dev)

Codes related to the following papers:
> Zachary J. Wegert, Ben Wilks, Ngamta Thamwattana, Vivien J. Challis, Santanu Koley, Michael H. Meylan (2025). "Wave energy conversion by floating and submerged piezoelectric bimorph plates". Ocean Engineering 354, 124792.

> Ben Wilks, Zachary J. Wegert, Michael H. Meylan, Vivien J. Challis and Ngamta Thamwattana (2026). "Water wave scattering by a submerged anisotropic elastic plate". Under consideration in JFM.

## Installation
The package can be installed using the Julia package manager: type `]` in the Julia REPL to enter package mode, and install using
```julia
pkg> add https://github.com/zjwegert/SemiAnalyticWECs.jl
```

## Usage
The numerical experiments for the above papers can be found in `/scripts/...`:
- `/scripts/Wegert_et_al_2025_OE/` contains the scripts for the Ocean Engineering paper.
- `/scripts/Wilks_et_al_2026_JFM/` contains the scripts for the JFM paper.