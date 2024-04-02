# PricingBattery.jl

This repository uses SDDP.jl to compute the average utility
we get when using a battery to trade electricity on the spot market.

## Installation

The package uses Julia's package manager to install
all the dependencies:

```
julia --project -e "using Pkg; Pkg.instantiate()"

```

## Generate the results
Use the script `generate_results.jl` to generate the results:

```
julia --project generate_results.jl

```

All the results are dumped as text files in the directory `./results`.

## Plot the results

Use the script `plot_results.jl` to plot the results:
```
julia --project plot_results.jl

```
