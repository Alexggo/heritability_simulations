# Heritability simulations

This repository contains a julia notebook for performing simulations of genetic and environmental effects on phenotypic traits.

# What is this script for?
Simulation of genetic variation: The script simulates parental genotypes for a specified number of loci, considering both additive and dominant effects. It then generates offspring genotypes based on parental genotypes and calculates resulting phenotypic values.
Environmental variation: Environmental effects on phenotypic traits are simulated by specifying a range of environmental values. Phenotypic values are calculated by incorporating both genetic and environmental effects.
Regression analysis: Linear regression models are fitted to the simulated data to explore the relationship between parental and offspring phenotypic values. The script calculates regression coefficients, intercepts, and R-squared values to quantify the strength of the relationship.
Visualization: The script generates visualizations, including scatter plots of parental-offspring phenotypic values, histograms of regression coefficients and intercepts, and a combined plot displaying regression lines for multiple simulations.
The script provides insights into the interaction between genetic and environmental factors in determining phenotypic variation, making it useful for educational purposes and research in quantitative genetics. And includes a model for epistasis and different environmental ranges.

# How can I run this?

1. Download the html file to your computer.
2. Run in binder or locally. I recommend to run the script locally in your computer.
3. Download the html file and open it in a browser.
4. Install [julia](https://julialang.org/downloads/) in your computer
5. Open the julia terminal and install Pluto:
```julia
using Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```
5. This code should open pluto in your browser.
6. Select the heritability notebook heritability_simulations.jl from the list.
7.  Toggle the different parameters to compute the heritability and generate plots.

For more details and usage instructions, please refer to the pdf, HTML or the script itself heritability1.jl (for julia/pluto).
