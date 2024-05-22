# Heritability simulations

This repository contains a julia notebook for performing simulations of genetic and environmental effects on phenotypic traits.

# What is this script for?
Simulation of genetic variation: The script simulates parental genotypes for a specified number of loci, considering both additive and dominant effects. It then generates offspring genotypes based on parental genotypes and calculates resulting phenotypic values.
Environmental variation: Environmental effects on phenotypic traits are simulated by specifying a range of environmental values. Phenotypic values are calculated by incorporating both genetic and environmental effects.
Regression analysis: Linear regression models are fitted to the simulated data to explore the relationship between parental and offspring phenotypic values. The script calculates regression coefficients, intercepts, and R-squared values to quantify the strength of the relationship.
Visualization: The script generates visualizations, including scatter plots of parental-offspring phenotypic values, histograms of regression coefficients and intercepts, and a combined plot displaying regression lines for multiple simulations.
The script provides insights into the interaction between genetic and environmental factors in determining phenotypic variation, making it useful for educational purposes and research in quantitative genetics. And includes a model for epistasis and different environmental ranges.

# Prerequirements:

This script requires two programs to be installed in your computer. Install these two programs first.

* [R](https://www.r-project.org/)

* [julia](https://julialang.org/downloads/)

# How can I run this?
1. Download the entire repository to your computer from GitHub.
2. Open your terminal and navigate to the repository that you downloaded.
3. Open the julia terminal and install [Rcall](https://juliainterop.github.io/RCall.jl/stable/installation/) and Pluto by typing the following to your terminal:
```julia
using Pkg

# You have to specify the directory where your R is located in your computer:
# You can figure this out by typing R.home() within the program R.
# In my case, R is installed in "C:/PROGRA~1/R/R-44~1.0"
ENV["R_HOME"] = "....directory of R home...."

Pkg.build("RCall")
Pkg.add("Pluto")
using Pluto
Pluto.run()
```
This code should open pluto in your browser.
3. Select the heritability notebook heritability_simulations.jl from the list.
4.  Toggle the different parameters to compute the heritability and generate plots.

For more details and usage instructions, please refer to the pdf, HTML or the script itself heritability1.jl (for julia/pluto).
