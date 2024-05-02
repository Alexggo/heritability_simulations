# Heritability simulations

This repository contains an R script that performs simulations of genetic and environmental effects on phenotypic traits. The script utilizes the tidyverse, broom, and patchwork packages for data manipulation, regression analysis, and plotting, respectively.

The main features of the script include:

Simulation of genetic variation: The script simulates parental genotypes for a specified number of loci, considering both additive and dominant effects. It then generates offspring genotypes based on parental genotypes and calculates resulting phenotypic values.
Environmental variation: Environmental effects on phenotypic traits are simulated by specifying a range of environmental values. Phenotypic values are calculated by incorporating both genetic and environmental effects.
Regression analysis: Linear regression models are fitted to the simulated data to explore the relationship between parental and offspring phenotypic values. The script calculates regression coefficients, intercepts, and R-squared values to quantify the strength of the relationship.
Visualization: The script generates visualizations, including scatter plots of parental-offspring phenotypic values, histograms of regression coefficients and intercepts, and a combined plot displaying regression lines for multiple simulations.
The script provides insights into the interaction between genetic and environmental factors in determining phenotypic variation, making it useful for educational purposes and research in quantitative genetics.

For more details and usage instructions, please refer to the provided R script.
