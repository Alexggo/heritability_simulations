using DataFrames
using Random
using Statistics


# Parameters
num_simulations = 100
num_groups = 30

number_loci = 6
min_num_hetero = 3

number_of_dom = 2
effect_A1 = +1
effect_A2 = 0
effect_d = 0.5

env_type = "Uniform"
env_sd = 2
env_effect = 1
epistasis = true
epistasis_level = 6
epistasis_value = +10
baseline = 10


if env_type == "Uniform"
  env_range = vcat(-3:-1, 1:3)
elseif env_type == "Normal"
  env_range = round.(randn(1000) .* env_sd, digits=2)
else
  println("error")
end

  
possible_num_homo = 0:(number_loci - min_num_hetero)
number_of_add = number_loci - number_of_dom
add_vec = fill("ADD", number_of_add)
dom_vec = fill("DOM", number_of_dom)
all_effect = vcat(add_vec, dom_vec)
all_effect = reshape(shuffle(all_effect),1,:)
mid_effect = mean([effect_A1, effect_A2])
effect_a1 = effect_A1 - mid_effect
effect_a2 = effect_A2 - mid_effect
minimum_pheno_value = baseline + env_effect * minimum(env_range) + 2 * number_loci * minimum([effect_A1, effect_A2])
  

# Function to sample genotypesfunction sample_genotype(num_loci)
function sample_genotype(num_loci)
    possible_num_homo = 0:(num_loci - min_num_hetero)
    hom = rand(possible_num_homo)
    het = num_loci - hom
    hom_vec = reshape(fill("HOM", hom),1,:)
    het_vec = reshape(fill("HET", het),1,:)
    all_loci = hcat(hom_vec, het_vec)
    all_loci = shuffle(all_loci)
    chromosome1 = reshape(rand(["A1", "A2"], num_loci), 1, :)
    chromosome2 = reshape(Vector{String}(undef, num_loci),1,:)
    
    for i in 1:length(all_loci)
      if all_loci[i] == "HOM"
          chromosome2[i] = chromosome1[i]
      elseif all_loci[i] == "HET"
          chromosome2[i] = (chromosome1[i] == "A1" ? "A2" : "A1")
      else
          println("error")
      end
  end
    
    genotype = vcat(chromosome1, chromosome2)
    return genotype
end

function get_phenotype(geno, baseline, effect_a1, effect_a2, effect_d)
  effect_val = []
  for e in 1:length(all_effect) # For every loci
      if all_effect[e] == "DOM" # If the locus is dominant
          if geno[1, e] == "A1" && geno[2, e] == "A1"
              push!(effect_val, 2 * (mid_effect + effect_a1))
          elseif (geno[1, e] == "A1" && geno[2, e] == "A2") || (geno[1, e] == "A2" && geno[2, e] == "A1")
              push!(effect_val, mid_effect + effect_d + mid_effect + effect_d)
          elseif geno[1, e] == "A2" && par2[2, e] == "A2"
              push!(effect_val, 2 * (mid_effect + effect_a2))
          else
              println("error")
          end
      elseif all_effect[e] == "ADD" # If the locus is additive
          if geno[1, e] == "A1" && geno[2, e] == "A1"
              push!(effect_val, 2 * (mid_effect + effect_a1))
          elseif (geno[1, e] == "A1" && geno[2, e] == "A2") || (geno[1, e] == "A2" && geno[2, e] == "A1")
              push!(effect_val, mid_effect + effect_a1 + mid_effect + effect_a2)
          elseif geno[1, e] == "A2" && geno[2, e] == "A2"
              push!(effect_val, 2 * (mid_effect + effect_a2))
          else
              println("error")
          end
      else
          println("error")
      end
  end

  # Epistasis term
  if epistasis == true # If epistasis is present
      epis_val = sum(effect_val .>= 1) >= epistasis_level
      geno_val = epis_val ? sum(effect_val) + epistasis_value : sum(effect_val)
  else
      # No epistasis
      geno_val = sum(effect_val)
  end

  env_val = rand(env_range) # Sample the environmental range
  pheno = baseline + env_effect * env_val + geno_val
  return [baseline, env_val, geno_val, pheno]
end

# Function to get offspring given two parental genotypes
function get_offspring(genotype1, genotype2)
  gamete1 = []
  for i in 1:size(genotype1, 2)
      if genotype1[1, i] == genotype1[2, i]
          push!(gamete1, genotype1[1, i])
      elseif genotype1[1, i] != genotype1[2, i]
          push!(gamete1, sample([genotype1[1, i], genotype1[2, i]], 1))
      end
  end

  gamete2 = []
  for i in 1:size(genotype2, 2)
      if genotype2[1, i] == genotype2[2, i]
          push!(gamete2, genotype2[1, i])
      elseif genotype2[1, i] != genotype2[2, i]
          push!(gamete2, sample([genotype2[1, i], genotype2[2, i]], 1))
      end
  end

  new_genotype = [gamete1'; gamete2']
  return new_genotype
end




# Function to calculate phenotype given genotype and other parameters
slope_vec = []
intercept_vec = []
R_squared_vec = []
plot_vec = []

for simulation in 1:num_simulations
  all_data = DataFrame()

  for group in 1:num_groups
  all_geno = Dict(
    "par1" => sample_genotype(number_loci),
    "par2" => sample_genotype(number_loci),
    "par3" => sample_genotype(number_loci),
    "par4" => sample_genotype(number_loci))

pheno_par1 = get_phenotype(all_geno["par1"], baseline, effect_a1, effect_a2, effect_d)
pheno_par2 =  get_phenotype(all_geno["par2"], baseline, effect_a1, effect_a2, effect_d)
pheno_par3 =  get_phenotype(all_geno["par3"], baseline, effect_a1, effect_a2, effect_d)
pheno_par4 =  get_phenotype(all_geno["par4"], baseline, effect_a1, effect_a2, effect_d)

[pheno_par1; pheno_par2]

pheno_table = DataFrame(pheno_par1=pheno_par1, pheno_par2=pheno_par2, pheno_par3=pheno_par3, pheno_par4=pheno_par4)
row.names!(pheno_table, ["par1", "par2", "par3", "par4"])


push!(all_data, pheno_table)




all_data_df = vcat(all_data...)
rename!(all_data_df, [:group, :largep_env1, :largep_geno1, :largep_pheno1,
                      :largep_env2, :largep_geno2, :largep_pheno2,
                      :large_midpar, :largeo_env1, :largeo_geno1, :largeo_pheno1,
                      :largeo_env2, :largeo_geno2, :largeo_pheno2,
                      :large_moff, :smallp_env1, :smallp_geno1, :smallp_pheno1,
                      :smallp_env2, :smallp_geno2, :smallp_pheno2,
                      :small_midpar, :smallo_env1, :smallo_geno1, :smallo_pheno1,
                      :smallo_env2, :smallo_geno2, :smallo_pheno2,
                      :small_midoff])

large_df = select(all_data_df, [:group, :largep_env1, :largep_geno1, :largep_pheno1,
                                :largep_env2, :largep_geno2, :largep_pheno2,
                                :large_midpar, :largeo_env1, :largeo_geno1, :largeo_pheno1,
                                :largeo_env2, :largeo_geno2, :largeo_pheno2,
                                :large_moff])
large_df[:type] = "large"

small_df = select(all_data_df, [:group, :smallp_env1, :smallp_geno1, :smallp_pheno1,
                                :smallp_env2, :smallp_geno2, :smallp_pheno2,
                                :small_midpar, :smallo_env1, :smallo_geno1, :smallo_pheno1,
                                :smallo_env2, :smallo_geno2, :smallo_pheno2,
                                :small_midoff])
small_df[:type] = "small"

all_data_df2 = vcat(large_df, small_df)


regression_lines_df = DataFrame(x=[], y=[], simulation=String[])
for i in 1:length(slope_vec)
  x_range = extrema(df[:midpar])
  y_range = slope_vec[i] .* x_range .+ intercept_vec[i]
  line_df = DataFrame(x=x_range, y=y_range, simulation="Simulation $i")
  append!(regression_lines_df, line_df)
end

using Gadfly
using StatsBase

# Plot the original points and regression lines
p1 = plot(df, x=:midpar, y=:midoff, Geom.point,
          Guide.xlabel("Midparent diameter (cm)"), Guide.ylabel("Midoffspring diameter (cm)"),
          Theme(panel_stroke=colorant"white", default_color=colorant"blue"),
          Geom.vline(xintercept=[10, minimum_pheno_value], style=[:solid :dash], color=[colorant"black" colorant"red"]),
          Geom.hline(yintercept=[10, minimum_pheno_value], style=[:solid :dash], color=[colorant"black" colorant"red"]))

# Print heritability
heredity_less_than_0 = count(x -> x <= 0, slope_vec) / length(slope_vec)
println("Heritability less than 0 = $heredity_less_than_0")

# Hexbin plot
scatter = plot(df2, x=:slope, y=:intercept, Geom.hexbin,
               Guide.xlabel("Slope"), Guide.ylabel("Intercept"),
               Theme(panel_stroke=colorant"white"))

# Histograms
hist_bottom = plot(layer(df2, x=:slope), Geom.histogram,
                   Guide.xlabel("Slope"), Guide.ylabel("Frequency"),
                   Theme(panel_stroke=colorant"white", default_color=colorant"blue"),
                   Coord.cartesian(xmin=0, xmax=1))
hist_right = plot(layer(df2, x=:intercept), Geom.histogram, Geom.transpose(),
                  Guide.xlabel("Frequency"), Guide.ylabel("Intercept"),
                  Theme(panel_stroke=colorant"white", default_color=colorant"blue"),
                  Coord.cartesian(ymin=0, ymax=1))
empty = plot(layer(x=[1], y=[1]), Geom.point, Geom.blank,
             Theme(panel_stroke=colorant"white"))

p2 = vstack(scatter, hstack(hist_right, hist_bottom, empty))

# Display plots
draw(PDF("p1.pdf", 6inch, 4inch), p1)
draw(PDF("p2.pdf", 6inch, 4inch), p2)