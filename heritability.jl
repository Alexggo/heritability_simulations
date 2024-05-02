using DataFrames
using Random
using Statistics
using GLM
using RCall
@rlibrary ggplot2

# Parameters
num_simulations = 500
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
epistasis = false
epistasis_level = 6
epistasis_value = +10
baseline = 10
mid_effect = mean([effect_A1, effect_A2])
effect_a1 = effect_A1 - mid_effect
effect_a2 = effect_A2 - mid_effect
possible_num_homo = 0:(number_loci - min_num_hetero)
number_of_add = number_loci - number_of_dom
add_vec = fill("ADD", number_of_add)
dom_vec = fill("DOM", number_of_dom)
all_effect = vcat(add_vec, dom_vec)

if env_type == "Uniform"
    env_range = vcat(-3:-1, 1:3)
elseif env_type == "Normal"
    env_range = round.(randn(1000) .* env_sd, digits=2)
else
println("error")
end

minimum_pheno_value = baseline + env_effect * minimum(env_range) + 2 * number_loci * minimum([effect_A1, effect_A2])

# Define functions:
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

function get_phenotype(geno, baseline, effect_a1, effect_a2, effect_d,effects)
effect_val = []
    for e in 1:length(effects) # For every loci
        if effects[e] == "DOM" # If the locus is dominant
            if geno[1, e] == "A1" && geno[2, e] == "A1"
                push!(effect_val, 2 * (mid_effect + effect_a1))
            elseif (geno[1, e] == "A1" && geno[2, e] == "A2") || (geno[1, e] == "A2" && geno[2, e] == "A1")
                push!(effect_val, mid_effect + effect_d + mid_effect + effect_d)
            elseif geno[1, e] == "A2" && geno[2, e] == "A2"
                push!(effect_val, 2 * (mid_effect + effect_a2))
            else
                println("error")
            end
        elseif effects[e] == "ADD" # If the locus is additive
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
    return [baseline env_val geno_val pheno]
end

# Function to get offspring given two parental genotypes
function get_offspring(genotype1, genotype2)
    gamete1 = []
    for i in 1:size(genotype1, 2)
        if genotype1[1, i] == genotype1[2, i]
            push!(gamete1, genotype1[1, i])
        elseif genotype1[1, i] != genotype1[2, i]
            push!(gamete1, rand([genotype1[1, i], genotype1[2, i]]))
        end
    end

    gamete2 = []
    for i in 1:size(genotype2, 2)
        if genotype2[1, i] == genotype2[2, i]
            push!(gamete2, genotype2[1, i])
        elseif genotype2[1, i] != genotype2[2, i]
            push!(gamete2, rand([genotype2[1, i], genotype2[2, i]]))
        end
    end

    gamete1=reshape(gamete1,1,:)
    gamete2=reshape(gamete2,1,:)

    new_genotype = vcat(gamete1, gamete2)
    return new_genotype
end


slope_vec = Vector{Float64}()
intercept_vec = Vector{Float64}()
R_squared_vec =Vector{Float64}()
plot_vec = []
all_sim = []
p1 = []
for simulation in 1:num_simulations
    println("Simulation: $simulation")
    all_data = [] # List with all_data per simulation

for group in 1:num_groups
    println("Group: $group")
      
    geno_par = Dict(
        "par1" => sample_genotype(number_loci),
        "par2" => sample_genotype(number_loci),
        "par3" => sample_genotype(number_loci),
        "par4" => sample_genotype(number_loci))

    pheno_par1 = get_phenotype(geno_par["par1"], baseline, effect_a1, effect_a2, effect_d,all_effect)
    pheno_par2 =  get_phenotype(geno_par["par2"], baseline, effect_a1, effect_a2, effect_d,all_effect)
    pheno_par3 =  get_phenotype(geno_par["par3"], baseline, effect_a1, effect_a2, effect_d,all_effect)
    pheno_par4 =  get_phenotype(geno_par["par4"], baseline, effect_a1, effect_a2, effect_d,all_effect)

    pheno_table=DataFrame(vcat(pheno_par1, pheno_par2, pheno_par3, pheno_par4),:auto)
    name_par = ["par1", "par2", "par3", "par4"]
    pheno_table = hcat(pheno_table, name_par,makeunique=true)
    rename!(pheno_table, [:baseline, :environment_val, :genotype_val, :phenotype_val, :name])
    sort!(pheno_table, :phenotype_val)

    small_parents = pheno_table[1:2, :]
    parentsmall_name1 = small_parents[1, :name]
    parentsmall_name2 = small_parents[2, :name]
    large_parents = pheno_table[3:4, :]
    parentlarge_name1 = large_parents[1, :name]
    parentlarge_name2 = large_parents[2, :name]

    geno_off = Dict(
        "off1_small" => get_offspring(geno_par[parentsmall_name1], geno_par[parentsmall_name2]),
        "off2_small" => get_offspring(geno_par[parentsmall_name1], geno_par[parentsmall_name2]),
        "off1_large" => get_offspring(geno_par[parentlarge_name1], geno_par[parentlarge_name2]),
        "off2_large" => get_offspring(geno_par[parentlarge_name1], geno_par[parentlarge_name2]))

    pheno_off1_small =  get_phenotype(geno_off["off1_small"], baseline, effect_a1, effect_a2, effect_d,all_effect)
    pheno_off2_small =  get_phenotype(geno_off["off2_small"], baseline, effect_a1, effect_a2, effect_d,all_effect)
    pheno_off1_large = get_phenotype(geno_off["off1_large"], baseline, effect_a1, effect_a2, effect_d,all_effect)
    pheno_off2_large =  get_phenotype(geno_off["off2_large"], baseline, effect_a1, effect_a2, effect_d,all_effect)

    all_data_df = DataFrame(
        simulation_n=simulation,
        group_n = group,
        largep_env1=large_parents[1,2],
        largep_geno1=large_parents[1,3],
        largep_pheno1=large_parents[1,4],
        largep_env2=large_parents[2,2],
        largep_geno2=large_parents[2,3],
        largep_pheno2=large_parents[2,4],
        midparent_large = mean(large_parents.phenotype_val),
        largeo_env1=pheno_off1_large[2],
        largeo_geno1=pheno_off1_large[3],
        largeo_pheno1=pheno_off1_large[4],
        largeo_env2=pheno_off2_large[2],
        largeo_geno2=pheno_off2_large[3],
        largeo_pheno2=pheno_off2_large[4],
        midoffspring_large = mean([pheno_off1_large[4] pheno_off2_large[4]]),
        smallp_env1=small_parents[1,2],
        smallp_geno1=small_parents[1,3],
        smallp_pheno1=small_parents[1,4],
        smallp_env2=small_parents[2,2],
        smallp_geno2=small_parents[2,3],
        smallp_pheno2=small_parents[2,4],
        midparent_small = mean(small_parents.phenotype_val),
        smallo_env1=pheno_off1_small[2],
        smallo_geno1=pheno_off1_small[3],
        smallo_pheno1=pheno_off1_small[4],
        smallo_env2=pheno_off2_small[2],
        smallo_geno2=pheno_off2_small[3],
        smallo_pheno2=pheno_off2_small[4],
        midoffspring_small = mean([pheno_off1_small[4], pheno_off2_small[4]]),
    )

    large_df = all_data_df[:,1:16]
    large_df = hcat(large_df, DataFrame(type = "large"))
    rename!(large_df,[:simulation,:group,:p_env1,:p_geno1,:p_pheno1,
    :p_env2,:p_geno2,:p_pheno2,
    :midpar,
    :o_env1,:o_geno1,:o_pheno1,
    :o_env2,:o_geno2,:o_pheno2,
    :midoff,:type])

    small_df = all_data_df[:,[1; 2; 17:30]]
    small_df = hcat(small_df, DataFrame(type = "small"))
    rename!(small_df,[:simulation,:group,:p_env1,:p_geno1,:p_pheno1,
    :p_env2,:p_geno2,:p_pheno2,
    :midpar,
    :o_env1,:o_geno1,:o_pheno1,
    :o_env2,:o_geno2,:o_pheno2,
    :midoff,:type])

    all_data_df2 = vcat(small_df,large_df)
    push!(all_data, all_data_df2)
end

# Concatenate into single dataframe:
all_g_df = vcat(all_data...)

mod1 = lm(@formula(midoff ~ midpar), all_g_df)
intercept = round(coef(mod1)[1],digits=3)
slope = round(coef(mod1)[2],digits=3)
r_squared = round(r2(mod1),digits=3)

push!(slope_vec, slope)
push!(intercept_vec, intercept)
push!(R_squared_vec, r_squared)

subtitle1 = "y~" * string(slope) * "*x + " * string(intercept) * ". R^2=" * string(r_squared)

@rput all_g_df
@rput subtitle1
R"""
p1_s <- ggplot(all_g_df,aes(x=midpar,y=midoff))+
  theme_minimal()+
  geom_point()+
  geom_abline()+
  ggtitle(subtitle1)
"""

p1_s =R"p1_s"

push!(p1, p1_s)

# Store list of dataframes:
push!(all_sim, all_g_df)

end

# Combine all data into single df:
all_sim_df = vcat(all_sim...)

mod2 = lm(@formula(midoff ~ midpar), all_sim_df)
intercept = round(coef(mod2)[1],digits=3)
slope = round(coef(mod2)[2],digits=3)
r_squared = round(r2(mod2),digits=3)

title2 = string(number_loci) * " loci, " * string(number_of_dom) * " dominant"
subtitle2 = "y~" * string(slope) * "*x + " * string(intercept) * ". R^2=" * string(r_squared)

# Initialize an empty DataFrame to store regression lines
regression_lines_df = DataFrame(x = Float64[], y = Float64[], simulation = String[])

# Loop through each slope and intercept combination
for i in 1:length(slope_vec)
    # Calculate the endpoints of the line using the range of x values
    x_range = extrema(all_sim_df.midpar)
    y_range = slope_vec[i] .* x_range .+ intercept_vec[i]

    xy_tuples = [(x_range[i], y_range[i]) for i in 1:2]
    df = DataFrame(x = [xy[1] for xy in xy_tuples], y = [xy[2] for xy in xy_tuples])
    df.simulation .= "Simulation $i"
   
    # Add this line to the DataFrame storing all regression lines
    append!(regression_lines_df, df)
end


# Plot the original points and regression lines
@rput all_sim_df
@rput regression_lines_df
@rput title2
@rput subtitle2
@rput minimum_pheno_value
R"""
p2<- ggplot(all_sim_df,aes(x = midpar, y = midoff)) +
  geom_point() +
  geom_line(data = regression_lines_df,
            aes(x = x, y = y, group = simulation),
            color = "red", alpha = 0.2) +
  geom_abline(slope=mean(slope_vec),
            intercept=mean(intercept_vec),
            col="blue",lwd=1.5)+
  xlab("Midparent diameter (cm)")+
  ylab("Midoffspring diameter (cm)")+
  theme_minimal() +
  labs(title=title2,
          subtitle=subtitle2)+
  geom_vline(aes(xintercept=10))+
  geom_hline(aes(yintercept=10))+
  geom_vline(aes(xintercept=minimum_pheno_value,col="red"))+
  geom_hline(aes(yintercept=minimum_pheno_value,col="red"))+
  theme(legend.position="none")
"""

p2 = R"p2"


  "Heritability less than 0=" * string(sum(slope_vec.<=0)/length(slope_vec))
  "Heritability less than 0.45=" * string(sum(slope_vec .<= 0.45) / length(slope_vec))
  "Average Heritability" * string(mean(slope_vec))

  df2= DataFrame(slope=slope_vec,intercept=intercept_vec)

@rput df2

R"""
scatter <- ggplot(df2,aes(x=slope,y=intercept))+
  theme_minimal()+
  geom_hex()

hist_bottom = ggplot()+
geom_histogram(data=df2,aes(slope))+
theme_minimal()
hist_right = ggplot()+
geom_histogram(data=df2,aes(intercept))+
coord_flip()+
theme_minimal()

empty <- ggplot() +
  geom_point(aes(1, 1), colour="white") +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
        
p3 <- grid.arrange(scatter, 
hist_right,
hist_bottom, 
empty,
ncol=2, nrow=2, widths=c(4, 1), heights=c(4, 1))
"""

my_plot = R"p3"