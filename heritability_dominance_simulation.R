library(tidyverse)
library(broom)
library(patchwork)
library(gridExtra)


num_simulations <- 100
num_groups <- 30

number_loci <- 6
min_num_hetero <- 3 # Minimum number of loci that are heterozygotic

number_of_dom <- 0 # Number of dominant loci
effect_A1 <- +1 # Value +1
effect_A2 <- 0 # Value 0

# Additive d=0, complete dominance d=a1, partial dominance mid<d<a1,
# Overdominance d>a1, underdominance. d<a1
effect_d <- 0.5

env_type <- "Uniform"
env_sd <- 2
env_effect <- 1 # Multiplier of the env_effect
epistasis <- TRUE # Turn to TRUE, to activate epistasis
epistasis_level <- 6 # If N loci have at least 1 A1
epistasis_value <- +10 # add X to the genotypic value
baseline <- 10


# Uniform distribution, with expected value 0
#env_range <-  # Normal distribution

if (env_type=="Uniform"){
  env_range <- c(-3:-1,1,3)
} else if (env_type=="Normal"){
  env_range <- round(rnorm(n = 1000,mean=0,sd=env_sd),2)
} else{print("error")}


possible_num_homo <- 0:(number_loci-min_num_hetero)
number_of_add <- number_loci - number_of_dom # Number of additive loci
add_vec <- rep("ADD",number_of_add)
dom_vec <- rep("DOM",number_of_dom)
all_effect <- c(add_vec,dom_vec)
all_effect <- sample(all_effect,size=number_loci,replace = FALSE)
mid_effect <- mean(c(effect_A1,effect_A2)) # Midpoint value
effect_a1 <- effect_A1-mid_effect # Value +a
effect_a2 <- effect_A2-mid_effect # Value -a
minimum_pheno_value <- baseline+env_effect*min(env_range)+2*number_loci*min(c(effect_A1,effect_A2))

slope_vec <- c()
intercept_vec <- c()
R.squared_vec <- c()
plot_vec <- list()
all_sim <- list()
for (simulation in 1:num_simulations){
all_data <- list()

for (group in 1:num_groups){

sample_genotype <- function(num_loci){
  hom <- sample(possible_num_homo,1) # How many loci are homozygote per genotype
  het <- num_loci-hom # How many loci are heterozygote per genotype
  hom_vec <- rep("HOM",hom)
  het_vec <- rep("HET",het)
  
  all_loci <- c(hom_vec,het_vec)
  all_loci <- sample(all_loci,size=num_loci,replace = FALSE) #vector with homo/heterozygote values
  chromosome1 <- sample(c("A1","A2"),num_loci,replace=TRUE) # random values for A1 or A2 in one chromosome
  chromosome2 <- c()
  for (i in seq_along(all_loci)){ # For every loci:
    if (all_loci[i]=="HOM"){ # If value is homozygote
      chromosome2[i] <- chromosome1[i] # The second chromosome has the same value
    } else if(all_loci[i]=="HET"){ # If value is heterozygote
      if (chromosome1[i]=="A1"){ # and the chromosome 1 is A1
        chromosome2[i] <- "A2" # chromosome 2 is A2
      } else if (chromosome1[i]=="A2"){ # and the chromosome 1 is A2
        chromosome2[i] <- "A1" # chromosome 2 is A1
      } 
    else{print("error")}
    }
    
  genotype <- rbind(chromosome1,chromosome2) # Bind both chromosomes to form a genotype
  return(genotype)
}
}

# Sample 4 parental genotypes:
par1 <- sample_genotype(number_loci) 
par2 <- sample_genotype(number_loci)
par3 <- sample_genotype(number_loci)
par4 <- sample_genotype(number_loci)
all_geno <- list(par1,par2,par3,par4)
names(all_geno) <- c("par1","par2","par3","par4")

# Given two parental genotypes, generate 2 gametes and combine them
get_offspring <- function(genotype1,genotype2){
  gamete1 <- c()
  for (i in 1:dim(genotype1)[2]){
    if (genotype1[1,i]==genotype1[2,i]){
      # If locus is homozygote, position is the same 
      gamete1[i] <- genotype1[1,i]
    } else if (genotype1[1,i]!=genotype1[2,i]){
      # If locus is heterozygote, roll a coin 
      gamete1[i] <- sample(c(genotype1[1,i],genotype1[2,i]),1)
    }
  }
  
  gamete2 <- c()
  for (i in 1:dim(genotype2)[2]){
    if (genotype2[1,i]==genotype2[2,i]){
      # If locus is homozygote, position is the same 
      gamete2[i] <- genotype2[1,i]
    } else if (genotype2[1,i]!=genotype2[2,i]){
      # If locus is heterozygote, roll a coin 
      gamete2[i] <- sample(c(genotype2[1,i],genotype2[2,i]),1)
    }
  }
  
  # Combine both gametes to make a new genotype
  new_genotype <- rbind(gamete1,gamete2)
  return(new_genotype)
}


# Get phenotype given a genotype, additive, dominance, environment and epistasis.
get_phenotype <- function(geno,baseline=baseline,
                          effect_a1,effect_a2,effect_d){
  effect_val <- c()
  for (e in seq_along(all_effect)){ # For every loci
    if (all_effect[e]=="DOM"){ # If the locus is dominant:
      effect_val[e] <- ifelse(par1[1,e]=="A1"&par2[2,e]=="A1",
                              2*(mid_effect+effect_a1),
                                  ifelse(par1[1,e]=="A1"&par2[2,e]=="A2"|
                                           par1[1,e]=="A2"&par2[2,e]=="A1",
                                         mid_effect+effect_d+mid_effect+effect_d,
                                         ifelse(par1[1,e]=="A2"&par2[2,e]=="A2",
                                                2*(mid_effect+effect_a2))))
      
      
    } else if(all_effect[e]=="ADD"){ # If the locus is additive:
      effect_val[e] <- ifelse(par1[1,e]=="A1"&par2[2,e]=="A1",
                              2*(mid_effect+effect_a1),
                              ifelse(par1[1,e]=="A1"&par2[2,e]=="A2"|
                                       par1[1,e]=="A2"&par2[2,e]=="A1",
                                     mid_effect+effect_a1+mid_effect+effect_a2,
                                     ifelse(par1[1,e]=="A2"&par2[2,e]=="A2",
                                            2*(mid_effect+effect_a2))))
    }  
      else{print("error")}
  }
  
  # Epistasis term:
  if(epistasis==TRUE){ #If epistasis is present
    epis_val <- ifelse(sum(effect_val>=1)>=epistasis_level,TRUE,FALSE)
    geno_val <- ifelse(epis_val,sum(effect_val)+epistasis_value,sum(effect_val))
  } else{
    # No epistasis:
    geno_val <- sum(effect_val)
  }
  

  env_val <- sample(env_range,size=1) # Sample the environmental range
  pheno <- baseline+env_effect*env_val+geno_val
  return(c(baseline,env_val,geno_val,pheno))
}

pheno_par1 <- get_phenotype(par1,baseline=baseline,
                            effect_a1=effect_a1,
                            effect_a2=effect_a2,
                            effect_d=effect_d)
pheno_par2 <- get_phenotype(par2,baseline=baseline,
                            effect_a1=effect_a1,
                            effect_a2=effect_a2,
                            effect_d=effect_d)
pheno_par3 <- get_phenotype(par3,baseline=baseline,
                            effect_a1=effect_a1,
                            effect_a2=effect_a2,
                            effect_d=effect_d)
pheno_par4 <- get_phenotype(par4,baseline=baseline,
                            effect_a1=effect_a1,
                            effect_a2=effect_a2,
                            effect_d=effect_d)

pheno_table <- rbind(pheno_par1,pheno_par2,
      pheno_par3,pheno_par4) 
name <- gsub("pheno_","",rownames(pheno_table))
pheno_table <- pheno_table |> as_tibble()
pheno_table$name <- name
colnames(pheno_table) <- c("baseline","environment_val",
                           "genotype_val","phenotype_val","name")

all_parents_phenotypes_sorted <- pheno_table |> 
  arrange(phenotype_val)

small_parents <- all_parents_phenotypes_sorted[1:2,]
large_parents <- all_parents_phenotypes_sorted[3:4,]


offspring1_large <- get_offspring(get(pull(large_parents[1,5])),get(pull(large_parents[2,5])))
offspring2_large <- get_offspring(get(pull(large_parents[1,5])),get(pull(large_parents[2,5])))
pheno_large_off1 <- get_phenotype(offspring1_large,baseline=baseline,
                                  effect_a1=effect_a1,
                                  effect_a2=effect_a2,
                                  effect_d=effect_d)
pheno_large_off2 <- get_phenotype(offspring2_large,baseline=baseline,
                                  effect_a1=effect_a1,
                                  effect_a2=effect_a2,
                                  effect_d=effect_d)
offspring1_small <- get_offspring(get(pull(small_parents[1,5])),get(pull(small_parents[2,5])))
offspring2_small <- get_offspring(get(pull(small_parents[1,5])),get(pull(small_parents[2,5])))
pheno_small_off1 <- get_phenotype(offspring1_small,baseline=baseline,
                                  effect_a1=effect_a1,
                                  effect_a2=effect_a2,
                                  effect_d=effect_d)
pheno_small_off2 <- get_phenotype(offspring2_small,baseline=baseline,
                                  effect_a1=effect_a1,
                                  effect_a2=effect_a2,
                                  effect_d=effect_d)
midparent_large <- mean(large_parents$phenotype_val)
midparent_small <- mean(small_parents$phenotype_val)

midoffspring_large <- mean(pheno_large_off1[4],pheno_large_off2[4])
midoffspring_small <- mean(pheno_small_off1[4],pheno_small_off2[4])

all_data[[group]] <- c(group,unlist(large_parents[1,c(-1,-5)]),
                     unlist(large_parents[2,c(-1,-5)]),
                     midparent_large,
                     pheno_large_off1[-1],
                     pheno_large_off2[-1],
                     midoffspring_large,
                     unlist(small_parents[1,c(-1,-5)]),
                     unlist(small_parents[2,c(-1,-5)]),
                     midparent_small,
                     pheno_small_off1[-1],
                     pheno_small_off2[-1],
                     midoffspring_small)
}

all_data_df <- do.call("rbind",all_data)
colnames(all_data_df) <- c("group","largep_env1","largep_geno1","largep_pheno1",
                           "largep_env2","largep_geno2","largep_pheno2",
                           "large_midpar",
                           "largeo_env1","largeo_geno1","largeo_pheno1",
                           "largeo_env2","largeo_geno2","largeo_pheno2",
                           "large_moff",
                           "smallp_env1","smallp_geno1","smallp_pheno1",
                           "smallp_env2","smallp_geno2","smallp_pheno2",
                           "small_midpar",
                           "smallo_env1","smallo_geno1","smallo_pheno1",
                           "smallo_env2","smallo_geno2","smallo_pheno2",
                           "small_midoff")

large_df <- all_data_df[,1:15] |> as.data.frame()
large_df$type <- "large"
colnames(large_df) <- c("group","p_env1","p_geno1","p_pheno1",
                        "p_env2","p_geno2","p_pheno2",
                        "midpar",
                        "o_env1","o_geno1","o_pheno1",
                        "o_env2","o_geno2","o_pheno2",
                        "midoff","type")
small_df <- all_data_df[,c(1,16:29)]|> as.data.frame()
small_df$type <- "small"
colnames(small_df) <- c("group","p_env1","p_geno1","p_pheno1",
                        "p_env2","p_geno2","p_pheno2",
                        "midpar",
                        "o_env1","o_geno1","o_pheno1",
                        "o_env2","o_geno2","o_pheno2",
                        "midoff","type")

all_data_df2 <- rbind(large_df,
                      small_df)

all_sim[[simulation]] <- all_data_df2
mod1 <- lm(midoff~midpar,all_data_df2)
glance(mod1)

slope <- pull(round(tidy(mod1)[2,2],3))
intercept <- pull(round(tidy(mod1)[1,2],3))
r.squared <- pull(round(glance(mod1)[1],3))


slope_vec[[simulation]] <- slope
intercept_vec[[simulation]] <- intercept
R.squared_vec[[simulation]] <- r.squared

subtitle1 <- paste0("y~",
       slope,
       "*x + ",
       intercept,
       ". R^2=",
       r.squared)
       

plot_vec[[simulation]] <- all_data_df2 |> 
  ggplot(aes(x=midpar,y=midoff))+
  theme_minimal()+
  geom_point()+
  geom_abline()+
  ggtitle(subtitle1)

}

# Data for 1000 simulations:

slope_vec <- unlist(slope_vec)
intercept_vec <- unlist(intercept_vec)
R.squared_vec <- unlist(R.squared_vec)
names(all_sim) <- paste0("Simulation",1:length(all_sim))
df <- bind_rows(all_sim, .id = "simulation_name")

mod2 <- lm(midoff~midpar,df)
glance(mod2)

slope <- pull(round(tidy(mod2)[2,2],3))
intercept <- pull(round(tidy(mod2)[1,2],3))
r.squared <- pull(round(glance(mod2)[1],3))

title1 <- paste0(number_loci," loci, ",number_of_dom, " dominant")

subtitle1 <- paste0("y~",
                 slope,
                 "*x + ",
                 intercept,
                 ". R^2=",
                 r.squared)

df |> 
  ggplot(aes(x=midpar,y=midoff))+
  geom_point()+
  geom_abline()+
  theme_minimal()+
  labs(title=title1,
          subtitle=subtitle1)


regression_lines_df <- data.frame()

# Loop through each slope and intercept combination
for (i in 1:length(slope_vec)) {
  # Calculate the endpoints of the line using the range of x values
  x_range <- range(df$midpar)
  y_range <- slope_vec[i] * x_range + intercept_vec[i]
  
  # Create a dataframe for this line
  line_df <- data.frame(x = x_range, y = y_range)
  
  # Add a column to identify the simulation
  line_df$simulation <- paste0("Simulation ", i)
  
  # Add this line to the dataframe storing all regression lines
  regression_lines_df <- rbind(regression_lines_df, line_df)
}

# Plot the original points and regression lines
p1 <- df |> 
  ggplot(aes(x = midpar, y = midoff)) +
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
  labs(title=title1,
          subtitle=subtitle1)+
  geom_vline(aes(xintercept=10))+
  geom_hline(aes(yintercept=10))+
  geom_vline(aes(xintercept=minimum_pheno_value,col="red"))+
  geom_hline(aes(yintercept=minimum_pheno_value,col="red"))+
  theme(legend.position = "none")
  

sum(slope_vec<=0.45)/length(slope_vec)

print(paste0("Heritability less than 0=",sum(slope_vec<=0)/length(slope_vec)))

df2 <- cbind(slope=slope_vec,intercept=intercept_vec) |> 
  as.data.frame()

scatter <- df2 |> 
  ggplot(aes(x=slope,y=intercept))+
  theme_minimal()+
  geom_hex()

hist_bottom <- ggplot()+geom_histogram(aes(df2$slope))+theme_minimal()
hist_right <- ggplot()+geom_histogram(aes(df2$intercept))+coord_flip()+theme_minimal()
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

p2 <- grid.arrange(scatter, 
             hist_right,
             hist_bottom, 
             empty,
             ncol=2, nrow=2, widths=c(4, 1), heights=c(4, 1))


