import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

num_simulations = 100
num_groups = 30
num_loci = 6
min_num_hetero = 3
number_of_dom = 0
effect_A1 = 1
effect_A2 = 0
effect_d = 0.5
env_type = "Uniform"
env_sd = 2
env_effect = 1
epistasis = False
epistasis_level = 6
epistasis_value = 10
baseline = 10