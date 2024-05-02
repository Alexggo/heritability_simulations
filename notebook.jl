### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 3b9e8fbc-06a6-11ef-3a3b-3d7aa4d6dcbe
begin
	using DataFrames
	using Random
	using RCall
	using Statistics
	using PlutoUI
end

# ╔═╡ a8104d0c-d0f8-465f-b008-b3f88ee566a4
@bind num_simulations Slider(1:300,20,true)

# ╔═╡ 05dfbf5f-0dc4-417f-b3e0-3aee1f664c9c
@bind num_groups Slider(1:30,20,true)

# ╔═╡ f99bb201-0e81-479f-9351-08338eee3aed
@bind min_num_hetero Slider(0:6,3,true)

# ╔═╡ cbb31dd4-2f5c-4b80-9b0c-acc5c9907592
@bind number_of_dom Slider(0:6,2,true)

# ╔═╡ abf49dab-2d3b-4d77-8343-cc823a4328a8
@bind effect_A1 Slider(-3:3,1,true)

# ╔═╡ ad4c0998-c9ea-434d-8e5d-3cb11dde3319
@bind effect_A2 Slider(-3:3,1,true)

# ╔═╡ b531fd9e-60e0-432e-83ca-e6e1e7ced580
@bind effect_d Slider(-2:2,0,true)

# ╔═╡ fe9a34ef-40f3-40cc-9dd7-4e760d83f389
@bind baseline Slider(-20:20,10,true)

# ╔═╡ 2228e6de-c06f-4ca8-b25b-f5bcb782da5c
@bind env_type Select(["Normal","Uniform"])

# ╔═╡ 5a354555-79f7-4286-afbb-49ced32a534b
@bind env_effect Slider(-3:3,1,true)

# ╔═╡ 87aad4db-c9d8-413e-83ba-5958507992a4
@bind env_sd Slider(-3.0:3.0,2.0,true)

# ╔═╡ 7e30da1e-a237-4f8b-b699-1bec2d3d04dc
@bind epistasis Select([true,false])

# ╔═╡ 64f83cb7-4b3a-499c-ade2-aaa0fb595dfd
@bind epistasis_level Slider(0.0:6.0,3.0,true)

# ╔═╡ bdcced23-2c32-4346-9d3d-1d44828d362f
@bind epistasis_value Slider(-10.0:10.0,10.0,true)

# ╔═╡ af08ac1d-c67e-41db-be93-131da3d1a63b
number_loci = 6

# ╔═╡ dca5dfaa-3ed9-4164-a0c3-d20a29c93763
possible_num_homo = 0:(number_loci - min_num_hetero)

# ╔═╡ 5452a1bb-2357-46d6-b8a9-965d28c5164a
number_of_add = number_loci - number_of_dom

# ╔═╡ debb787f-8dd5-4334-bd98-0fbf4d361c0c
add_vec = fill("ADD", number_of_add)

# ╔═╡ a8899725-c94d-430a-a35e-281095b27f70
dom_vec = fill("DOM", number_of_dom)

# ╔═╡ 25fb7e46-3b15-4814-8fed-ffa829b1adfb
all_effect1 = vcat(add_vec, dom_vec)

# ╔═╡ 447e8a57-3327-45ed-85af-0d0df2e0db1e
all_effect = shuffle(all_effect1)

# ╔═╡ 39d737e9-151a-499c-9b39-ac46d9096777
mid_effect = mean([effect_A1, effect_A2])

# ╔═╡ e1322b7f-7e76-4d5f-88ba-e308909ff1c2
effect_a1 = effect_A1 - mid_effect

# ╔═╡ e8c19eec-ffe9-47bf-a488-724bd4b1ec3d
effect_a2 = effect_A2 - mid_effect

# ╔═╡ b5c4ddb1-2e0e-47a9-bb34-8ea9fd35a0ad

if env_type == "Uniform"
  env_range = vcat(-3:-1, 1:3)
elseif env_type == "Normal"
  env_range = round.(randn(1000) .* env_sd, digits=2)
else
  println("error")
end

# ╔═╡ 758e05fa-5e7d-4cc2-97d0-c25862af4f13
minimum_pheno_value = baseline + env_effect * minimum(env_range) + 2 * number_loci * minimum([effect_A1, effect_A2])

# ╔═╡ 941ee4e4-686f-4ab1-b73a-5aa71ba9da13
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

# ╔═╡ 52db1c95-204a-4f7b-8d6e-2890ec9bde32
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
  return [baseline env_val geno_val pheno]
end

# ╔═╡ e17d9637-3173-4d35-b041-60c0a6549739
function get_offspring(genotype1, genotype2)
  gamete1 = []
  for i in 1:size(genotype1, 2)
      if genotype1[1, i] == genotype1[2, i]
          gamete1[i] = genotype1[1,i]
      elseif genotype1[1, i] != genotype1[2, i]
          gamete1[i] = rand([genotype1[1, i], genotype1[2, i]],1)
      end
  end

  gamete2 = []
  for i in 1:size(genotype2, 2)
      if genotype2[1, i] == genotype2[2, i]
          gamete2[i] = genotype2[1,i]
      elseif genotype2[1, i] != genotype2[2, i]
          gamete2[i] = rand([genotype2[1, i], genotype2[2, i]],1)
      end
  end

  new_genotype = [gamete1 gamete2]
  return new_genotype
end

# ╔═╡ 082686db-6582-432c-8af6-990619908b4a
# ╠═╡ disabled = true
#=╠═╡
i=2
  ╠═╡ =#

# ╔═╡ 16c01522-189e-43bd-b2fe-343ee4880793
  all_geno = Dict(
    "par1" => sample_genotype(number_loci),
    "par2" => sample_genotype(number_loci),
    "par3" => sample_genotype(number_loci),
    "par4" => sample_genotype(number_loci))

# ╔═╡ c84f4266-15bc-4a90-bc88-6ac4aff30999
genotype1=all_geno["par1"]

# ╔═╡ cd16d9b0-5af7-4831-9baf-17e85706e078
genotype2=all_geno["par2"]

# ╔═╡ 2a797e95-867e-49a3-8193-a47481cc90ec
get_offspring(all_geno["par1"],all_geno["par2"])

# ╔═╡ 64322ab5-9e1d-411a-8a9c-eb4a32e649b9
begin
	pheno_par1 = get_phenotype(all_geno["par1"], baseline, effect_a1, effect_a2, effect_d)
	pheno_par2 =  get_phenotype(all_geno["par2"], baseline, effect_a1, effect_a2, effect_d)
	pheno_par3 =  get_phenotype(all_geno["par3"], baseline, effect_a1, effect_a2, effect_d)
	pheno_par4 =  get_phenotype(all_geno["par4"], baseline, effect_a1, effect_a2, effect_d)
end

# ╔═╡ f749e5cb-61b4-4aac-ae76-66f73fa7b686
pheno_table=DataFrame(vcat(pheno_par1, pheno_par2, pheno_par3, pheno_par4),:auto)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DataFrames = "~1.6.1"
PlutoUI = "~0.7.59"
RCall = "~0.14.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "2b51acb1fe9bec812d6699ec0be4d92c9a9ab25a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "1568b28f91293458345dabba6a5ea3f183250a61"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.8"

    [deps.CategoricalArrays.extensions]
    CategoricalArraysJSONExt = "JSON"
    CategoricalArraysRecipesBaseExt = "RecipesBase"
    CategoricalArraysSentinelArraysExt = "SentinelArrays"
    CategoricalArraysStructTypesExt = "StructTypes"

    [deps.CategoricalArrays.weakdeps]
    JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SentinelArrays = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
    StructTypes = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "51cab8e982c5b598eea9c8ceaced4b58d9dd37c9"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.10.0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.RCall]]
deps = ["CategoricalArrays", "Conda", "DataFrames", "DataStructures", "Dates", "Libdl", "Missings", "Preferences", "REPL", "Random", "Requires", "StatsModels", "WinReg"]
git-tree-sha1 = "846b2aab2d312fda5e7b099fc217c661e8fae27e"
uuid = "6f49c342-dc21-5d91-9882-a32aef131414"
version = "0.14.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "5cf6c4583533ee38639f73b880f35fc85f2941e0"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.3"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WinReg]]
git-tree-sha1 = "cd910906b099402bcc50b3eafa9634244e5ec83b"
uuid = "1b915085-20d7-51cf-bf83-8f477d6f5128"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═3b9e8fbc-06a6-11ef-3a3b-3d7aa4d6dcbe
# ╠═a8104d0c-d0f8-465f-b008-b3f88ee566a4
# ╠═05dfbf5f-0dc4-417f-b3e0-3aee1f664c9c
# ╠═f99bb201-0e81-479f-9351-08338eee3aed
# ╠═cbb31dd4-2f5c-4b80-9b0c-acc5c9907592
# ╠═abf49dab-2d3b-4d77-8343-cc823a4328a8
# ╠═ad4c0998-c9ea-434d-8e5d-3cb11dde3319
# ╠═b531fd9e-60e0-432e-83ca-e6e1e7ced580
# ╠═fe9a34ef-40f3-40cc-9dd7-4e760d83f389
# ╠═2228e6de-c06f-4ca8-b25b-f5bcb782da5c
# ╠═5a354555-79f7-4286-afbb-49ced32a534b
# ╠═87aad4db-c9d8-413e-83ba-5958507992a4
# ╠═7e30da1e-a237-4f8b-b699-1bec2d3d04dc
# ╠═64f83cb7-4b3a-499c-ade2-aaa0fb595dfd
# ╠═bdcced23-2c32-4346-9d3d-1d44828d362f
# ╠═af08ac1d-c67e-41db-be93-131da3d1a63b
# ╠═dca5dfaa-3ed9-4164-a0c3-d20a29c93763
# ╠═5452a1bb-2357-46d6-b8a9-965d28c5164a
# ╠═debb787f-8dd5-4334-bd98-0fbf4d361c0c
# ╠═a8899725-c94d-430a-a35e-281095b27f70
# ╠═25fb7e46-3b15-4814-8fed-ffa829b1adfb
# ╠═447e8a57-3327-45ed-85af-0d0df2e0db1e
# ╠═39d737e9-151a-499c-9b39-ac46d9096777
# ╠═e1322b7f-7e76-4d5f-88ba-e308909ff1c2
# ╠═e8c19eec-ffe9-47bf-a488-724bd4b1ec3d
# ╠═b5c4ddb1-2e0e-47a9-bb34-8ea9fd35a0ad
# ╠═758e05fa-5e7d-4cc2-97d0-c25862af4f13
# ╟─941ee4e4-686f-4ab1-b73a-5aa71ba9da13
# ╟─52db1c95-204a-4f7b-8d6e-2890ec9bde32
# ╠═e17d9637-3173-4d35-b041-60c0a6549739
# ╠═c84f4266-15bc-4a90-bc88-6ac4aff30999
# ╠═cd16d9b0-5af7-4831-9baf-17e85706e078
# ╠═082686db-6582-432c-8af6-990619908b4a
# ╠═2a797e95-867e-49a3-8193-a47481cc90ec
# ╠═16c01522-189e-43bd-b2fe-343ee4880793
# ╠═64322ab5-9e1d-411a-8a9c-eb4a32e649b9
# ╠═f749e5cb-61b4-4aac-ae76-66f73fa7b686
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
