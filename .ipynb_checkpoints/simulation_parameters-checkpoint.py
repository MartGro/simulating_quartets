from simulate_tree import constant
import scipy


#################
#Things to try
################
possible_substitution_model = [
    "random",
    "wag",
]
possible_alpha_range = [(0.05, 1),(0.99, 1),(0.49, 0.5),(0.0049,0.005)]

possible_profile = ["all.freq","equi.freq"]

possible_profile_resampler = [('dirichlet', 10),('dirichlet', 10000)]

possible_heterogeneous_branch_ratio = ["random",0,1,0.5]

possible_rate_swap_ratio = [constant(0),constant(1),constant(0.1),"random"]

profile_swap_model = [scipy.stats.randint(0, 20),constant(0),constant(10),constant(20)]






##################################
#General settings
#######################

param_dict_template ={
    "name":"template",
    "taxon_count_model":constant(20),
    "internal_branch_model":scipy.stats.uniform(0.02, 1),
    "external_branch_model":scipy.stats.uniform(0.02, 1),
    "site_count_range":(200, 201),
    "substitution_model":"random",
    "alpha_range":(0.05, 1),
    "profile":'all.freq',
    "profile_resampler":('dirichlet', 10),
    "heterogeneous_branch_ratio":'random',
    "rate_swap_ratio":'random',
    "profile_swap_model":scipy.stats.randint(0, 20)

}

wag_dict = {"substitution_model":"wag"}

no_het_dict = {    
    "heterogeneous_branch_ratio":0.0,
    "rate_swap_ratio":constant(0),
    "profile_swap_model":constant(0)
    }

equi_profile_dict = {
    "profile":'equi.freq',
}


big_alpha_dict = {
    "alpha_range":(0.99, 1)
}


#####################################
#Actual model parameter sets
###################

##
training1_params = {
    "name":"training1_base_case",
    "taxon_count_model":constant(20),
    "internal_branch_model":scipy.stats.uniform(0.02, 1),
    "external_branch_model":scipy.stats.uniform(0.02, 1),
    "site_count_range":(200, 201),
    "substitution_model":"random",
    "alpha_range":(0.05, 1),
    "profile":"all.freq",
    "profile_resampler":('dirichlet', 10),
    "heterogeneous_branch_ratio":0.9,
    "rate_swap_ratio":'random',
    "profile_swap_model":scipy.stats.randint(10, 20),
}



##
training1_wag_params = training1_params.copy()
training1_wag_params.update(wag_dict)
training1_wag_params.update({
    "name":"training1_wag",
})

##
training1_wag_no_het_params = training1_params.copy()
training1_wag_no_het_params.update({
    "name":"training1_wag_no_het",})
training1_wag_no_het_params.update(wag_dict)
training1_wag_no_het_params.update(no_het_dict)


##
training1_equi_profile = training1_params.copy()
training1_equi_profile.update(equi_profile_dict)
training1_equi_profile.update({
    "name":"training1_equi_profile",})

#######################
#######################

training2_params = {
    "name":"training2_base_case",
    "taxon_count_model":constant(20),
    "internal_branch_model":scipy.stats.gamma(a=1, scale=0.1/1),
    "external_branch_model":scipy.stats.gamma(a=2, scale=0.2/2),
    "site_count_range":(200, 201),
    "substitution_model":"random",
    "alpha_range":(0.05, 1),
    "profile":"all.freq",
    "profile_resampler":('dirichlet', 10),
    "heterogeneous_branch_ratio":0.9,
    "rate_swap_ratio":'random',
    "profile_swap_model":scipy.stats.randint(10, 20),
}



#######################
#######################


training3_params = {
    "name":"training3_base_case",
    "taxon_count_model":constant(20),
    "internal_branch_model":scipy.stats.uniform(0.02, 1),
    "external_branch_model":scipy.stats.uniform(0.02, 1),
    "site_count_range":(200, 201),
    "substitution_model":"random",
    "alpha_range":(0.05, 1),
    "profile":"all.freq",
    "profile_resampler":('dirichlet', 10),
    "heterogeneous_branch_ratio":'random',
    "rate_swap_ratio":'random',
    "profile_swap_model":scipy.stats.randint(0, 20),
}