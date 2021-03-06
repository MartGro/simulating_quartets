import ete3
import evosimz.tree as tree
import evosimz.sequence as sequence
import scipy
import subprocess
import os
from ete3 import Tree
from phylogeny_utilities import utilities
import pickle
import itertools
import generate_quartets


__all__ = [
    'constant',
    'get_pdist_from_fasta',
    'create_phylip_string',
    'create_fasta_string',
    'create_parameter_string',
    'write_fasta_file',
    'write_phylip_file',
    'write_original_newick_string',
    'write_parameter_file',
    'run_iqtree',
    'read_iq_tree_newick',
    'make_simulator_from_parameters',
    'print_param_dict',
    'simulate_tree_from_parameters',
    'get_iqtree_accuracy',
    
    
    
    
]


class constant():
    def __init__(self,constant):
        self.constant = constant
        self.args = constant
    def rvs(self,*args):
        return self.constant
    def __str__(self):
        return str(self.constant)
    
    
    
#Input
#:fasta_filepath: path to a fasta file
#Output
#Dictionary with average_p_distance, maximum_p_distance, number_of_pairs(# of sequences choose 2), number of columns that never change, total length of the alignment
def get_pdist_from_fasta(fasta_filepath):
    #get average pdistance
    avg_p,equal_columns,tot_columns,num_pairs = utilities.get_avg_pdistance_of_fasta(fasta_filepath)    
    #Get the maximum p distance
    max_p= utilities.get_avg_pdistance_of_fasta(fasta_filepath,True)
    
    pdist_dict = {"avg_p":avg_p,"max_p":max_p,"num_pairs":num_pairs,"equal_columns":equal_columns,"tot_columns":tot_columns}
    return pdist_dict




#create a .phylip formatted string for 20 sequences and a lenght of 200
def create_phylip_string(simulated_tree):
    return_str = " 20 200\n"
    for leaf in simulated_tree:
        return_str+="{}    {}\n".format(leaf.name,leaf.sequence.decode())
    return return_str


def create_fasta_string(simulated_tree):
    return_str = ""
    for leaf in simulated_tree:
        return_str+=">{}\n{}\n".format(leaf.name,leaf.sequence.decode())
    return return_str

def create_parameter_string(parameters):
    s = ""
    s+= "name:                        {}\n".format(parameters["name"])
    s+= "substitution_model:          {}\n".format(parameters["substitution_model"])
    s+= "alpha_range:                 {}\n".format(parameters["alpha_range"])
    s+= "profile:                     {}\n".format(parameters["profile"])
    s+= "profile_resampler:           {}\n".format(parameters["profile_resampler"])
    s+= "heterogeneous_branch_ratio:  {}\n".format(parameters["heterogeneous_branch_ratio"])
    s+= "rate_swap_ratio:             {}\n".format(str(parameters["rate_swap_ratio"]))
    s+= "profile_swap_model:          {}         #(number of swaps)\n\n\n".format(str(parameters["profile_swap_model"].args))
    return s



def write_fasta_file(dir_path,file_name,fasta_str):
    file_path = dir_path+file_name+".fasta"
    with open(file_path,"w") as f:
        f.write(fasta_str)
    return file_path

def write_phylip_file(dir_path,file_name,phylip_str):
    file_path = dir_path+file_name+".phy"
    with open(file_path,"w") as f:
        f.write(phylip_str)

def write_original_newick_string(dir_path,file_name,newick_string):
    file_path = dir_path+file_name+".orig_tree"
    with open(file_path,"w") as f:
        f.write(newick_string)
        
def write_parameter_file(dir_path,file_name,parameter_dict):
    param_string = create_parameter_string(parameter_dict)
    file_path = dir_path+file_name+".params"
    with open(file_path,"w") as f:
        f.write(param_string)
        
        
def dump_results(dir_path,file_name,result_dict):
    file_path = dir_path+file_name+".result_pickle"
    with open(file_path,"wb") as f:
        pickle.dump(result_dict,f)


#This depends on the environment -> could be made independent
def run_iqtree(fasta_filepath,model=None,specify_cores=None):
    parameters = ["../iqtree-2.0.5-Linux/bin/iqtree2", "-s",fasta_filepath,"-redo"]
    if model != None:
        parameters += ["-m",model]
    if specify_cores != None:
        parameters += ["-nt",specify_cores]
    subprocess.call(parameters)

        
def read_iq_tree_newick(fasta_filepath):
    tree_filename = fasta_filepath+".treefile"
    with open(tree_filename,"r") as f:
        tree_newick = f.read()
    return tree_newick

def make_simulator_from_parameters(parameters):
    simulator = tree.TreeSimulator(
        taxon_count_model=parameters["taxon_count_model"],
        internal_branch_model=parameters["internal_branch_model"],
        external_branch_model=parameters["external_branch_model"],
        sequence_simulator=sequence.HeterogeneousProteinSequenceSimulator(
            site_count_range=parameters["site_count_range"],
            substitution_model=parameters["substitution_model"],
            alpha_range=parameters["alpha_range"],
            profile=parameters["profile"],
            profile_resampler=parameters["profile_resampler"],
            heterogeneous_branch_ratio=parameters["heterogeneous_branch_ratio"],
            rate_swap_ratio=parameters["rate_swap_ratio"],
            profile_swap_model=parameters["profile_swap_model"],
        ))
    return simulator



#def prune_quartets_from_tree(simulated_tree):
#    sequence_list = []
#    for leaf in simulated_tree:
#        sequence_list.append((leaf.name,leaf.sequence.decode())
#                             
#    for combination in itertools.combinations(sequence_list,4):
#        tree_copy = simulated_tree.copy()
#        tree_
#    return sequence_list
                             





def simulate_tree_from_parameters(parameter_dict,dir_path,file_name):
    simulator=make_simulator_from_parameters(parameter_dict)
    simulated_tree = simulator.generate()
    
    fasta_str = create_fasta_string(simulated_tree)
    fasta_filepath = write_fasta_file(dir_path,file_name,fasta_str)
    
    original_newick = simulated_tree.write()
    write_original_newick_string(dir_path,file_name,original_newick)
    
    phylip_str = create_phylip_string(simulated_tree)
    write_phylip_file(dir_path,file_name,phylip_str)
    
    write_parameter_file(dir_path,file_name,parameter_dict)
    
    #p-distance evaluation file from evaluation data
    pdist_dict = get_pdist_from_fasta(fasta_filepath)
    return (pdist_dict,fasta_filepath,simulated_tree)



def create_quartet_pickles(simulated_tree,dir_path):
    quartet_sampler = generate_quartets.given_tree_quartet_sampler_from_full_tree(simulated_tree)
    dir_path = dir_path+"/quartet_files"
    os.makedirs(dir_path)
    generate_quartets.create_pickle_from_quartet_sampler(quartet_sampler,dir_path)



def get_iqtree_accuracy(parameter_dict,dir_path,file_name,model = None,specify_cores = None):
    #simulate a tree
    pdist_dict,fasta_filepath,simulated_tree = simulate_tree_from_parameters(parameter_dict,dir_path,file_name)
    create_quartet_pickles(simulated_tree,dir_path)
    
    run_iqtree(fasta_filepath,model,specify_cores)
    
    iq_tree_newick = read_iq_tree_newick(fasta_filepath)
    iq_tree_tree = Tree(iq_tree_newick)    
    
    comparison = simulated_tree.compare(iq_tree_tree,unrooted=True)
    avg_p = pdist_dict["avg_p"]
    max_p = pdist_dict["max_p"]
    return_dict = {"rf":comparison["rf"],"size":comparison["effective_tree_size"],"max_rf":comparison["max_rf"],"avg_p":avg_p,"max_p":max_p}
    dump_results(dir_path,file_name,return_dict)
    return return_dict
    
    

