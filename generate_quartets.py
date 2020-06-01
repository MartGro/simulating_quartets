from evosimz.tree import BaseTreeSimulator
import evosimz
import random
import re
import pickle
import itertools
import math

import pathlib


path_fasta_mapping = pathlib.Path("/insert/path/to/fastafile")
path_pickles = "/home/martin/LRZ_Sync_Share/Spring20/Comp_Phylo/Project/phydl/data/trees/RFP_Pruned_All"

#####
#
#####
def mapping_from_fasta(path):
    with open(path,"r") as file: 
        content = file.read()
        
    mapping_dict = {}
    
    splits = content.split("\n")
    splits = [i.rstrip() for i in splits if len(i)>0]
    number_of_sequences = int(len(splits)/2)
    assert number_of_sequences*2 == len(splits)
    for i in range(0,number_of_sequences*2,2):
        key = splits[i]
        key_filtered = re.findall(r"[1-9]\d?",key)
        mapping_dict[key_filtered[0]]= splits[i+1]
    return mapping_dict




class given_tree_quartet_sampler_from_fasta(BaseTreeSimulator):
    def __init__(self, given_phylotree,sequence_mapping,random = False):
        """Create a quartet tree subsampler
        """
        super().__init__()
        self.big_tree = given_phylotree
        self.sequence_mapping = sequence_mapping
        self.random = random
        self.state = 0
        self.iterator = itertools.cycle(itertools.combinations(self.big_tree.get_leaves(),4))
        self.num_of_leafs = len(self.big_tree.get_leaves())
        self.num_of_quartets = int(math.factorial(self.num_of_leafs)/(math.factorial(4)*math.factorial(self.num_of_leafs-4)))

        #self.lba_ratio = lba_ratio
        
    def generate(self):
        tree = self.big_tree.copy()
        
        if self.random == True:
            leaves = random.sample(tree.get_leaves(),4)
        else:
            leaves = list(next(self.iterator))
            leaves = random.sample(leaves,4)
            leafnames = [l.name for l in leaves]
            leaves = [tree.search_nodes(name=ln)[0] for ln in leafnames ]
            self.state += 1

        tree.prune(leaves)
        tree = self.assign_sequences_to_leafs(tree)
        tree.unroot()
        return tree

    
    def assign_sequences_to_leafs(self,tree):
        for leaf in tree: 
            leaf.add_features(sequence = self.sequence_mapping[leaf.name])
        return tree
    
    
    
#this one assumes that there are already sequences at the leaves
class given_tree_quartet_sampler_from_full_tree(BaseTreeSimulator):
    def __init__(self, given_phylotree,random = False):
        """Create a quartet tree subsampler
        """
        super().__init__()
        self.big_tree = given_phylotree
        self.random = random
        self.state = 0
        self.iterator = itertools.cycle(itertools.combinations(self.big_tree.get_leaves(),4))
        self.num_of_leafs = len(self.big_tree.get_leaves())
        self.num_of_quartets = int(math.factorial(self.num_of_leafs)/(math.factorial(4)*math.factorial(self.num_of_leafs-4)))
        
    def generate(self):
        tree = self.big_tree.copy()
        
        if self.random == True:
            leaves = random.sample(tree.get_leaves(),4)
        else:
            leaves = list(next(self.iterator))
            leaves = random.sample(leaves,4)
            leafnames = [l.name for l in leaves]
            leaves = [tree.search_nodes(name=ln)[0] for ln in leafnames ]
            self.state += 1

        tree.prune(leaves)
        tree.unroot()
        return tree

    
    def assign_sequences_to_leafs(self,tree):
        for leaf in tree: 
            leaf.add_features(sequence = self.sequence_mapping[leaf.name])
        return tree
    
    
    
    
        

def create_pickle_from_quartet_sampler(quartet_sampler,path,num = None,random = False):
    #pathlib.Path(path)
    if num == None:
        num = quartet_sampler.num_of_quartets
        
    for i in range(num):
        with open(path+"/{:06d}.pickle".format(i),"wb") as file:
            t = quartet_sampler.generate()
            pickle.dump(t, file, protocol=pickle.HIGHEST_PROTOCOL)