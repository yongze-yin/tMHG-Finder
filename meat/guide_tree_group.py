import string
from typing import List
import dendropy
import math
import copy
from Bio import Phylo
from collections import defaultdict

def tree_height(node):
    """
    Given a dendropy root node, calculate the level height of the tree.
    """
    return 0 if node.is_leaf() else max([1+tree_height(child) for child in node.child_nodes()])

def shortest_reroot(guide_tree_path: string) -> string:
    """
    Given a directory to a guide tree, reroot on every node and pick the one with
    the least level height, return the newick tree string. 
    """
    guide_tree = dendropy.Tree.get(path = guide_tree_path, schema = 'newick', rooting='force-rooted')
    best_tree, least_height = None, math.inf
    for node in guide_tree.postorder_node_iter():
        guide_tree_node_copy = copy.deepcopy(guide_tree)
        guide_tree_node_copy.reroot_at_node(copy.deepcopy(node))
        # print("taxa amount", len([str(taxon).replace("'","") for taxon in guide_tree_node_copy.taxon_namespace]))
        height = tree_height(guide_tree_node_copy.seed_node)
        if height < least_height:
            best_tree = guide_tree_node_copy
            least_height = height
    # print('best tree taxa amount',len([str(taxon).replace("'","") for taxon in best_tree.taxon_namespace]))
    # return best_tree.as_string(schema='newick')
    return best_tree

def internal_node_leaves(internal):
    """
    Inputting a dendropy tree node, return a concatenation ("|") of leaf taxa under the left child, and a concatenation("|")
    of leaf taxa under the right child. The leaf child leaves and right child leaves are concatenated after by a ","
    """
    children = internal.child_nodes()
    right_child = children[0]
    left_child = children[1]
    right_child_leaves = [str(n.taxon).replace("'","") for n in right_child.leaf_nodes()]
    left_child_leaves = [str(n.taxon).replace("'","") for n in left_child.leaf_nodes()]
    return ','.join(["|".join(right_child_leaves),"|".join(left_child_leaves)])
    # print(len(internal.leaf_nodes()),internal, [str(n.taxon).replace("'","") for n in internal.leaf_nodes()])
    # print(len(left_child.leaf_nodes()),left_child, [str(n.taxon).replace("'","") for n in left_child.leaf_nodes()])
    # print(len(right_child.leaf_nodes()),right_child, [str(n.taxon).replace("'","") for n in right_child.leaf_nodes()])

def initial_taxa_internal(tree):
    """
    Traverse every internal node and collect the leaves belong to left child and right child.
    Output an inital list of leaf pairs to be partitioned at each internal node. And also output
    the set of taxa included in the initial tree.
    """
    internal_node_list = tree.internal_nodes(exclude_seed_node=False)
    time_node = defaultdict(lambda: [])
    taxa = [str(taxon).replace("'","") for taxon in tree.taxon_namespace]
    for internal in internal_node_list:
        leaf_string = internal_node_leaves(internal)
        level = max([branch.count('|') for branch in leaf_string.split(',')])
        below_taxa_count = 2 + leaf_string.count("|")
        time_node[level].append(leaf_string)
    remaining_pair = []
    for level in sorted(time_node.keys()):
        remaining_pair += time_node[level]
    visited_node_MHG = {taxon:None for taxon in taxa}
    return visited_node_MHG, remaining_pair
        
def give_me_the_next_visit(visited_node_MHG, remaining_pair):
    """
    show the next available/to be partitioned internal node
    take the two children node out of the vistied_node_MHG;
    Dont forget to update them after finish partitioning
    Output: boolean, dict, updated_visited_node_MHG, updated_remaining_pair

    boolean: whether the next pair is ready to MHG or not
    dict: to be parititoned two children nodes with their current MHG
    updated_*: if boolean if true, remove the nodes from the record; else return the record.
    """
    for node in remaining_pair:
        two_nodes = set(node.split(','))
        if two_nodes.issubset(set(visited_node_MHG.keys())):
            ready_MHG = {}
            remaining_pair.remove(node)
            for n in two_nodes:
                ready_MHG[n] = copy.deepcopy(visited_node_MHG[n])
            return True, node, ready_MHG, visited_node_MHG, remaining_pair
    return False, None, None, visited_node_MHG, remaining_pair