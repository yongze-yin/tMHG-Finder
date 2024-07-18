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

def shortest_reroot(guide_tree_path: string, reroot = True):
    """
    Given a directory to a guide tree, reroot on every node and pick the one with
    the least level height, return the newick tree obj. 
    """
    guide_tree = dendropy.Tree.get(path = guide_tree_path, schema = 'newick', rooting='force-rooted')
    if not reroot:
        return guide_tree
    tip_cnt = len(guide_tree.leaf_nodes())
    best_tree, least_height = guide_tree, math.inf
    for node in guide_tree.postorder_node_iter():
        guide_tree_node_copy = copy.deepcopy(guide_tree)
        guide_tree_node_copy.reroot_at_node(copy.deepcopy(node))
        if len(guide_tree_node_copy.leaf_nodes())!= tip_cnt:
            continue
        height = tree_height(guide_tree_node_copy.seed_node)
        if height < least_height:
            best_tree = guide_tree_node_copy
            least_height = height
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
    merged_children_string = ','.join(["|".join(right_child_leaves),"|".join(left_child_leaves)])
    if len(children) == 2:
        # Tree is bifurcating
        return [merged_children_string]
    else:
        # Tree is trifurcating
        other_child = children[2]
        other_child_leaves = [str(n.taxon).replace("'","") for n in other_child.leaf_nodes()]
        brother_merged_leaves = right_child_leaves + left_child_leaves
        super_merged_children_string = ','.join(["|".join(other_child_leaves),"|".join(brother_merged_leaves)])
        return [merged_children_string, super_merged_children_string]
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
        leaf_string_list = internal_node_leaves(internal)
        for leaf_string in leaf_string_list:
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
        two_nodes = node.split(',')
        left_n, right_n = two_nodes[0], two_nodes[1]
        left_n_tup, right_n_tup = tuple(sorted(left_n.split("|"))), tuple(sorted(right_n.split("|")))
        sorted_tuple_key_dict = {tuple(sorted(key.split("|"))):key for key in visited_node_MHG}
        if left_n_tup in sorted_tuple_key_dict and right_n_tup in sorted_tuple_key_dict:
            ready_MHG = {}
            remaining_pair.remove(node)
            ready_MHG[sorted_tuple_key_dict[left_n_tup]] = copy.deepcopy(visited_node_MHG[sorted_tuple_key_dict[left_n_tup]])
            ready_MHG[sorted_tuple_key_dict[right_n_tup]] = copy.deepcopy(visited_node_MHG[sorted_tuple_key_dict[right_n_tup]])
            return True, node, ready_MHG, visited_node_MHG, remaining_pair
    return False, None, None, visited_node_MHG, remaining_pair