#!/usr/bin/env python3

# The MIT License (MIT)
# Copyright (c) 2019 Thomas Huetter.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''
Program: Analyze DNA sequences based on the manhattan distance and the ranks of
so-called nucleotide vectors of a query- and a reference group.
'''

import argparse
import csv
import os
from math import log
from Bio import Phylo
from Bio import AlignIO

#### GLOBAL DICTS
# Maps a dna sequence character to a vector index.
char_to_vector_idx = {
    'a': 0,
    'c': 1,
    'g': 2,
    't': 3,
    '-': 4
}

# Maps a combination of two dna sequence characters to a vector index.
two_char_to_vector_idx = {
    'aa': 0, 'ac': 1, 'ag': 2, 'at': 3, 'a-': 4,
    'ca': 5, 'cc': 6, 'cg': 7, 'ct': 8, 'c-': 9,
    'ga': 10, 'gc': 11, 'gg': 12, 'gt': 13, 'g-': 14,
    'ta': 15, 'tc': 16, 'tg': 17, 'tt': 18, 't-': 19,
    '-a': 20, '-c': 21, '-g': 22, '-t': 23, '--': 24
}


def nuc_vectors_one_positions(list_of_alignments):
    """
        Given a list of alignments, this function returns a list of lists,
        containing the nucleotide vectors for each position. We count the
        number of occurences per base (a,c,g,t,-) and store them in a list
        [|a|,|t|,|c|,|g|,|-|].
    """

    # Stop if there are no alignments in the list.
    if not list_of_alignments:
        raise SystemError("'nuc_vectors_one_positions': no alignment in \
            list_of_alignments")

    # Holds a list with nucleotide vectors per alignment position.
    nuc_vectors = []

    # Assuming that all alignments are of the same size.
    for pos in range(len(list_of_alignments[0])):
        nuc_vector = [0, 0, 0, 0, 0]
        # Sum up the number of different bases over all alignments at position
        # pos.
        for alignment in list_of_alignments:
            # Nucleotide at current position in lower case.
            nucleotide = alignment.seq[pos].lower()
            # If statement needed to fix input characters other a, c, g, t, -.
            # Treat all others as '-'.
            if nucleotide in char_to_vector_idx:
                nuc_vector[char_to_vector_idx[nucleotide]] +=1
            else:
                nuc_vector[char_to_vector_idx['-']] +=1
        nuc_vectors.append(nuc_vector)

    return nuc_vectors


def nuc_vectors_two_positions(list_of_alignments, combinations):
    """
        Returns the base combination of all positions in the
        alignment.
    """

    # Stop if there are no alignmentsn or combinations in the list.
    if not list_of_alignments:
        raise SystemError("'nuc_vectors_two_positions_k': no alignment in \
            list_of_alignments")
    if not combinations:
        raise SystemError("'nuc_vectors_two_positions_k': no combinations \
            in combinations")

    # Holds a list with nucleotide vectors per alignment position combination.
    nuc_vectors = []

    # Assuming that all alignments are of the same size.
    for i, j in combinations:
        # Given 5 possible characters, there are 25 combinations of two
        # positions.
        nuc_vector = 25 * [0]
        # Sum up the number of different bases over all alignments at
        # position pos.
        for alignment in list_of_alignments:
            # Combined nucleotide (i, j); nucleotide at position i
            # in lower case.
            nucleotide_i = alignment.seq[i].lower()
            # Fix for different input characters than a, c, g, t, -. Treat
            # all others as '-'.
            base_string = ""
            if nucleotide_i not in char_to_vector_idx: base_string += "-"
            else: base_string += nucleotide_i

            # Combined nucleotide (i, j); nucleotide at position j
            # in lower case.
            nucleotide_j = alignment.seq[j].lower()
            if nucleotide_j not in char_to_vector_idx: base_string += "-"
            else: base_string += nucleotide_j

            nuc_vector[two_char_to_vector_idx[base_string]] += 1
        nuc_vectors.append([(i, j), nuc_vector])

    return nuc_vectors


def assign_pre_id_to_inner_nodes(phylo_tree):
    """
        Replace the name of the inner nodes of a given phylogenetic tree with
        its preorder number in the tree.
    """
    idx = 0
    for node in phylo_tree.find_clades(terminal=False, order='preorder'):
        node.name = '%d' % (idx)
        idx += 1

    return phylo_tree



def get_all_terminals_of_inner_node(phylo_tree, node_id):
    """
        Given a phylogenetic tree and a preorder ID of an inner node, return
        a list of all leaf nodes of the subtree rooted at the inner node.
    """

    # Since the node IDs are unique, the loop is only iterated once.
    leaf_nodes = list()
    for clade in phylo_tree.find_clades(name=node_id, order='preorder'):
        leaf_nodes += clade.get_terminals(order='postorder')

    return leaf_nodes


def get_alignments_of_group(group_node_ids, alignments):
    """
        Given a list of node IDs and a list of alignments, return a list with
        all alignments that belong the nodes in the group.
    """

    # Create an inverted list of alignments names from all alignments.
    inverted_list_alignments = {}
    alignment_idx = 0
    for alignment in alignments:
        inverted_list_alignments[alignment.name] = alignment_idx
        alignment_idx += 1

    # Lookup the inverted list to get the index of the alignment based on the
    # node names specified in the group node ids.
    group_alignments = list()
    for node in group_node_ids:
        if node.name in inverted_list_alignments:
            group_alignments.append(
                alignments[inverted_list_alignments[node.name]])

    return group_alignments


def list_difference(list1, list2):
    """
        Given two lists with alignments list1 and list2, return a new list
        (new_list2) that contains all elements of list2 without elements of
        list1.
    """

    # Create an inverted list of alignments names in list1 to allow fast
    # lookups.
    inverted_list1 = {}
    for alignment in list1:
        inverted_list1[alignment.name] = 1

    # Copy only elements of list2 that are not in list1, by looking up the
    # inverted list.
    new_list2 = list()
    for alignment in list2:
        if alignment.name not in inverted_list1:
            new_list2.append(alignment)

    return new_list2


def norm_manhattan_distance(nuc_vector1, nuc_vector2):
    """
        Given two lists containing nucleotide vectors, this functions returns
        the normalized Manhattan distance between these two vectors.
    """

    return (sum(abs(a-b) for a,b in zip(nuc_vector1, nuc_vector2))) / (sum(nuc_vector1) + sum(nuc_vector2))


def get_rank(nuc_vector, one_pos=True):
    """
        Given a nucleotide vector, return the rank and the base that occurs the
        most.
    """

    rank = 0
    max_base_nr = 0
    max_base_pos = 0

    for element in range(len(nuc_vector)):
        if nuc_vector[element] != 0:
            rank += 1
            if nuc_vector[element] > max_base_nr:
                max_base = position_to_base(element, one_pos)
                max_base_pos = element
                max_base_nr = nuc_vector[element]

    return rank, position_to_base(max_base_pos, one_pos)


def position_to_base(pos, one_pos=True):
    """
        Given a base vector index, return the according base.
    """
    if one_pos:
        if pos == 0: return 'a'
        elif pos == 1: return 'c'
        elif pos == 2: return 'g'
        elif pos == 3: return 't'
        elif pos == 4: return '-'
        else: raise SystemError("'position_to_base': one_pos index out of bounds")
    else:
        if pos == 0: return 'aa'
        elif pos == 1: return 'ac'
        elif pos == 2: return 'ag'
        elif pos == 3: return 'at'
        elif pos == 4: return 'a-'
        elif pos == 5: return 'ca'
        elif pos == 6: return 'cc'
        elif pos == 7: return 'cg'
        elif pos == 8: return 'ct'
        elif pos == 9: return 'c-'
        elif pos == 10: return 'ga'
        elif pos == 11: return 'gc'
        elif pos == 12: return 'gg'
        elif pos == 13: return 'gt'
        elif pos == 14: return 'g-'
        elif pos == 15: return 'ta'
        elif pos == 16: return 'tc'
        elif pos == 17: return 'tg'
        elif pos == 18: return 'tt'
        elif pos == 19: return 't-'
        elif pos == 20: return '-a'
        elif pos == 21: return '-c'
        elif pos == 22: return '-g'
        elif pos == 23: return '-t'
        elif pos == 24: return '--'
        else: raise SystemError("'position_to_base': two_pos index out of bounds")


def two_position_to_base(pos):
    """
        Given a base vector index, return the according base.
    """

    if pos == 0: return 'a'
    elif pos == 1: return 'c'
    elif pos == 2: return 'g'
    elif pos == 3: return 't'
    elif pos == 4: return '-'
    else: raise SystemError("'position_to_base': index out of bounds")


def assign_indices(sorted_one_pos_ranking):
    """
        Given the output of the one positional analysis, this function returns
        the index of each position. As defined by Waegele and Mayer (2007), we
        only consider positions that are unique within the query group and
        distinguish three indices: binary, asymmetric, and noisy.
        We identify them based on the discriminative power, the rank of the
        query group, and the rank of the reference group.
    """

    binary = list()
    asymmetric = list()
    noisy = list()

    # Assign an index for all positions. pos_info is structured as follows
    # [pos, disc. power, (query rank, query max base), (ref rank, ref max base)]
    # An example for pos_info is [434, 1.0, "(1, 'g')", "(1, 'a')"].
    for pos_info in sorted_one_pos_ranking:
        # We only consider positions with unique bases in the query group (query
        # rank = 1).
        if pos_info[2][0] == 1 and pos_info[2][1] != '-':
            # A discriminative power of 1 is either binary or asymmetric.
            if pos_info[1] == 1:
                # A position is binary if the rank of the reference group is 1.
                # Otherwise, the position is asymmetric.
                if pos_info[3][0] == 1:
                    pos_info.append("binary")
                    binary.append(pos_info)
                else:
                    pos_info.append("asymmetric")
                    asymmetric.append(pos_info)
            else:
                # In case that the discriminative power is less than 1, the
                # position is noisy. In case that the reference rank is 1, the
                # reference and query group have the same unique base. Hence, it
                # is not noisy.
                if pos_info[3][0] > 1:
                    pos_info.append("noisy")
                    noisy.append(pos_info)

    return binary, asymmetric, noisy


def one_positional_calculation(query_group, reference_group):
    """
        Analyse each position of by the Manhattan distance of their base
        vectors and the number of different bases.
    """

    one_pos_index = list()
    nuc_vectors_query_group = nuc_vectors_one_positions(query_group)
    nuc_vectors_ref_group = nuc_vectors_one_positions(reference_group)

    for pos in range(len(nuc_vectors_query_group)):
        one_pos_index.append([pos + 1,
            norm_manhattan_distance(nuc_vectors_query_group[pos],
                nuc_vectors_ref_group[pos]),
            get_rank(nuc_vectors_query_group[pos]),
            get_rank(nuc_vectors_ref_group[pos])
        ])

    return one_pos_index

def number_of_different_bases(nuc_vector_query_group, nuc_vector_ref_group):
    """
        Returns the number of different bases at a specific position in the
        reference as well as the control group.
    """

    counter = 0
    for pos in range(len(nuc_vector_query_group)):
        if (nuc_vector_query_group[pos] + nuc_vector_ref_group[pos]) > 0:
            counter += 1

    return counter


def two_positional_calculation(query_group, ref_group, combinations):
    """
        Analyse the combination of two positions by the Manhattan distance of
        their nucleotide vectors and the rank of th.
    """

    binary = list()
    asymmetric = list()
    noisy = list()

    nuc_vectors_query_group = nuc_vectors_two_positions(query_group,
        combinations)
    nuc_vectors_ref_group = nuc_vectors_two_positions(ref_group, combinations)

    # There are as many entries in the two nucleotide vectors as combinations.
    for i in range(len(combinations)):
        # nuc_vectors_query_group[i] = [(4,67), [0,4,3,12,...]]
        disc_power = norm_manhattan_distance(nuc_vectors_query_group[i][1],
            nuc_vectors_ref_group[i][1])
        query_rank, query_max_base = get_rank(nuc_vectors_query_group[i][1],
            False)
        ref_rank, ref_max_base = get_rank(nuc_vectors_ref_group[i][1], False)

        pos_info = [(combinations[i][0]+1, combinations[i][1]+1), disc_power,
            (query_rank, query_max_base), (ref_rank, ref_max_base)]

        # We only consider positions with unique bases in the query group.
        if query_rank == 1 and query_max_base != '-':
            if disc_power == 1:
                # Binary if the reference group has only one but another base.
                # Otherwise, asymmetric.
                if ref_rank == 1:
                    pos_info.append("binary")
                    binary.append(pos_info)
                else:
                    pos_info.append("asymmetric")
                    asymmetric.append(pos_info)
            else:
                # In case that the discriminative power is less than 1, the
                # position is noisy.
                if ref_rank > 1:
                    pos_info.append("noisy")
                    noisy.append(pos_info)

                # Otherwise, there is no index for this position.

    return binary, asymmetric, noisy


def shannon_entropy_calculation(alignments):
    """
        Computes the Shannon entropy for each position in the given list of
        alignments.
    """

    nuc_vectors = nuc_vectors_one_positions(alignments)
    number_of_alignments = len(alignments)
    entropy_per_pos = list()

    for pos in range(len(nuc_vectors)):
        # Entropy H = E[I] = sum[z in Z] -pz * log2(pz)
        # Z = [a, c, g, t, -]
        entropy = 0.0
        for base_idx in range(5):
            if nuc_vectors[pos][base_idx] > 0:
                p = nuc_vectors[pos][base_idx] / number_of_alignments
                entropy += -p * log(p, 2)
        entropy_per_pos.append([pos+1, entropy])

    return entropy_per_pos


def averages_calculation(shannon_entropy, k):
    """
        Computes the averages for a group of k shannon entropy values.
    """

    entropies = [i[1] for i in shannon_entropy]
    averages = list()

    # For all entropy values.
    for val in range(len(entropies)):
        # Compute the lower- and upperbound.
        lb = val - k
        ub = val + k

        # Check bounds.
        if lb < 0:
            lb = 0
        if ub > len(entropies):
            ub = len(entropies)

        # Calculate the average of the given entropy values.
        averages.append(sum(entropies[lb:ub]) / len(entropies[lb:ub]))

    return averages


def index_analysis(alignments, phylo_tree, query_group_id, ref_group_id,
                    filename=None, order_pos=False):
    """
        Wrapper for index calculation. Assembles the alignment groups and writes
        the index result in a csv file.
    """

    # Perform one positional analysis and index positions based on
    # discriminative power and ranks of the query and reference group.
    sorted_one_pos_ranking = one_positional_analysis(alignments, phylo_tree,
        query_group_id, ref_group_id, None, order_pos)

    # Assign an index for each position.
    binary, asymmetric, noisy = assign_indices(sorted_one_pos_ranking)

    # Write result to csv file.
    if filename is not None:
        write_csv_file(filename, ["position", "disriminative power",
            "query rank", "reference_rank", "index"],
            binary + asymmetric + noisy)

    return binary, asymmetric, noisy


def one_positional_analysis(alignments, phylo_tree, query_group_id,
                            ref_group_id, filename=None, order_pos=False):
    """
        Wrapper for one positional analysis. Assembles the alignment groups
        and writes the one positional analysis results in a csv file.
    """

    # Assemble the alignments according to the specified groups.
    alignments_query_group = get_alignments_of_group(
        get_all_terminals_of_inner_node(phylo_tree, query_group_id), alignments)
    alignments_ref_group = get_alignments_of_group(
        get_all_terminals_of_inner_node(phylo_tree, ref_group_id), alignments)
    alignments_ref_group = list_difference(alignments_query_group,
                                                alignments_ref_group)

    # Perform the one positional analysis per position.
    pos_ranking = one_positional_calculation(
        alignments_query_group, alignments_ref_group)

    # Sort descending by Manhatten distance and ascending by query and reference
    # rank or by position.
    if order_pos:
        sorted_pos_ranking = sorted(pos_ranking, key = lambda y: (y[0]))
    else:
        sorted_pos_ranking = sorted(pos_ranking,
            key = lambda y: (-y[1], y[2][0], y[3][0]))

    # Write result to csv file.
    if filename is not None:
        write_csv_file(filename, ["position", "disriminative power",
            "query rank", "reference_rank"], sorted_pos_ranking)

    return sorted_pos_ranking


def signature_nucleotide_detection(alignments, phylo_tree, query_group_id,
                                    ref_group_id, threshold=1, filename=None,
                                    order_pos=False, filter_sig_nuc=False):
    """
        Wrapper for the signature nucleotide detection algorithm. Given a query
        and reference group, all positions are ranked by manhattan distance and
        query and reference rank. Further they are indexed, filtered, and
        ordered based on the in the input parameters. Given a threshold for the
        k-window, all noisy positions in this k-window are considered in
        combinations of positions in order to make them asymmetric. These
        asymmetric combinations are included in the final ranking.
    """

    # Assemble the alignments according to the specified groups.
    alignments_query_group = get_alignments_of_group(
        get_all_terminals_of_inner_node(phylo_tree, query_group_id), alignments)
    alignments_ref_group = get_alignments_of_group(
        get_all_terminals_of_inner_node(phylo_tree, ref_group_id), alignments)
    alignments_ref_group = list_difference(alignments_query_group,
                                                alignments_ref_group)

    # Perform one positional analysis and index positions based on
    # discriminative power and ranks of the query and reference group.
    one_pos_ranking = one_positional_analysis(alignments, phylo_tree,
        query_group_id, ref_group_id, None, order_pos)

    # Assign an index for each position.
    binary, asymmetric, noisy = assign_indices(one_pos_ranking)

    # If a threshold is given, perform two position analysis.
    # Do not if there are too little noisy positions.
    if threshold > 1 and len(noisy) >= 2:

        # Sort noisy list by position.
        sorted_noisy = sorted(noisy, key = lambda y: (y[0]))

        # All possible combinations of single noisy positions.
        combinations = list()
        for i in range(len(sorted_noisy)):
            for j in range(i+1, len(sorted_noisy)):
                if sorted_noisy[j][0] - sorted_noisy[i][0] >= threshold:
                    break
                combinations.append((sorted_noisy[i][0]-1, sorted_noisy[j][0]-1))

        # Compute the index for two position combinations. Consider only single
        # noisy positions and look for combinations that become asymmetric.
        b, a, n = two_positional_calculation(alignments_query_group,
                alignments_ref_group, combinations)

        # Extend the single asymmetric positions with the combined ones.
        one_pos_ranking.extend(a)

    # Extended one pos ranking, including two positional results.
    ranking = one_pos_ranking

    # Filter the result such that only indexed positions are returned.
    if filter_sig_nuc:
        ranking = [x for x in one_pos_ranking if len(x) is 5]

    # Order either by position or ranking.
    if order_pos:
        ranking_sorted = sorted(ranking,
            key=lambda y: y if isinstance(y, int) else y[0])
    else:
        ranking_sorted = sorted(ranking,
            key = lambda y: (-y[1], y[2][0], y[3][0]))

    # Write result to csv file.
    if filename is not None:
        write_csv_file(filename, ["position", "disriminative power",
            "query rank", "reference_rank", "index"], ranking_sorted)

    return ranking_sorted


def shannon_entropy_analysis(alignments, phylo_tree, group_id = None,
                            filename=None):
    """
        Wrapper for shannon entropy. Computes and writes the shannon entropy
        for each position in the reference group.
    """

    # If a group id is given, assemble the alignment based on the group.
    # Otherwise, take the full alignment.
    if group_id is not None:
        alignment = get_alignments_of_group(get_all_terminals_of_inner_node(
            phylo_tree, group_id), alignments)
    else:
        alignment = alignments

    # Perform the shannon entropy analysis per position.
    shannon_entropy = shannon_entropy_calculation(alignment)

    # Calculate the mean for each position based on the next and previous k
    # values.
    averages = averages_calculation(shannon_entropy, 10)

    for val in range(len(shannon_entropy)):
        shannon_entropy[val].append(averages[val])

    # Write result to csv file.
    if filename is not None:
        write_csv_file(filename, ["position", "shannon_entropy", "average"],
            shannon_entropy)

    return shannon_entropy


def write_csv_file(filename, header, data):
    """
        Writes the header and the data items into a file with a given filename.
    """

    writer = csv.writer(open(filename, "w"))
    writer.writerow(header)
    writer.writerows(data)


def get_depth(phylo_tree):
    """
        Returns the depth of a tree.
    """
    depth = 0
    for terminal_node in phylo_tree.get_terminals(order='preorder'):
        path_length = len(phylo_tree.get_path(target=terminal_node))
        if path_length > depth:
            depth = path_length

    return depth


def main():

    # Parse input arguments.  (required=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", type=str,
        help="Path to input file that contains the alignments.")
    parser.add_argument("--alignment_format", type=str, default="fasta",
        help="Specify the alignment file format. Default is set to 'fasta'.")
    parser.add_argument("--phylotree", type=str,
        help="Path to input file that contains the phylogenetic tree.")
    parser.add_argument("--phylotree_format", type=str, default="newick",
        help="Specify the phylogenetic tree file format. Default is set to \
        'newick'.")
    parser.add_argument("--query_group_id", type=str, help="Preorder ID of \
        an inner node of the phylogenetic tree. All species (leaf nodes) in \
        that subtree are considered in the query group of the analysis.")
    parser.add_argument("--reference_group_id", type=str, help="Preorder ID of \
        an inner node of the phylogenetic tree. All species (leaf nodes) in \
        that subtree are considered in the reference group of the analysis.")
    parser.add_argument("--print_tree", default=False, action="store_true",
        help="Print the tree using the Biopython draw function.")
    parser.add_argument("--signature_nucleotides",default=False,
        action="store_true", help="Perform the signature nucleotide analysis.")
    parser.add_argument("--k_window", type=int, default=1,
        help="Combine two noisy positions in order to find more asymmetric \
        pairs resulting in a higher number of unique characteristics. The \
        required parameter k is used for the so-called k-window. Only noisy \
        positions in a range of k are considered to be combined for analysis.")
    parser.add_argument("--filter", default=False, action="store_true",
        help="Only save and/or display signature nucleotides.")
    parser.add_argument("--save", default=False, action="store_true",
        help="Write the results of the analysis to a file.")
    parser.add_argument("--pos_order", default=False, action="store_true",
        help="Set this flag to order the result by position instead of the \
        ranking.")
    parser.add_argument("--shannon_entropy",default=False, action="store_true",
        help="Compute the shannon entropy for each position in the reference \
        group.")
    args = parser.parse_args()

    # Used for csv output files.
    filename = os.path.splitext(os.path.basename(args.alignment))[0]

    # Parse input files and assign preorder ID to inner nodes of the
    # phylogenetic tree.
    alignments = AlignIO.read(args.alignment, args.alignment_format)
    phylo_tree = Phylo.read(args.phylotree, args.phylotree_format)
    assign_pre_id_to_inner_nodes(phylo_tree)

    # Perform signature nucleotide detection.
    if args.signature_nucleotides:
        sig_nuc_filename = None
        if args.save:
            sig_nuc_filename = filename + "_signature_nucleotides.csv"

        signature_nucleotide_detection(alignments, phylo_tree,
            args.query_group_id, args.reference_group_id, args.k_window,
            sig_nuc_filename, args.pos_order, args.filter)

    # Perform shannon entropy calculation.
    if args.shannon_entropy:
        entropy_filename = None
        if args.save:
            entropy_filename = filename + "_entropy.csv"

        shannon_entropy_analysis(alignments, phylo_tree,
            args.query_group_id, entropy_filename)

    # Print the input phylogenetic tree.
    if args.print_tree:
        Phylo.draw(phylo_tree)



if __name__ == "__main__":
    main()
