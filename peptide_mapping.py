# Copyright (c) 2025 Christian M. Zmasek
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to
# whom the Software is furnished to do so, subject to the
# following conditions:
#
# The above copyright notice and this permission notice shall
# be included in all copies or substantial portions of the
# Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
import argparse as ap
import os
import re
import sys

import fasta_parser
import target_match


class PeptideMapping(object):
    VERSION = '1.0.0'

    ID_RE = re.compile(">\\s*(.+)")
    GAP_RE = re.compile("[-\\s]+")
    WS_RE = re.compile("\\s+")
    COMBINE_WHITESPACE_RE = re.compile(r"\s+")

    @staticmethod
    def perform_mapping_windowed(query, target, query_start, query_end, window_start, window_end):
        if window_end > len(target) - 1 or window_start < 0:
            raise ValueError("Start or end are out of bounds.")
        if window_end <= window_start:
            raise ValueError("Start must be smaller than end.")

        return PeptideMapping.perform_mapping(query, target[window_start: window_end + 1], query_start, query_end, 12,
                                              0.1)

    @staticmethod
    def hamming_distance(s1, s2):
        """Calculate the Hamming distance between two strings."""
        if len(s1) != len(s2):
            raise ValueError("The strings must have the same length.")
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    @staticmethod
    def inverted_normalized_hamming_distance(s1, s2):
        return PeptideMapping.calc_inverted_normalized_hamming_distance(PeptideMapping.hamming_distance(s1, s2),
                                                                        len(s1))

    @staticmethod
    def calc_inverted_normalized_hamming_distance(dist, length):
        return 1 - dist / length

    @staticmethod
    def calculate_score(dist, query_center_coord, target_center_coord, coord_distance_weight):
        return dist + coord_distance_weight * abs(query_center_coord - target_center_coord)

    @staticmethod
    def perform_mapping(query, target, query_start, query_end, max_dist, coord_distance_weight):
        if len(target) <= len(query):
            raise ValueError("Query must be shorter than target.")
        # least_dist = len(query) + 1
        best_matches = []
        query_len = len(query)
        for target_start in range(len(target) - query_len + 1):
            t = target[target_start:target_start + query_len]
            dist = PeptideMapping.hamming_distance(t, query)
            # if dist <= least_dist:
            if dist <= max_dist:
                # if dist < least_dist:
                #    best_matches = []
                #    least_dist = dist
                target_end = target_start + query_len - 1

                query_center_coord = query_start + (query_end - query_start) / 2
                target_center_coord = target_start + (target_end - target_start) / 2
                score = PeptideMapping.calculate_score(dist, query_center_coord, target_center_coord,
                                                       coord_distance_weight)

                best_matches.append(target_match.TargetMatch(query_start, query_end, query, target_start, target_end, t,
                                                             dist, score))

        return best_matches

    @staticmethod
    def transpose(in_table, seq_file, protein_name, outfile, max_dist, coord_distance_weight, verbose):

        if os.path.isfile(outfile):
            print(outfile + ' already exists')
            sys.exit()

        seqs = fasta_parser.parse_fasta_file(seq_file, remove_gaps=True)
        seq = seqs[0]
        target_seq = seq.get_seq()
        of = open(outfile, 'w')

        PeptideMapping.write_header(of)

        with open(in_table, 'r') as file:
            for line in file:
                line = line.strip()
                if not line.startswith("#"):
                    values = line.split("\t")
                    query_seq = values[0]
                    orf = values[1]
                    protein = values[2]
                    query_from = int(values[3]) - 1
                    query_to = int(values[4]) - 1

                    if protein == protein_name:
                        m = PeptideMapping.perform_mapping(query_seq, target_seq, query_from, query_to, max_dist,
                                                           coord_distance_weight)
                        if len(m) > 0:
                            m.sort()
                            m0 = m[0]
                            norm = PeptideMapping.calc_inverted_normalized_hamming_distance(m0.get_distance(),
                                                                                            m0.calc_length())

                            coord_dist = abs(
                                (query_from + (query_to - query_from) / 2) - (
                                        m0.get_target_start() + (m0.get_target_end() - m0.get_target_start()) / 2))

                            of.write(query_seq)
                            of.write("\t")
                            of.write(protein_name)
                            of.write("\t")
                            of.write(orf)
                            of.write("\t")
                            of.write(str(query_from))
                            of.write("\t")
                            of.write(str(query_to))
                            of.write("\t")
                            of.write(m0.get_target_sequence())
                            of.write("\t")
                            of.write(str(m0.get_target_start()))
                            of.write("\t")
                            of.write(str(m0.get_target_end()))
                            of.write("\t")
                            of.write(str(norm))
                            of.write("\t")
                            of.write(str(m0.get_distance()))
                            of.write("\t")
                            of.write(str(m0.get_score()))
                            of.write("\t")
                            of.write(str(coord_dist))
                            of.write("\t")
                            of.write(str(len(m)))
                            of.write("\n")
                        else:
                            of.write(query_seq)
                            of.write("\t")
                            of.write(protein_name)
                            of.write("\t")
                            of.write(orf)
                            of.write("\t")
                            of.write(str(query_from))
                            of.write("\t")
                            of.write(str(query_to))
                            of.write("\t")
                            of.write("")
                            of.write("\t")
                            of.write("")
                            of.write("\t")
                            of.write("")
                            of.write("\t")
                            of.write("")
                            of.write("\t")
                            of.write("")
                            of.write("\t")
                            of.write("")
                            of.write("\t")
                            of.write("")
                            of.write("\t")
                            of.write(str(len(m)))
                            of.write("\n")

                        if verbose:
                            for i, x in enumerate(m):
                                norm = PeptideMapping.calc_inverted_normalized_hamming_distance(x.get_distance(),
                                                                                                x.calc_length())
                                print(str(i) + ": ")
                                print("Hamming distance     : " + str(x.get_distance()))
                                print("Norm Hamming distance: " + str(norm))
                                print("Score                : " + str(x.get_score()))
                                print(x.to_aln())
                                print("----------------------------------------")
        of.close()

    @staticmethod
    def write_header(of):
        of.write("#QUERY SEQ")
        of.write("\t")
        of.write("PROTEIN")
        of.write("\t")
        of.write("ORF")
        of.write("\t")
        of.write("QUERY FROM")
        of.write("\t")
        of.write("QUERY FROM")
        of.write("\t")
        of.write("TARGET SEQ")
        of.write("\t")
        of.write("TARGET FROM")
        of.write("\t")
        of.write("TARGET TO")
        of.write("\t")
        of.write("NORM HAMMING DIST")
        of.write("\t")
        of.write("DIST")
        of.write("\t")
        of.write("SCORE")
        of.write("\t")
        of.write("COORD DIST")
        of.write("\t")
        of.write("MATCH COUNT")
        of.write("\n")


if __name__ == "__main__":
    # Example:
    # % peptide_mapping -name S -md 12 -w 0.1 NL63_S.fasta SARS2.txt Spike_SARS2_to_NL63.txt

    argument_parser = ap.ArgumentParser(prog='peptide_mapping',
                                        description='transposition/mapping of peptide sequences')

    argument_parser.add_argument(dest='in_seq', help='fasta target sequence file (example \'Spike_NL63.fasta\')',
                                 type=str)

    argument_parser.add_argument(dest='in_tab',
                                 help='table listing query peptides, tab separated (example \'SARS2.txt\')',
                                 type=str)

    argument_parser.add_argument(dest='out_file', help='out file (example \'Spike_SARS2_to_NL63.txt\')',
                                 type=str)

    argument_parser.add_argument('-md', dest='maximal_distance', help='maximal distance (default: 12)', type=int,
                                 default=12)

    argument_parser.add_argument('-w', dest='coord_distance_weight', help='coordinate distance weight (default: 0.1)',
                                 type=float,
                                 default=0.1)

    argument_parser.add_argument('-verbose', dest='verbose', help='verbose', type=bool,
                                 default=False)

    argument_parser.add_argument('-name', dest='protein_name', help='protein name (example \'S\')', type=str,
                                 required=True
                                 )
    argument_parser.add_argument('--version', action='version', version='%(prog)s ' + PeptideMapping.VERSION)

    args = argument_parser.parse_args()

    in_seq = args.in_seq
    in_tab = args.in_tab
    out_file = args.out_file
    md = args.maximal_distance
    pn = args.protein_name
    verb = args.verbose
    coord_dist_weight = args.coord_distance_weight

    print('Version               : ' + PeptideMapping.VERSION)
    print('Protein name          : ' + str(pn))
    print('Maximal distance      : ' + str(md))
    print('Coord. distance weight: ' + str(coord_dist_weight))
    print('Target sequence file  : ' + str(in_seq))
    print('Query peptides file   : ' + str(in_tab))
    print('Output                : ' + str(out_file))
    print()

    PeptideMapping.transpose(in_tab, in_seq, pn, out_file, md, coord_dist_weight, verb)
