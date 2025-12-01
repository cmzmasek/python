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

import molseq


class CleanFasta(object):
    VERSION = '1.0.3'

    ID_RE = re.compile(">\\s*(.+)")
    GAP_RE = re.compile("[-\\s]+")
    WS_RE = re.compile("\\s+")
    COMBINE_WHITESPACE_RE = re.compile(r"\s+")

    @staticmethod
    def stream_fasta(in_stream, remove_gaps=False):
        seq_id = None
        seq = []
        for line in in_stream:
            line = line.strip()
            if line:
                if line.startswith(">"):
                    if seq_id:
                        yield molseq.MolSeq(seq_id, "".join(seq))
                    seq_id = CleanFasta.ID_RE.search(line).group(1)
                    seq = []
                else:
                    if remove_gaps:
                        line = re.sub(CleanFasta.GAP_RE, "", line)
                    else:
                        line = re.sub(CleanFasta.WS_RE, "", line)
                    seq.append(line)
        if seq_id:
            yield molseq.MolSeq(seq_id, "".join(seq))

    @staticmethod
    def clean_mol_seqs(infile, outfile, min_length, min_ratio, aa, unique_ids, max_length=-1):

        if os.path.isfile(outfile):
            print(outfile + ' already exists')
            sys.exit()

        f0 = open(infile)
        f1 = open(outfile, 'w')

        total = 0
        ignored_irr_chars = 0
        ignored_length_too_short = 0
        ignored_length_too_long = 0
        ignored_name = 0
        ignored_numbers = 0
        ignored_identical_id = 0
        passed = 0
        ids = set()

        for seq in CleanFasta.stream_fasta(f0, True):
            total += 1
            seq_name = seq.get_seq_id().strip()
            if len(seq_name) > 0:
                if re.search(r'\d', seq.get_seq()):
                    ignored_numbers += 1
                else:
                    length = seq.get_length()
                    if length < min_length:
                        ignored_length_too_short += 1
                    elif min_length < max_length < length:
                        ignored_length_too_long += 1
                    else:
                        if aa:
                            reg = length - seq.count_irregular_chars_aa()
                        else:
                            reg = seq.count_regular_chars_na()
                        r = reg / length

                        if r >= min_ratio:
                            nn = CleanFasta.COMBINE_WHITESPACE_RE.sub(' ', seq_name).strip()
                            if unique_ids and nn in ids:
                                ignored_identical_id += 1
                            else:
                                ids.add(nn)
                                seq.set_seq_id(nn)
                                f1.write(seq.to_fasta_wrapped(80))
                                f1.write('\n')
                                passed += 1
                        else:
                            ignored_irr_chars += 1


            else:
                ignored_name += 1
                print('Ignored because empty id:')
                print(str(seq))

        f0.close()
        f1.close()
        print('Version                   : ' + CleanFasta.VERSION)
        print('Input                     : ' + str(total))
        print('Ignored no name           : ' + str(ignored_name))
        print('Ignored numbers in seq    : ' + str(ignored_numbers))
        print('Ignored length (too short): ' + str(ignored_length_too_short))
        print('Ignored length (too long) : ' + str(ignored_length_too_long))
        print('Ignored irreg chars       : ' + str(ignored_irr_chars))
        print('Ignored identical ids     : ' + str(ignored_identical_id))
        print('Passed                    : ' + str(passed))
        print('Wrote to                  : ' + str(outfile))


if __name__ == "__main__":
    argument_parser = ap.ArgumentParser(prog='clean_fasta')

    argument_parser.add_argument(dest='in_file', help='fasta in file (example \'sequences.fasta\')',
                                 type=str)

    argument_parser.add_argument(dest='out_file', help='fasta out file (example \'sequences_29400_09999.fasta\')',
                                 type=str)

    argument_parser.add_argument('-ml', dest='minimal_length', help='minimal length (default: 20)', type=int,
                                 default=20)

    argument_parser.add_argument('-mal', dest='maximal_length', help='maximal length', type=int,
                                 default=-1)

    argument_parser.add_argument('-r', dest='ratio', help='valid char ratio(default: 0.99)', type=float, default=0.99)

    argument_parser.add_argument('-t', dest='type', help='aa or na (default: aa)', type=str, default='aa')

    argument_parser.add_argument('-u', dest='unique_ids', help='t or f (default: t)', type=str, default='t')

    argument_parser.add_argument('--version', action='version', version='%(prog)s ' + CleanFasta.VERSION)

    args = argument_parser.parse_args()

    in_file = args.in_file
    out_file = args.out_file
    ml = args.minimal_length
    mal = args.maximal_length
    ra = args.ratio
    t_str = args.type
    u_str = args.unique_ids

    amino = True
    if t_str.lower() == 'na':
        amino = False

    u_id = True
    if u_str.lower() == 'f':
        u_id = False

    print('Minimal length        : ' + str(ml))
    if mal > ml:
        print('Maximal length        : ' + str(mal))
    print('Valid char ratio      : ' + str(ra))
    if amino:
        print('Type                  : AA')
    else:
        print('Type                  : NA')
    if u_id:
        print('Only unique ids       : yes')
    else:
        print('Only unique ids        : no (allow non-unique sequence names')
    print()

    CleanFasta.clean_mol_seqs(
        in_file,
        out_file,
        ml,
        ra,
        amino,
        u_id,
        mal
    )
