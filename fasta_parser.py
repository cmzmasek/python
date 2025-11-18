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

import re

from molseq_class import MolSeq

# using uppercase to indicate constant
LABEL_RE = re.compile(">\s*(.+)")
GAP_RE = re.compile("[-\s]+")
WS_RE = re.compile("[\s]+")


def parse_fasta_file(file, remove_gaps=False):
    with open(file, 'r') as infile:
        # memorize what has been read in current_label and current_seq
        current_label = ''
        current_seq = ''
        MolSeq_list = []

        for line in infile:
            line = line.strip()

            if line:
                if line.startswith('>'):
                    # if this is a new seq but not the 1st seq, create a MolSeq object of the previous seq
                    if current_seq:
                        MolSeq_list.append(MolSeq(current_label, current_seq))

                        # then reset current_label and current_seq
                        current_label = LABEL_RE.search(line).group(1)
                        current_seq = ''

                    # for the 1st seq in the file
                    else:
                        current_label = LABEL_RE.search(line).group(1)

                else:
                    if remove_gaps:
                        line = re.sub(GAP_RE, "", line)
                    else:
                        line = re.sub(WS_RE, "", line)
                    current_seq += line

        # for the last seq in the file
        MolSeq_list.append(MolSeq(current_label, current_seq))

        return MolSeq_list


if __name__ == '__main__':
    MolSeq_list1 = parse_fasta_file('test.fasta')
    for seq in MolSeq_list1:
        print(seq.to_fasta())

    MolSeq_list2 = parse_fasta_file('test.fasta', remove_gaps=True)
    for seq in MolSeq_list2:
        print(seq.to_fasta())
