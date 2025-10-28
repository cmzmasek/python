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

import textwrap


class MolSeq(object):
    """A molecular sequence.

        Attributes:
            seq_id: A string representing the sequence identifier or name.
            seq: A string representing of the molecular sequence.
    """

    def __init__(self, seq_id, seq):
        self.__seq_id = str(seq_id).strip()
        self.__seq = str(seq).strip()

    def get_seq_id(self):
        return self.__seq_id

    def set_seq_id(self, id):
        self.__seq_id = id

    def get_seq(self):
        return self.__seq

    def get_length(self):
        return len(self.__seq)

    def to_fasta(self):
        return ">{}\n{}".format(self.get_seq_id(), self.get_seq())

    def to_fasta_wrapped(self, w):
        x = textwrap.fill(self.get_seq(), width=w)
        return ">{}\n{}".format(self.get_seq_id(), x)

    def count_regular_chars_na(self):
        a = self.get_seq().count('a') + self.get_seq().count('A')
        c = self.get_seq().count('c') + self.get_seq().count('C')
        g = self.get_seq().count('g') + self.get_seq().count('G')
        t = self.get_seq().count('t') + self.get_seq().count('T')
        return a + c + g + t

    def count_irregular_chars_aa(self):
        return (self.get_seq().count('_') + self.get_seq().count('-') +
                self.get_seq().count('?') + self.get_seq().count('X') + self.get_seq().count('x') +
                self.get_seq().count('*') + self.get_seq().count('.'))

    def __str__(self):
        return self.to_fasta_wrapped(60)

    def __len__(self):
        return self.get_length()

    def __repr__(self):
        return "{}:{}:{}".format(self.__class__.__name__, self.get_seq_id(), self.get_seq())


if __name__ == "__main__":
    seq0 = MolSeq(" abcd ", " acgttgtca")

    print(seq0.to_fasta())
    print(seq0.get_length())
    print(str(seq0))
    print(repr(seq0))
