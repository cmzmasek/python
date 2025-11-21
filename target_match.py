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


class TargetMatch(object):

    # Used by PeptideMapping

    def __init__(self, query_start, query_end, query_sequence, target_start, target_end, target_sequence, distance,
                 score=0.0):
        if target_end <= target_start or query_end <= query_start:
            raise ValueError("Start must be smaller than end.")
        if len(query_sequence) != len(target_sequence):
            raise ValueError("Query and target sequences must have same length.")
        if target_end - target_start != query_end - query_start:
            raise ValueError("Query and target must have same length by coordinates.")
        if target_end - target_start + 1 != len(target_sequence):
            raise ValueError("Target sequence length does not match length by coordinates")
        if query_end - query_start + 1 != len(query_sequence):
            raise ValueError("Query sequence length does not match length by coordinates")
        self.__query_start = int(query_start)
        self.__query_end = int(query_end)
        self.__target_start = int(target_start)
        self.__target_end = int(target_end)
        self.__distance = float(distance)
        self.__score = float(score)
        self.__query_sequence = str(query_sequence)
        self.__target_sequence = str(target_sequence)

    def get_query_start(self):
        return self.__query_start

    def get_query_end(self):
        return self.__query_end

    def get_target_start(self):
        return self.__target_start

    def get_target_end(self):
        return self.__target_end

    def get_distance(self):
        return self.__distance

    def get_score(self):
        return self.__score

    def get_target_sequence(self):
        return self.__target_sequence

    def get_query_sequence(self):
        return self.__target_sequence

    def calc_length(self):
        return 1 + self.get_query_end() - self.get_query_start()

    def to_fasta(self):
        return ">{}\n{}".format(
            str(self.get_target_start()) + "-" + str(self.get_target_end()) + "_" + str(self.get_distance()),
            self.get_target_sequence())

    def to_aln(self):
        return "{}\n{}\n{}\n{}".format(
            str(self.get_query_start()) + "-" + str(self.get_query_end()),
            str(self.get_query_sequence()),
            str(self.get_target_sequence()),
            str(self.get_target_start()) + "-" + str(self.get_target_end())
        )

    def __lt__(self, other):
        return self.get_score() < other.get_score()

    def __str__(self):
        return self.to_fasta()

    def __len__(self):
        return self.calc_length()

    def __repr__(self):
        return "{}:{}:{}:{}:{}:{}:{}:{}".format(self.__class__.__name__, self.get_query_start(), self.get_query_end(),
                                                self.get_target_start(), self.get_target_end(),
                                                self.get_query_sequence(), self.get_target_sequence(),
                                                self.get_distance())


if __name__ == "__main__":
    t = TargetMatch(3, 7, "abgty", 13, 17, "abxxx", 11, 0.1)

    print(t.to_fasta())
    print(t.get_score())
    print(t.calc_length())
    print(str(t))
    print(repr(t))
    print(t.to_aln())

    t2 = TargetMatch(3, 7, "abgty", 13, 17, "abxxx", 11, 0.2)

    print(t < t2)
