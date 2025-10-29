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

    def __init__(self, start, end, score, sequence):
        if end <= start:
            raise ValueError("Start must be smaller than end.")
        self.__start = int(start)
        self.__end = int(end)
        self.__score = int(score)
        self.__sequence = str(sequence).strip()

    def get_start(self):
        return self.__start

    def get_end(self):
        return self.__end

    def get_score(self):
        return self.__score

    def get_sequence(self):
        return self.__sequence

    def calc_length(self):
        return 1 + self.get_end() - self.get_start()

    def to_fasta(self):
        return ">{}\n{}".format(str(self.get_start()) + "-" + str(self.get_end()) + "_" + str(self.get_score()),
                                self.get_sequence())

    def __str__(self):
        return self.to_fasta()

    def __len__(self):
        return self.calc_length()

    def __repr__(self):
        return "{}:{}:{}:{}:{}".format(self.__class__.__name__, self.get_start(), self.get_end(), self.get_sequence(),
                                       self.get_score())


if __name__ == "__main__":
    t = TargetMatch(3, 5, 0.1, "abgty")

    print(t.to_fasta())
    print(t.calc_length())
    print(str(t))
    print(repr(t))
