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

import target_match


class PeptideMapping(object):
    VERSION = '0.0.0'

    ID_RE = re.compile(">\\s*(.+)")
    GAP_RE = re.compile("[-\\s]+")
    WS_RE = re.compile("\\s+")
    COMBINE_WHITESPACE_RE = re.compile(r"\s+")

    @staticmethod
    def map(query, target):
        best_score = len(query) + 1
        best_matches = None
        query_len = len(query)
        for i in range(len(target) - query_len + 1):
            t = target[i:i + query_len]
            d = PeptideMapping.hamming_distance(t, query)
            # print(t + ": " + str(d))
            if d <= best_score:
                if d < best_score:
                    best_matches = []
                    best_score = d
                best_matches.append(target_match.TargetMatch(i, i + query_len - 1, d, t))

        return best_matches

    @staticmethod
    def map_windowed(query, target, window_start, window_end):
        return map(query, target[window_start: window_end + 1])

    @staticmethod
    def hamming_distance(s1, s2):
        """Calculate the Hamming distance between two strings."""
        if len(s1) != len(s2):
            raise ValueError("The strings must have the same length.")
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    @staticmethod
    def inverted_normalized_hamming_distance(s1, s2):
        return 1 - PeptideMapping.hamming_distance(s1, s2) / len(s1)


if __name__ == "__main__":
    print(PeptideMapping.inverted_normalized_hamming_distance("abcd", "axyx"))

    print(PeptideMapping.inverted_normalized_hamming_distance("abcd", "abcd"))

    m = PeptideMapping.map("query", "queryquery")
    print(m)

    m = PeptideMapping.map("abc", "abccbacbcbabc")
    print(m)

    m = PeptideMapping.map("abcd", "mnapnmpnmnpnonappn")
    print(m)
