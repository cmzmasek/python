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
    def perform_mapping_windowed(query, target, query_start, query_end, window_start, window_end):
        if window_end > len(target) - 1 or window_start < 0:
            raise ValueError("Start or end are out of bounds.")
        if window_end <= window_start:
            raise ValueError("Start must be smaller than end.")

        return perform_mapping(query, target[window_start: window_end + 1], query_start, query_end)

    @staticmethod
    def hamming_distance(s1, s2):
        """Calculate the Hamming distance between two strings."""
        if len(s1) != len(s2):
            raise ValueError("The strings must have the same length.")
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    @staticmethod
    def inverted_normalized_hamming_distance(s1, s2):
        return 1 - PeptideMapping.hamming_distance(s1, s2) / len(s1)


def perform_mapping(query, target, query_start, query_end):
    if len(target) <= len(query):
        raise ValueError("Query must be shorter than target.")
    best_score = len(query) + 1
    best_matches = None
    query_len = len(query)
    for i in range(len(target) - query_len + 1):
        t = target[i:i + query_len]
        dist = PeptideMapping.hamming_distance(t, query)
        # print(t + ": " + str(d))
        if dist <= best_score:
            if dist < best_score:
                best_matches = []
                best_score = dist
            best_matches.append(target_match.TargetMatch(query_start, query_end, i, i + query_len - 1, dist, query, t))

    return best_matches


if __name__ == "__main__":
    print(PeptideMapping.inverted_normalized_hamming_distance("abcd", "axyx"))

    print(PeptideMapping.inverted_normalized_hamming_distance("abcd", "abcd"))

    m = perform_mapping("query", "queryquery", 0, 4)
    for x in m:
        print(str(x.get_score()) + ":")
        print(x.to_aln())

    m = perform_mapping("abc", "abccbacbcbabc", 0, 2)
    for x in m:
        print(str(x.get_score()) + ":")
        print(x.to_aln())

    m = perform_mapping("abcd", "mnapnmpnmnpnonappn", 0, 3)

    for x in m:
        print(str(x.get_score()) + ":")
        print(x.to_aln())

    m = PeptideMapping.perform_mapping_windowed("abcd", "mnapnmpnmnpnonappn", 0, 3, 0, 10)

    for x in m:
        print(str(x.get_score()) + ":")
        print(x.to_aln())
