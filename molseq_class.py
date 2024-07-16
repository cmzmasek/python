class MolSeq(object):
    """
    for molecular sequences
    :parameter:
        seq_label: a string representing the sequence.
                    Input seq_label could be an integer, which needs to be converted to a str.
        seq: a string representing the molecular sequence
    """

    def __init__(self, label, seq):
        self.__label = str(label).strip()
        self.__seq = str(seq).strip()

    def __str__(self):
        return self.to_fasta(60)

    def __len__(self):
        return self.get_length()

    def to_fasta(self, chars_per_line=60):
        """
        :param chars_per_line: an optional argument to limit numbers of chars per line.
                                False: no limit on chars_per_line
        :return: a string of fasta formatted sequence
        """
        if chars_per_line == False:
            return f">{self.__label}\n{self.__seq}"
        else:
            self.__fasta_line = ''
            for i in range(self.get_length()):
                if (i + 1) % chars_per_line == 0:
                    self.__fasta_line += self.__seq[i] + '\n'
                else:
                    self.__fasta_line += self.__seq[i]
            return f">{self.__label}\n{self.__fasta_line}"

    def get_label(self):
        return self.__label

    def set_label(self, new_label):
        self.__label = new_label

    def get_seq(self):
        return self.__seq

    def get_length(self):
        return len(self.__seq)


if __name__ == '__main__':
    seq1 = MolSeq(111,
                  'abbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddcccabbddddccc')

    print('\nchars_per_line=default')
    print(seq1.to_fasta())
    print('\nchars_per_line=False')
    print(seq1.to_fasta(False))
    print('\nchars_per_line=30')
    print(seq1.to_fasta(30))
    print('\nprint method')
    print(seq1)
    print(seq1.get_length())
    print(len(seq1))
    print('current label:')
    print(seq1.get_label(y='x'))
    print('resetting label to new_label')
    seq1.set_label('new_label')
    print('new label:')
    print(seq1.get_label())
    print('get seq:')
    print(seq1.get_seq())
