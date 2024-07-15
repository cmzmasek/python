import argparse as ap
import os
import re
import sys

from molseq_class import MolSeq

# using uppercase to indicate constant
LABEL_RE = re.compile(">\\s*(.+)")
GAP_RE = re.compile("[-\\s]+")
WS_RE = re.compile("\\s+")


class MsaRename(object):
    VERSION = '1.0.0'

    @staticmethod
    def rename(msa_file, id_map_file, outfile):
        if os.path.isfile(outfile):
            print(outfile + ' already exists')
            sys.exit(-1)
        id_name_dic = {}
        mapping_counter = 0
        with open(id_map_file, 'r') as f:
            for line in f:
                line = line.strip()
                s = line.split('\t')
                #s = line.split(' ', 1)
                if s[0] in id_name_dic:
                    print("Error: Duplicate id: " + s[0])
                    sys.exit(-1)
                if s[1] in id_name_dic.values():
                    print("Error: Duplicate name: " + s[1])
                    #sys.exit(-1)

                id_name_dic[s[0]] = s[1]
                mapping_counter += 1
        print('Read in ' + str(mapping_counter) + ' mappings')
        seqs = MsaRename.parse_fasta_file(msa_file, False)

        counter = 0
        with open(outfile, 'x') as o:
            for seq in seqs:
                if not seq.get_label() in id_name_dic:
                    print("Error: No mapping for: " + seq.get_label())
                    sys.exit(-1)
                seq.set_label(id_name_dic[seq.get_label()])
                o.write(seq.to_fasta())
                o.write('\n')
                counter += 1
        print('Mapped ' + str(counter) + ' names')

    @staticmethod
    def parse_fasta_file(file, remove_gaps=False):
        with open(file, 'r') as infile:
            # memorize what has been read in current_label and current_seq
            current_label = ''
            current_seq = ''
            mol_seq_list = []

            for line in infile:
                line = line.strip()

                if line:
                    if line.startswith('>'):
                        # if this is a new seq but not the 1st seq, create a MolSeq object of the previous seq
                        if current_seq:
                            mol_seq_list.append(MolSeq(current_label, current_seq))

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
            mol_seq_list.append(MolSeq(current_label, current_seq))

            return mol_seq_list


if __name__ == '__main__':
    argument_parser = ap.ArgumentParser(prog='msa_rename')

    argument_parser.add_argument(dest='fasta_file', help='fasta file (example \'mafft.fasta\')',
                                 type=str)

    argument_parser.add_argument(dest='map_file', help='map file (example \'sars2.nim\')',
                                 type=str)

    argument_parser.add_argument(dest='out_file', help='output file (example \'mafft_fn.fasta\')',
                                 type=str)
    argument_parser.add_argument('--version', action='version', version='%(prog)s ' + MsaRename.VERSION)

    args = argument_parser.parse_args()

    in_file = args.fasta_file
    map_file = args.map_file
    out_file = args.out_file

    MsaRename.rename(in_file, map_file, out_file)
