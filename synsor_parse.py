import argparse as ap
import os
import sys
import re


class SynsorParse(object):
    VERSION = '1.0.0'

    @staticmethod
    def parse(infile, outfile):

        if os.path.isfile(outfile):
            print(outfile + ' already exists')
            sys.exit()

        f1 = open(outfile, 'w')

        total = 0
        not_syn = 0

        with open(infile, 'r') as file:
            for line in file:
                if not line.startswith("seqId"):
                    total += 1
                    x = line.split('\t')
                    seqid = x[0]
                    prob = x[2]
                    pred_class = x[4]
                    if pred_class == "Not synthetic":
                        not_syn += 1
                    seqid = re.sub(r"\s+", "_", seqid)
                    f1.write(seqid + '\t')
                    f1.write(pred_class + '\t')
                    f1.write(prob + '\n')
        f1.close()

        print("Total  : " + str(total))
        print("Not syn: " + str(not_syn))


if __name__ == "__main__":
    argument_parser = ap.ArgumentParser(prog='synsor_parse')

    argument_parser.add_argument(dest='in_file',
                                 type=str)

    argument_parser.add_argument(dest='out_file', type=str)

    args = argument_parser.parse_args()

    in_file = args.in_file
    out_file = args.out_file

    SynsorParse.parse(
        in_file,
        out_file)
