import argparse as ap
import os
import re
import sys


class LinExtract(object):
    VERSION = '1.0.0'

    @staticmethod
    def extract(infile, outfile):

        if os.path.isfile(outfile):
            print(outfile + ' already exists')
            sys.exit()

        f1 = open(outfile, 'w')

        with open(infile, 'r') as file:
            saw_org = False
            for line in file:
                if saw_org:
                    if line.startswith(" "):
                        line = re.sub(r';\s*', '\t', line.strip())
                        if line.endswith('.'):
                            line = line[:-1]
                            f1.write(line + '\n')
                            saw_org = False
                        else:
                            f1.write(line)
                    else:
                        saw_org = False
                elif line.startswith("VERSION"):
                    x = line.split()
                    f1.write(x[1].strip())
                    f1.write("\t")
                elif line.startswith("  ORGANISM"):
                    saw_org = True
                    tax = line[12:]
                    f1.write(tax.strip())
                    f1.write("\t")

        f1.close()


if __name__ == "__main__":
    argument_parser = ap.ArgumentParser(prog='lin_extract')

    argument_parser.add_argument(dest='in_file',
                                 type=str)

    argument_parser.add_argument(dest='out_file', type=str)

    args = argument_parser.parse_args()

    in_file = args.in_file
    out_file = args.out_file

    LinExtract.extract(
        in_file,
        out_file)
