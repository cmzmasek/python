import re
import molseq


class GbProteinX(object):
    ID_RE = re.compile(">\s*(.+)")  # using uppercase to indicate constant
    GAP_RE = re.compile("[-\s]+")
    WS_RE = re.compile("[\s]+")

    NAME_RE = re.compile('lcl\|([A-Z0-9]+)\.\d+_prot')

    @staticmethod
    def stream_fasta(in_stream, remove_gaps=False):
        seq_id = None
        seq = []
        for line in in_stream:
            line = line.strip()
            if line:  # empty strings are false
                if line.startswith(">"):
                    if seq_id:
                        yield molseq.MolSeq(seq_id, "".join(seq))
                    seq_id = GbProteinX.ID_RE.search(line).group(1)
                    seq = []
                else:
                    if remove_gaps:
                        line = re.sub(GbProteinX.GAP_RE, "", line)
                    else:
                        line = re.sub(GbProteinX.WS_RE, "", line)
                    seq.append(line)
        if seq_id:
            yield molseq.MolSeq(seq_id, "".join(seq))

    @staticmethod
    def extract(infile, outfile, verbose=True):
        f0 = open(infile)
        f1 = open(outfile, 'w')
        total = 0
        ignored = 0
        kept = 0
        for mol_seq in GbProteinX.stream_fasta(f0, True):
            total += 1
            name = mol_seq.get_seq_id()
            # if '[gene=S]' in name or 'surface glycoprotein' in name:
            if 'gene=S' in name:
                kept += 1
                irreg = mol_seq.calculate_irregular_chars_aa()
                length = mol_seq.get_length()
                if irreg > 0 or length != 1273:
                    print(name)
                    print('Length         : ' + str(length))
                    print('Irregular chars: ' + str(irreg))
                    continue
                m1 = GbProteinX.NAME_RE.match(name)
                if m1:
                    mol_seq.set_seq_id(m1.group(1))
                else:
                    print(name)
                    exit()

                f1.write(mol_seq.to_fasta_wrapped(80))
                f1.write('\n')
            else:
                ignored += 1
        f0.close()
        f1.close()
        print('Total  : ' + str(total))
        print('Kept   : ' + str(kept))
        print('Ignored: ' + str(ignored))
        print()


if __name__ == "__main__":
    GbProteinX.extract('sequence_8.fasta', 'sequence_8_S.fasta', True)
