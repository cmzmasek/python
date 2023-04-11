import re
import sys

import molseq
import csv


class SeqAnnotate(object):
    ID_RE = re.compile(">\s*(.+)")
    GAP_RE = re.compile("[-\s]+")
    WS_RE = re.compile("[\s]+")
    ACCN_RE = re.compile("accn\|(\S+)")

    @staticmethod
    def stream_fasta(in_stream, remove_gaps=False):
        seq_id = None
        seq = []
        for line in in_stream:
            line = line.strip()
            if line:
                if line.startswith(">"):
                    if seq_id:
                        yield molseq.MolSeq(seq_id, "".join(seq))
                    seq_id = SeqAnnotate.ID_RE.search(line).group(1)
                    seq = []
                else:
                    if remove_gaps:
                        line = re.sub(SeqAnnotate.GAP_RE, "", line)
                    else:
                        line = re.sub(SeqAnnotate.WS_RE, "", line)
                    seq.append(line)
        if seq_id:
            yield molseq.MolSeq(seq_id, "".join(seq))

    @staticmethod
    def run(seq_file, annotation_file, out_file, min_length, min_ratio):

        annotations = {}
        with open(annotation_file) as af:
            reader = csv.reader(af)
            for row in reader:
                annotations[row[43]] = row

        f0 = open(seq_file)
        f1 = open(out_file, 'w')
        total = 0
        ignored_irr_chars = 0
        ignored_length = 0
        kept = 0
        wrong_segment = 0

        for mol_seq in SeqAnnotate.stream_fasta(f0, True):
            total += 1
            reg = mol_seq.count_regular_chars_na()
            length = mol_seq.get_length()
            r = reg / length
            if length >= min_length:
                if r >= min_ratio:
                    gb_id = SeqAnnotate.ACCN_RE.search(mol_seq.get_seq_id()).group(1)
                    if gb_id in annotations:
                        a = annotations[gb_id]
                        if a[20] == 'L':
                            kept += 1
                            # Genome Id           0
                            # Genome              1
                            # Species            12
                            # Strain             15
                            # Segment            20
                            # GB acc             43
                            # Collection Year    68
                            # Isolation Country  70
                            # HostName           74
                            mol_seq.set_seq_id(
                                a[1] + '|' + a[12] + '|' + a[15] + '|' + a[20] + '|' + a[43] + '|' + a[68] + '|' + a[
                                    70] + '|' + a[74]);
                            f1.write(mol_seq.to_fasta_wrapped(80))
                            f1.write('\n')
                        else:
                            wrong_segment += 1
                    else:
                        print('Error: Annotation for ' + gb_id + ' not found')
                        sys.exit(-1)
                else:
                    ignored_irr_chars += 1
            else:
                ignored_length += 1

        f0.close()
        f1.close()

        print('Genomes: Total               : ' + str(total))
        print('Genomes: Ignored Length      : ' + str(ignored_length))
        print('Genomes: Ignored Irreg Chars : ' + str(ignored_irr_chars))
        print('Not L segment                : ' + str(wrong_segment))
        print('Genomes: Kept                : ' + str(kept))
        print()


if __name__ == "__main__":
    SeqAnnotate.run(
        '/home/lambda/Dropbox/WORK/JCVI/PHYLOGENETICS/Arenaviridae/BVBRC/Arenaviridae_5417_genomes_DNA.fasta',
        '/home/lambda/Dropbox/WORK/JCVI/PHYLOGENETICS/Arenaviridae/BVBRC/Arenaviridae_5417_genomes_DNA.csv',
        '/home/lambda/Dropbox/WORK/JCVI/PHYLOGENETICS/Arenaviridae/BVBRC/Arenaviridae_5417_genomes_DNA_L_4000_0999',
        4000, 0.999)
