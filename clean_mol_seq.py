import re
import molseq


class CleanMolSeq(object):
    ID_RE = re.compile(">\s*(.+)")
    GAP_RE = re.compile("[-\s]+")
    WS_RE = re.compile("[\s]+")

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
                    seq_id = CleanMolSeq.ID_RE.search(line).group(1)
                    seq = []
                else:
                    if remove_gaps:
                        line = re.sub(CleanMolSeq.GAP_RE, "", line)
                    else:
                        line = re.sub(CleanMolSeq.WS_RE, "", line)
                    seq.append(line)
        if seq_id:
            yield molseq.MolSeq(seq_id, "".join(seq))

    @staticmethod
    def read_protein_fasta_file(infile, outfile, min_length):
        genomes = set()
        ignored_length = 0
        ignored_irr = 0
        ignored_name = 0
        total = 0
        f0 = open(infile)
        f1 = open(outfile, 'w')
        for mol_seq in CleanMolSeq.stream_fasta(f0, True):
            total += 1
            length = mol_seq.get_length()
            if length >= min_length:
                irr = mol_seq.count_irregular_chars_aa()
                if irr < 1:
                    name = mol_seq.get_seq_id()
                    s = name.split('|')
                    if len(s) > 2:
                        genome_acc = s[1]
                        protein_acc = s[2]
                        f1.write(genome_acc + '\t' + protein_acc + '\n')
                        genomes.add(genome_acc)
                    else:
                        ignored_name += 1
                else:
                    ignored_irr += 1
            else:
                ignored_length += 1
        f0.close()
        f1.close()
        print('Protein: Min Length               : ' + str(min_length))
        print('Protein: Total                    : ' + str(total))
        print('Protein: Ignored Length           : ' + str(ignored_length))
        print('Protein: Irreg Chars              : ' + str(ignored_irr))
        print('Protein: Ignored Ill-Formated Name: ' + str(ignored_name))
        print('Protein: Returned                 : ' + str(len(genomes)))

        return genomes

    @staticmethod
    def extract_from_protein_fasta_file(infile, outfile, keep_proteins_genome_acc):
        f0 = open(infile)
        f1 = open(outfile, 'w')
        total = 0
        for mol_seq in CleanMolSeq.stream_fasta(f0, True):
            name = mol_seq.get_seq_id()
            s = name.split('|')
            genome_acc = s[1]
            if genome_acc in keep_proteins_genome_acc:
                total += 1
                f1.write(mol_seq.to_fasta_wrapped(80))
                f1.write('\n')
        f0.close()
        f1.close()
        print('Protein: Sequences stored         : ' + str(total))

    @staticmethod
    def clean_mol_seqs(infile, outfile, min_length, min_ratio, verbose, protein_file=None, protein_outfile=None,
                       protein_outfile_fasta=None,
                       protein_min_length=0):
        genomes = None
        if protein_file:
            genomes = CleanMolSeq.read_protein_fasta_file(protein_file, protein_outfile, protein_min_length)
        f0 = open(infile)
        f1 = open(outfile, 'w')
        total = 0
        ignored_irr_chars = 0
        ignored_length = 0
        ignored_name = 0
        ignored_no_protein = 0
        kept = 0
        keep_proteins_genome_acc = set()
        for mol_seq in CleanMolSeq.stream_fasta(f0, True):
            name_lwr = mol_seq.get_seq_id().lower()
            if '|severe_acute_respiratory_syndrome_related_coronavirus' in name_lwr or '2019_ncov' in name_lwr or 'hcov_19' in name_lwr or 'sars_cov_2' in name_lwr or 'sars_cov2' in name_lwr:
                total += 1
                reg = mol_seq.count_regular_chars_na()
                length = mol_seq.get_length()
                r = reg / length

                if length >= min_length:
                    if r >= min_ratio:
                        if genomes:
                            s = mol_seq.get_seq_id().split('|')
                            genome_acc = s[1]
                            if genome_acc not in genomes:
                                ignored_no_protein += 1
                                continue
                            keep_proteins_genome_acc.add(genome_acc)
                        kept += 1
                        f1.write(mol_seq.to_fasta_wrapped(80))
                        f1.write('\n')
                    else:
                        ignored_irr_chars += 1
                else:
                    ignored_length += 1
            else:
                ignored_name += 1
        f0.close()
        f1.close()
        if genomes:
            CleanMolSeq.extract_from_protein_fasta_file(protein_file, protein_outfile_fasta, keep_proteins_genome_acc)

        print('Genomes: Ignored Name        : ' + str(ignored_name))
        print('Genomes: Total (Correct Name): ' + str(total))
        print('Genomes: Ignored Length      : ' + str(ignored_length))
        print('Genomes: Ignored Irreg Chars : ' + str(ignored_irr_chars))
        if genomes:
            print('Genomes: Ignored No Protein  : ' + str(ignored_no_protein))
        print('Genomes: Kept                : ' + str(kept))
        print()


if __name__ == "__main__":
    CleanMolSeq.clean_mol_seqs(
        'C:/Users/czmasek/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_4_JUNE_21/SARS2.fasta',
        'C:/Users/czmasek/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_4_JUNE_21/SARS2_29400_09999.fasta',
        29400,
        0.9999,
        True,
        'C:/Users/czmasek/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_4_JUNE_21/ProteinFastaResults_S.fasta',
        'C:/Users/czmasek/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_4_JUNE_21/SARS2_29400_09999_prot_out.txt',
        'C:/Users/czmasek/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_4_JUNE_21/SARS2_29400_09999_prot_out.fasta',
        1250)
