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
        irreg_count = dict()
        i0 = 0
        i50 = 0
        i100 = 0
        i150 = 0
        i200 = 0
        i250 = 0
        i300 = 0
        i350 = 0
        i400 = 0
        i450 = 0
        i500 = 0
        i550 = 0
        i600 = 0
        i650 = 0
        i700 = 0
        i750 = 0
        i800 = 0
        i850 = 0
        i900 = 0
        i950 = 0
        i1000 = 0
        i1000p = 0

        for mol_seq in CleanMolSeq.stream_fasta(f0, True):
            name_lwr = mol_seq.get_seq_id().lower()
            if '|severe_acute_respiratory_syndrome_related_coronavirus' in name_lwr or '2019_ncov' in name_lwr or 'hcov_19' in name_lwr or 'sars_cov_2' in name_lwr or 'sars_cov2' in name_lwr:
                total += 1
                reg = mol_seq.count_regular_chars_na()
                length = mol_seq.get_length()
                irreg = length - reg
                r = reg / length
                if length >= min_length:
                    if irreg in irreg_count:
                        irreg_count[irreg] += 1
                    else:
                        irreg_count[irreg] = 1

                    if irreg == 0:
                        i0 += 1
                    elif irreg < 50:
                        i50 += 1
                    elif irreg < 100:
                        i100 += 1
                    elif irreg < 150:
                        i150 += 1
                    elif irreg < 200:
                        i200 += 1
                    elif irreg < 250:
                        i250 += 1
                    elif irreg < 300:
                        i300 += 1
                    elif irreg < 350:
                        i350 += 1
                    elif irreg < 400:
                        i400 += 1
                    elif irreg < 450:
                        i450 += 1
                    elif irreg < 500:
                        i500 += 1
                    elif irreg < 550:
                        i550 += 1
                    elif irreg < 600:
                        i600 += 1
                    elif irreg < 650:
                        i650 += 1
                    elif irreg < 700:
                        i700 += 1
                    elif irreg < 750:
                        i750 += 1
                    elif irreg < 800:
                        i800 += 1
                    elif irreg < 850:
                        i850 += 1
                    elif irreg < 900:
                        i900 += 1
                    elif irreg < 950:
                        i950 += 1
                    elif irreg < 1000:
                        i1000 += 1
                    else:
                        i1000p += 1

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

        print(str(i0))
        print(str(i50))
        print(str(i100))
        print(str(i150))
        print(str(i200))
        print(str(i250))
        print(str(i300))
        print(str(i350))
        print(str(i400))
        print(str(i450))
        print(str(i500))
        print(str(i550))
        print(str(i600))
        print(str(i650))
        print(str(i700))
        print(str(i750))
        print(str(i800))
        print(str(i850))
        print(str(i900))
        print(str(i950))
        print(str(i1000))
        print(str(i1000p))

        f0.close()
        f1.close()
        if genomes:
            CleanMolSeq.extract_from_protein_fasta_file(protein_file, protein_outfile_fasta,
                                                        keep_proteins_genome_acc)

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
       '/home/lambda/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_6_AUG_21/SARS2.fasta',
       '/home/lambda/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_6_AUG_21/SARS2_29400_09999.fasta',
       29400,
       0.9999,
       True,
       '/home/lambda/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_6_AUG_21/ProteinFastaResults_S.fasta',
       '/home/lambda/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_6_AUG_21/SARS2_29400_09999_prot_out.txt',
       '/home/lambda/Dropbox/WORK/JCVI/SARS_COV_2_REF_TREE/I_6_AUG_21/SARS2_29400_09999_prot_out.fasta',
       1250)

