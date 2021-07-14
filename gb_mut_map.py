import re

RE_1 = re.compile("(.+?)\\.\d+\s+(.+?)\\.\d+")


def parse(genome_to_prot_map_file, prot_to_mut_map_file, outfile):
    prot_to_genome_dict = {}
    with open(genome_to_prot_map_file, 'r') as in_file:
        for line in in_file:
            line = line.strip()
            m1 = RE_1.match(line)
            if m1:
                prot_to_genome_dict[m1.group(2)] = m1.group(1)

    out_file = open(outfile, 'w')
    with open(prot_to_mut_map_file, 'r') as in_file:
        counter = 0
        for line in in_file:
            line = line.strip()
            line_split = line.split('\t')
            if len(line_split) > 6:
                prot = line_split[0]
                mut = line_split[6]
                if prot in prot_to_genome_dict:
                    print(str(counter) + ': ' + prot + '->' + prot_to_genome_dict[prot] + '->' + mut)
                    out_file.write(prot_to_genome_dict[prot] + '\t' + mut)
                    out_file.write('\n')
                    counter += 1

    out_file.close()


if __name__ == '__main__':
    parse('genbank_features_map_i2.txt', 'spike.var2.genbank', 'spike_var_genbank_genome_to_mut_i2.txt')
