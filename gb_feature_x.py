import re

FEATURE_RE = re.compile(">\s*Feature\s+(?:gb|dbj)\|(.+)\|", re.IGNORECASE)
TARGET_RE = re.compile("product\s+surface\s+glycoprotein", re.IGNORECASE)
PROTEIN_ID_RE = re.compile("protein_id\s+(?:gb|dbj)\|(.+)\|", re.IGNORECASE)


def parse(infile, outfile):
    out_file = open(outfile, 'w')
    with open(infile, 'r') as in_file:
        current_feature = None
        saw_target = False
        counter = 0
        for line in in_file:
            line = line.strip()
            if line:
                feat_m = FEATURE_RE.match(line)
                if feat_m:
                    current_feature = feat_m.group(1)
                    saw_target = False
                else:
                    if current_feature:
                        target_m = TARGET_RE.match(line)
                        if target_m:
                            saw_target = True
                        else:
                            protein_m = PROTEIN_ID_RE.match(line)
                            if saw_target is True and protein_m:
                                protein_id = protein_m.group(1)
                                print(str(counter))
                                print(current_feature + ' ' + protein_id)
                                out_file.write(current_feature + '\t' + protein_id)
                                out_file.write('\n')
                                counter += 1
                                saw_target = False
                                current_feature = None

    out_file.close()


if __name__ == '__main__':
    parse('sequence_i2.txt', 'genbank_features_map_i2.txt')
