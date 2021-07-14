import re

RE_1 = re.compile("\[ncbi\]:(.+)")


def parse(prot_to_mut_map_file, outfile, max):
    out_file = open(outfile, 'w')
    with open(prot_to_mut_map_file, 'r') as in_file:
        counter = 0
        for line in in_file:
            line = line.strip()
            m1 = RE_1.match(line)
            if m1:
                gb_id = m1.group(1)
                if gb_id != 'Isolate':
                    counter += 1
                    out_file.write(gb_id)
                    if counter == max:
                        out_file.write('\n\n\n')
                        counter = 0
                    else:
                        out_file.write(', ')

    out_file.close()


if __name__ == '__main__':
    parse('ext_node_ids.txt', 'ext_node_ids_reorg.txt', 200)
