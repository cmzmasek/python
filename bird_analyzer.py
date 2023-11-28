import os.path
import sys

import pandas as pd


# Last modified 2023/11/27
# Roshni Bhattacharya
# Christian M. Zmasek

class BirdAnalyzer(object):
    VERSION = '0.0.9'

    DEBUG = False
    DEBUG_2 = False

    TAXONOMY_FILE_SCIENTIFIC_NAME = 'scientific name'
    TAXONOMY_FILE_COMMON_NAME = 'English name'
    TAXONOMY_FILE_ORDER = 'order'
    TAXONOMY_FILE_FAMILY = 'family'

    TAXONOMY_MAP_FILE_ORDER = 'ORDER'
    TAXONOMY_MAP_FILE_FAMILY = 'FAMILY'
    TAXONOMY_MAP_FILE_GENERA = 'GENERA'
    TAXONOMY_MAP_FILE_FEATURE = 'FEATURE'

    ANNOTATION_FILE_HOST_NAME = 'Host Name'
    ANNOTATION_FILE_COMMON_HOST_NAME = 'Host Common Name'
    ANNOTATION_FILE_GB_ACC = 'GenBank Accessions'
    ANNOTATION_FILE_GENOME_ID = 'Genome ID'

    NAN = 'nan'

    OTHER_MAMMALS = ['dog', 'horse']
    SWINE = ['pig', 'swine', 'sus scrofa domesticus', 'sus scrofa']

    @staticmethod
    def read_taxonomy_file(file_name):
        source_df = pd.read_csv(file_name, encoding='unicode_escape', dtype={
            BirdAnalyzer.TAXONOMY_FILE_SCIENTIFIC_NAME: str,
            BirdAnalyzer.TAXONOMY_FILE_COMMON_NAME: str,
            BirdAnalyzer.TAXONOMY_FILE_ORDER: str,
            BirdAnalyzer.TAXONOMY_FILE_FAMILY: str
        })
        return source_df

    @staticmethod
    def read_taxonomy_to_feature_map_file(file_name):
        source_df = pd.read_csv(file_name, encoding='unicode_escape', dtype={
            BirdAnalyzer.TAXONOMY_MAP_FILE_ORDER: str,
            BirdAnalyzer.TAXONOMY_MAP_FILE_FAMILY: str,
            BirdAnalyzer.TAXONOMY_MAP_FILE_GENERA: str,
            BirdAnalyzer.TAXONOMY_MAP_FILE_FEATURE: str
        })
        return source_df

    @staticmethod
    def make_taxonomy_to_feature_map(taxonomy_master_file_name, taxonomy_to_feature_map_file_name):
        taxonomy_master_table_df = BirdAnalyzer.read_taxonomy_file(taxonomy_master_file_name)
        taxonomy_to_feature_map_df = BirdAnalyzer.read_taxonomy_to_feature_map_file(taxonomy_to_feature_map_file_name)

        order_to_feature = {}
        family_to_feature = {}
        genus_to_feature = {}

        for index, row in taxonomy_to_feature_map_df.iterrows():

            order = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_ORDER]).lower().strip()
            family = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_FAMILY]).lower().strip()
            genera = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_GENERA]).lower().strip()
            feature = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_FEATURE]).lower().strip()
            if order != BirdAnalyzer.NAN and feature != BirdAnalyzer.NAN:
                order_to_feature[order] = feature
            if family != BirdAnalyzer.NAN and feature != BirdAnalyzer.NAN:
                family_to_feature[family] = feature
            if genera != BirdAnalyzer.NAN and feature != BirdAnalyzer.NAN:
                if genera.find(',') > 0:
                    for g in genera.split(','):
                        genus_to_feature[g.strip()] = feature
                else:
                    genus_to_feature[genera] = feature

        if BirdAnalyzer.DEBUG:
            print(genus_to_feature)
            print(family_to_feature)
            print(order_to_feature)

        taxonomy_to_feature = {}

        for index, row in taxonomy_master_table_df.iterrows():
            sn = str(row[BirdAnalyzer.TAXONOMY_FILE_SCIENTIFIC_NAME]).lower().strip()
            cn = str(row[BirdAnalyzer.TAXONOMY_FILE_COMMON_NAME]).lower().strip()
            order = str(row[BirdAnalyzer.TAXONOMY_FILE_ORDER]).lower().strip()
            family = str(row[BirdAnalyzer.TAXONOMY_FILE_FAMILY]).lower().strip()

            if family.find('(') > 0:
                family = family[0: family.index('(')].strip()

            if sn != BirdAnalyzer.NAN:
                genus = sn[0: sn.find(' ')].strip()
            else:
                genus = BirdAnalyzer.NAN

            if sn != BirdAnalyzer.NAN:
                if genus != BirdAnalyzer.NAN and genus in genus_to_feature:
                    taxonomy_to_feature[sn] = genus_to_feature[genus]
                elif family != BirdAnalyzer.NAN and family in family_to_feature:
                    taxonomy_to_feature[sn] = family_to_feature[family]
                elif order != BirdAnalyzer.NAN and order in order_to_feature:
                    taxonomy_to_feature[sn] = order_to_feature[order]

            if cn != BirdAnalyzer.NAN:
                if genus != BirdAnalyzer.NAN and genus in genus_to_feature:
                    taxonomy_to_feature[cn] = genus_to_feature[genus]
                elif family != BirdAnalyzer.NAN and family in family_to_feature:
                    taxonomy_to_feature[cn] = family_to_feature[family]
                elif order != BirdAnalyzer.NAN and order in order_to_feature:
                    taxonomy_to_feature[cn] = order_to_feature[order]

        return taxonomy_to_feature

    @staticmethod
    def run(taxonomy_master_file_name, taxonomy_to_feature_map_file_name, annotation_file, outfile):

        if os.path.isfile(outfile):
            print(outfile + ' already exists')
            sys.exit()

        o = open(outfile, 'w')

        taxonomy_to_feature = BirdAnalyzer.make_taxonomy_to_feature_map(taxonomy_master_file_name,
                                                                        taxonomy_to_feature_map_file_name)
        if BirdAnalyzer.DEBUG:
            print(taxonomy_to_feature)

        annotation_df = pd.read_csv(annotation_file, encoding='unicode_escape', dtype={
            BirdAnalyzer.ANNOTATION_FILE_HOST_NAME: str,
            BirdAnalyzer.ANNOTATION_FILE_COMMON_HOST_NAME: str,
            BirdAnalyzer.ANNOTATION_FILE_GB_ACC: str,
            BirdAnalyzer.ANNOTATION_FILE_GENOME_ID: str,
        }, low_memory=False)

        mapped_on_scientific_name = 0
        mapped_on_common_name = 0
        human = 0
        swine = 0
        other_mammals = 0
        total = 0
        not_mapped = 0

        for index, row in annotation_df.iterrows():
            total += 1
            host_name = str(row[BirdAnalyzer.ANNOTATION_FILE_HOST_NAME]).lower().strip()
            host_common_name = str(row[BirdAnalyzer.ANNOTATION_FILE_COMMON_HOST_NAME]).lower().strip()
            gb_acc = str(row[BirdAnalyzer.ANNOTATION_FILE_GB_ACC]).strip()
            genome_id = str(row[BirdAnalyzer.ANNOTATION_FILE_GENOME_ID]).strip()

            if BirdAnalyzer.DEBUG:
                print(host_name + ', ' + host_common_name)

            if host_common_name == 'human' or host_name == 'homo sapiens':
                human += 1
                o.write(gb_acc + '\t' + genome_id + '\t' + 'human' + '\n')
            elif host_common_name in BirdAnalyzer.SWINE or host_name in BirdAnalyzer.SWINE:
                swine += 1
                o.write(gb_acc + '\t' + genome_id + '\t' + 'swine' + '\n')
            elif host_common_name in BirdAnalyzer.OTHER_MAMMALS or host_name in BirdAnalyzer.OTHER_MAMMALS:
                other_mammals += 1
                o.write(gb_acc + '\t' + genome_id + '\t' + 'other_mammals' + '\n')
            else:
                if host_name in taxonomy_to_feature:
                    mapped_on_scientific_name += 1
                    if BirdAnalyzer.DEBUG_2:
                        print(host_name + " -> " + taxonomy_to_feature[host_name])
                    o.write(gb_acc + '\t' + genome_id + '\t' + taxonomy_to_feature[host_name] + '\n')
                elif host_common_name in taxonomy_to_feature:
                    mapped_on_common_name += 1
                    if BirdAnalyzer.DEBUG_2:
                        print(host_common_name + " -> " + taxonomy_to_feature[host_common_name])
                    o.write(gb_acc + '\t' + genome_id + '\t' + taxonomy_to_feature[host_common_name] + '\n')
                else:
                    # print(host_name + ', ' + host_common_name)
                    not_mapped += 1

        o.close()
        print()
        print("Total                          : " + str(total))
        print("Human                          : " + str(human))
        print("Swine                          : " + str(swine))
        print("Other mammals                  : " + str(other_mammals))
        print("Birds mapped on scientific_name: " + str(mapped_on_scientific_name))
        print("Birds mapped on common name    : " + str(mapped_on_common_name))
        print("Not mapped                     : " + str(not_mapped))


if __name__ == "__main__":
    # Uncomment this for testing in PyCharm:
    BirdAnalyzer.run('/Users/czmasek/Dropbox/WORK/JCVI/DL/ebird_oct_2022.csv',
                     '/Users/czmasek/Dropbox/WORK/JCVI/DL/ebird_to_feature_map.csv',
                     '/Users/czmasek/Dropbox/WORK/JCVI/DL/Flu_A_Complete.csv',
                     '/Users/czmasek/Dropbox/WORK/JCVI/DL/outfile.tsv')

# argument_parser = ap.ArgumentParser(prog='bird_analyzer',
#                                     description='mapping of taxonomy to feature')
#
# argument_parser.add_argument(dest='taxonomy_file', help='taxonomy master table (example \'ebird_oct_2022.csv\')',
#                              type=ap.FileType('r'))
# argument_parser.add_argument(dest='feature_map_file',
#                              help='taxonomy to feature map file (example \'ebird_to_feature_map.csv\')',
#                              type=ap.FileType('r'))
# argument_parser.add_argument(dest='annotation_file', help='annotation file (example \'Flu_A_Complete.csv\')',
#                              type=ap.FileType('r'))
# argument_parser.add_argument(dest='out_file', help='output file (example \'outfile.tsv\')',
#                              type=str)
# argument_parser.add_argument('--version', action='version', version='%(prog)s ' + BirdAnalyzer.VERSION)
#
# args = argument_parser.parse_args()
#
# tf = args.taxonomy_file
# fm = args.feature_map_file
# af = args.annotation_file
# of = args.out_file
#
# BirdAnalyzer.run(tf, fm, af, of)
