import pandas as pd


# Version 0.0.6
# Last modified 2023/15/06
# Roshni Bhattacharya
# Christian M. Zmasek

class BirdAnalyzer(object):
    DEBUG = False

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
    def make_taxonomy_to_feature_map(ebird_master_table_file_name, ebird_to_feature_map_file_name):
        taxonomy_master_table_df = BirdAnalyzer.read_taxonomy_file(ebird_master_table_file_name)
        taxonomy_to_feature_map_df = BirdAnalyzer.read_taxonomy_to_feature_map_file(ebird_to_feature_map_file_name)

        order_to_feature = {}
        family_to_feature = {}
        genus_to_feature = {}

        for index, row in taxonomy_to_feature_map_df.iterrows():

            order = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_ORDER])
            family = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_FAMILY])
            genera = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_GENERA])
            feature = str(row[BirdAnalyzer.TAXONOMY_MAP_FILE_FEATURE])
            if order != 'nan' and feature != 'nan':
                order_to_feature[order] = feature
            if family != 'nan' and feature != 'nan':
                family_to_feature[family] = feature
            if genera != 'nan' and feature != 'nan':
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
            sn = str(row[BirdAnalyzer.TAXONOMY_FILE_SCIENTIFIC_NAME])
            cn = str(row[BirdAnalyzer.TAXONOMY_FILE_COMMON_NAME])
            order = str(row[BirdAnalyzer.TAXONOMY_FILE_ORDER])
            family = str(row[BirdAnalyzer.TAXONOMY_FILE_FAMILY])

            if family.find('(') > 0:
                family = family[0: family.index('(')].strip()

            if sn != 'nan':
                genus = sn[0: sn.find(' ')].strip()
            else:
                genus = 'nan'

            if sn != 'nan':
                if genus != 'nan' and genus in genus_to_feature:
                    taxonomy_to_feature[sn] = genus_to_feature[genus]
                elif family != 'nan' and family in family_to_feature:
                    taxonomy_to_feature[sn] = family_to_feature[family]
                elif order != 'nan' and order in order_to_feature:
                    taxonomy_to_feature[sn] = order_to_feature[order]

            if cn != 'nan':
                if genus != 'nan' and genus in genus_to_feature:
                    taxonomy_to_feature[cn] = genus_to_feature[genus]
                elif family != 'nan' and family in family_to_feature:
                    taxonomy_to_feature[cn] = family_to_feature[family]
                elif order != 'nan' and order in order_to_feature:
                    taxonomy_to_feature[cn] = order_to_feature[order]

        return taxonomy_to_feature

    @staticmethod
    def run(taxonomy_master_file_name, taxonomy_to_feature_map_file_name, annotation_file):

        taxonomy_to_feature = BirdAnalyzer.make_taxonomy_to_feature_map(taxonomy_master_file_name,
                                                                        taxonomy_to_feature_map_file_name)
        if BirdAnalyzer.DEBUG:
            print(taxonomy_to_feature)

        annotation_df = pd.read_csv(annotation_file, encoding='unicode_escape', dtype={
            BirdAnalyzer.ANNOTATION_FILE_HOST_NAME: str,
            BirdAnalyzer.ANNOTATION_FILE_COMMON_HOST_NAME: str
        }, low_memory=False)

        mapped_on_scientific_name = 0
        mapped_on_common_name = 0
        human = 0
        swine = 0
        total = 0
        not_mapped = 0

        for index, row in annotation_df.iterrows():
            total += 1
            host_name = str(row[BirdAnalyzer.ANNOTATION_FILE_HOST_NAME])
            host_common_name = str(row[BirdAnalyzer.ANNOTATION_FILE_COMMON_HOST_NAME])
            if BirdAnalyzer.DEBUG:
                print(host_name + ', ' + host_common_name)

            if host_common_name.lower() == 'human' or host_name.lower() == 'homo sapiens':
                human += 1
            elif host_common_name.lower() == 'pig' or host_common_name.lower() == 'swine' \
                    or host_name.lower() == 'swine' or host_name.lower() == 'sus scrofa domesticus' \
                    or host_name.lower() == 'sus scrofa':
                swine += 1
            else:
                if host_name in taxonomy_to_feature:
                    mapped_on_scientific_name += 1
                    print(host_name + " -> " + taxonomy_to_feature[host_name])
                elif host_common_name in taxonomy_to_feature:
                    mapped_on_common_name += 1
                    print(host_common_name + " -> " + taxonomy_to_feature[host_common_name])
                else:
                    not_mapped += 1

        print()
        print("Total                          : " + str(total))
        print("Human                          : " + str(human))
        print("Swine                          : " + str(swine))
        print("Birds mapped on scientific_name: " + str(mapped_on_scientific_name))
        print("Birds mapped on common name    : " + str(mapped_on_common_name))
        print("Not mapped                     : " + str(not_mapped))


if __name__ == "__main__":
    BirdAnalyzer.run('/Users/czmasek/Dropbox/WORK/JCVI/DL/ebird_oct_2022.csv',
                     '/Users/czmasek/Dropbox/WORK/JCVI/DL/ebird_to_feature_map.csv',
                     '/Users/czmasek/Dropbox/WORK/JCVI/DL/Flu_A_Complete.csv')
