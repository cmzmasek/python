import pandas as pd
import re


# Version 0.0.1
# Last modified 2023/11/06
# Christian M. Zmasek

class BirdAnalyzer(object):

    @staticmethod
    def read_ebird_file(file_name):
        source_df = pd.read_csv(file_name, encoding='unicode_escape', dtype={
            'scientific name': str,
            'English name': str,
            'order': str,
            'family': str
        })
        return source_df

    @staticmethod
    def read_family_to_feature_map_file(file_name):
        source_df = pd.read_csv(file_name, encoding='unicode_escape', dtype={
            'ORDER': str,
            'FAMILY': str,
            'GENERA': str,
            'FEATURE': str
        })
        return source_df

    @staticmethod
    def read(ebird_master_table_file_name, ebird_to_feature_map_file_name):
        ebird_master_table_df = BirdAnalyzer.read_ebird_file(ebird_master_table_file_name)
        ebird_to_feature_map_df = BirdAnalyzer.read_family_to_feature_map_file(ebird_to_feature_map_file_name)

        order_to_feature = {}
        family_to_feature = {}
        genus_to_feature = {}

        for index, row in ebird_to_feature_map_df.iterrows():
            order = str(row['ORDER'])
            family = str(row['FAMILY'])
            genera = str(row['GENERA'])
            feature = str(row['FEATURE'])
            if order != 'nan' and feature != 'nan':
                order_to_feature[order] = feature
            if family != 'nan' and feature != 'nan':
                family_to_feature[family] = feature
            if genera != 'nan' and feature != 'nan':
                if genera.find(',') > 0:
                    for g in genera.split(','):
                        genus_to_feature[g] = feature
                else:
                    genus_to_feature[genera] = feature

        species_to_feature = {}

        for index, row in ebird_master_table_df.iterrows():
            sn = str(row['scientific name'])
            cn = str(row['English name'])
            order = str(row['order'])
            family = str(row['family'])

            if family.find('(') > 0:
                family = family[0: family.index('(')].strip()
            genus = sn[0: sn.find(' ')].strip()

            if sn != 'nan':
                if order != 'nan' and order in order_to_feature:
                    species_to_feature[sn] = order_to_feature[order]
                    if family != 'nan' and family in family_to_feature:
                        species_to_feature[sn] = family_to_feature[family]
                        if genus in genus_to_feature:
                            species_to_feature[sn] = genus_to_feature[genus]

        print(species_to_feature)


if __name__ == "__main__":
    BirdAnalyzer.read('/Users/czmasek/Dropbox/WORK/JCVI/DL/ebird_oct_2022.csv',
                      '/Users/czmasek/Dropbox/WORK/JCVI/DL/ebird_to_feature_map.csv')
