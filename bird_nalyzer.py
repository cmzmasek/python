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
            'order': str,
            'family': str
        })

    @staticmethod
    def read_family_to_feature_map_file(file_name):
        source_df = pd.read_csv(file_name, encoding='unicode_escape', dtype={
            'ORDER': str,
            'FAMILY': str,
            'GENERA': str,
            'FEATURE': str
        })


if __name__ == "__main__":
    pass
