import pandas as pd
import re


# Version 0.0.1
# Last modified 2023/11/06
# Christian M. Zmasek


@staticmethod
    def read_ebird_file(file_name):
        source_df = pd.read_csv(file_name, encoding='unicode_escape', dtype={
            "scientific name": str,
            "family": str
        })
    def read_family_to_feature_map_file():






