import pandas as pd
import re


class HostAnalyzer(object):
    SOURCE_SCI_NAME = 'HOST_NAME'
    SOURCE_COMMON_NAME = 'HOST_COMMON_NAME'
    TARGET_SCI_NAME = 'SCI_NAME'
    TARGET_COMMON_NAME = 'PRIMARY_COM_NAME'
    TARGET_ORDER1 = 'ORDER1'
    TARGET_FAMILY = 'FAMILY'
    TARGET_HOST = 'HOST'
    TARGET_MENU1 = 'MENU1'
    TARGET_MENU2 = 'MENU2'
    TARGET_MENU3 = 'MENU3'

    PAREN_RE = re.compile("(.+?)\s*\((.+)\)")
    COMMA_RE = re.compile("(.+?)\s*,\s*(.+)")

    @staticmethod
    def write_match(file, source_sci_name, source_common_name, match_type, df):
        file.write(
            source_sci_name + '\t' + source_common_name + '\t' + match_type + '\t' + str(
                df.iloc[0][HostAnalyzer.TARGET_SCI_NAME]) + '\t' +
            str(df.iloc[0][HostAnalyzer.TARGET_COMMON_NAME]) + '\t' + str(
                df.iloc[0][HostAnalyzer.TARGET_ORDER1]) + '\t' +
            str(df.iloc[0][HostAnalyzer.TARGET_FAMILY]) + '\t' + str(df.iloc[0][HostAnalyzer.TARGET_HOST]) + '\t' +
            str(df.iloc[0][HostAnalyzer.TARGET_MENU1]) + '\t' + str(df.iloc[0][HostAnalyzer.TARGET_MENU2]) + '\t' +
            str(df.iloc[0][HostAnalyzer.TARGET_MENU3]) + '\n')

    @staticmethod
    def write_no_match(file, source_sci_name, source_common_name):
        file.write(
            source_sci_name + '\t' + source_common_name + '\n')

    @staticmethod
    def read(source_file_name, target_file_name):
        match_outfile = open("matches.txt", "w")
        no_match_outfile = open("no_matches.txt", "w")

        match_outfile.write(
            'SOURCE_SCI_NAME\tSOURCE_COMMON_NAME\tMATCH_TYPE\t' + HostAnalyzer.TARGET_SCI_NAME + '\t'
            + HostAnalyzer.TARGET_COMMON_NAME + '\t' + HostAnalyzer.TARGET_ORDER1 + '\t' + HostAnalyzer.TARGET_FAMILY + '\t' +
            HostAnalyzer.TARGET_HOST + '\t' + HostAnalyzer.TARGET_MENU1 + '\t' + HostAnalyzer.TARGET_MENU2
            + '\t' + HostAnalyzer.TARGET_MENU3 + '\n')

        no_match_outfile.write('SOURCE_SCI_NAME\tSOURCE_COMMON_NAME\n')

        source_df = pd.read_csv(source_file_name, encoding='unicode_escape', dtype={
            HostAnalyzer.SOURCE_SCI_NAME: str,
            HostAnalyzer.SOURCE_COMMON_NAME: str
        })

        target_df = pd.read_csv(target_file_name, sep='\t', encoding='unicode_escape', dtype={
            HostAnalyzer.TARGET_SCI_NAME: str,
            HostAnalyzer.TARGET_COMMON_NAME: str,
            HostAnalyzer.TARGET_ORDER1: str,
            HostAnalyzer.TARGET_FAMILY: str,
            HostAnalyzer.TARGET_HOST: str,
            HostAnalyzer.TARGET_MENU1: str,
            HostAnalyzer.TARGET_MENU2: str,
            HostAnalyzer.TARGET_MENU3: str
        })

        match_counter = 0
        no_match_counter = 0
        multiple_match_counter = 0
        match_sci_to_sci = 0
        match_common_to_sci = 0
        match_sci_to_common = 0
        match_common_to_common = 0
        #
        match_sci_wo_paren_to_sci = 0
        match_sci_wo_paren_to_common = 0
        match_sci_in_paren_to_sci = 0
        match_sci_in_paren_to_common = 0
        match_sci_before_comma_to_sci = 0
        match_sci_before_comma_to_common = 0
        match_sci_after_comma_to_sci = 0
        match_sci_after_comma_to_common = 0

        for index, row in source_df.iterrows():
            source_sci_name = str(row[HostAnalyzer.SOURCE_SCI_NAME]).strip()
            source_common_name = str(row[HostAnalyzer.SOURCE_COMMON_NAME]).strip()
            paren_found = re.search(HostAnalyzer.PAREN_RE, source_sci_name)
            comma_found = re.search(HostAnalyzer.COMMA_RE, source_sci_name)
            match = False

            if paren_found:
                source_sci_name_wo_paren = paren_found.group(1).strip()
                source_sci_name_in_paren = paren_found.group(2).strip()
                sci_wo_paren_to_sci_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_SCI_NAME] == source_sci_name_wo_paren]
                sci_in_paren_to_sci_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_SCI_NAME] == source_sci_name_in_paren]
                sci_wo_paren_to_common_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_COMMON_NAME] == source_sci_name_wo_paren]
                sci_in_paren_to_common_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_COMMON_NAME] == source_sci_name_in_paren]
                if not sci_wo_paren_to_sci_df.empty:
                    match = True
                    match_sci_wo_paren_to_sci += 1
                    if len(sci_wo_paren_to_sci_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_before_paren_to_sci',
                                                 sci_wo_paren_to_sci_df)
                if not sci_in_paren_to_sci_df.empty:
                    match = True
                    match_sci_in_paren_to_sci += 1
                    if len(sci_in_paren_to_sci_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_in_paren_to_sci',
                                                 sci_in_paren_to_sci_df)
                if not sci_wo_paren_to_common_df.empty:
                    match = True
                    match_sci_wo_paren_to_common += 1
                    if len(sci_wo_paren_to_common_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_before_paren_to_common',
                                                 sci_wo_paren_to_common_df)
                if not sci_in_paren_to_common_df.empty:
                    match = True
                    match_sci_in_paren_to_common += 1
                    if len(sci_in_paren_to_common_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_in_paren_to_common',
                                                 sci_in_paren_to_common_df)
            if comma_found:
                source_sci_name_before_comma = comma_found.group(1).strip()
                source_sci_name_after_comma = comma_found.group(2).strip()
                sci_before_comma_to_sci_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_SCI_NAME] == source_sci_name_before_comma]
                sci_after_comma_to_sci_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_SCI_NAME] == source_sci_name_after_comma]
                sci_before_comma_to_common_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_COMMON_NAME] == source_sci_name_before_comma]
                sci_after_comma_to_common_df = target_df.loc[
                    target_df[HostAnalyzer.TARGET_COMMON_NAME] == source_sci_name_after_comma]
                if not sci_before_comma_to_sci_df.empty:
                    match = True
                    match_sci_before_comma_to_sci += 1
                    if len(sci_before_comma_to_sci_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_before_comma_to_sci',
                                                 sci_before_comma_to_sci_df)
                if not sci_after_comma_to_sci_df.empty:
                    match = True
                    match_sci_after_comma_to_sci += 1
                    if len(sci_after_comma_to_sci_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_after_comma_to_sci',
                                                 sci_after_comma_to_sci_df)
                if not sci_before_comma_to_common_df.empty:
                    match = True
                    match_sci_before_comma_to_common += 1
                    if len(sci_before_comma_to_common_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_before_comma_to_common',
                                                 sci_before_comma_to_common_df)
                if not sci_after_comma_to_common_df.empty:
                    match = True
                    match_sci_after_comma_to_common += 1
                    if len(sci_after_comma_to_common_df.index) > 1:
                        multiple_match_counter += 1
                    else:
                        HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name,
                                                 'sci_after_comma_to_common',
                                                 sci_after_comma_to_common_df)
            sci_to_sci_df = target_df.loc[target_df[HostAnalyzer.TARGET_SCI_NAME] == source_sci_name]
            common_to_sci_df = target_df.loc[target_df[HostAnalyzer.TARGET_SCI_NAME] == source_common_name]
            sci_to_common_df = target_df.loc[target_df[HostAnalyzer.TARGET_COMMON_NAME] == source_sci_name]
            common_to_common_df = target_df.loc[target_df[HostAnalyzer.TARGET_COMMON_NAME] == source_common_name]

            if not sci_to_sci_df.empty:
                match = True
                match_sci_to_sci += 1
                if len(sci_to_sci_df.index) > 1:
                    multiple_match_counter += 1
                else:
                    HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name, 'sci_to_sci',
                                             sci_to_sci_df)
            if not common_to_sci_df.empty:
                match = True
                match_common_to_sci += 1
                if len(common_to_sci_df.index) > 1:
                    multiple_match_counter += 1
                else:
                    HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name, 'common_to_sci',
                                             common_to_sci_df)
            if not sci_to_common_df.empty:
                match = True
                match_sci_to_common += 1
                if len(sci_to_common_df.index) > 1:
                    multiple_match_counter += 1
                else:
                    HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name, 'sci_to_common',
                                             sci_to_common_df)
            if not common_to_common_df.empty:
                match = True
                match_common_to_common += 1
                if len(common_to_common_df.index) > 1:
                    multiple_match_counter += 1
                else:
                    HostAnalyzer.write_match(match_outfile, source_sci_name, source_common_name, 'common_to_common',
                                             common_to_common_df)
            if match:
                match_counter += 1
            else:
                no_match_counter += 1
                HostAnalyzer.write_no_match(no_match_outfile, source_sci_name, source_common_name)

        match_outfile.close()
        no_match_outfile.close()

        print('match                 : ' + str(match_counter))
        print('match sci to sci      : ' + str(match_sci_to_sci))
        print('match common to sci   : ' + str(match_common_to_sci))
        print('match sci to common   : ' + str(match_sci_to_common))
        print('match common to common: ' + str(match_common_to_common))
        print('match sci before paren to sci   : ' + str(match_sci_wo_paren_to_sci))
        print('match sci before paren to common: ' + str(match_sci_wo_paren_to_common))
        print('match sci in paren to sci       : ' + str(match_sci_in_paren_to_sci))
        print('match sci in paren to common    : ' + str(match_sci_in_paren_to_common))
        print('match sci before comma to sci   : ' + str(match_sci_before_comma_to_sci))
        print('match sci before comma to common: ' + str(match_sci_before_comma_to_common))
        print('match sci after comma to sci    : ' + str(match_sci_after_comma_to_sci))
        print('match sci after comma to common : ' + str(match_sci_after_comma_to_common))

        print('no match: ' + str(no_match_counter))
        print('multiple match counter: ' + str(multiple_match_counter))


if __name__ == "__main__":
    HostAnalyzer.read('C:/Users/czmasek/Dropbox/WORK/JCVI/HOSTS/BVBRC_host_lookup.csv',
                      'C:/Users/czmasek/Dropbox/WORK/JCVI/HOSTS/OrderFamSpec-MENUS1-3.txt'
                      )
