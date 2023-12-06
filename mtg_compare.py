import argparse as ap
import os
import sys

import tcg_card_parsers as parser


class MtgCompare(object):
    VERSION = '1.0.0'

    @staticmethod
    def compare(scryfall_infile, manabox_infile, outfile_base):
        scryfall_cards_list = parser.CardParsers.parse_scryfall_spreadsheet(scryfall_infile)
        manabox_cards_list = parser.CardParsers.parse_manabox_export(manabox_infile)

        print('Scryfall (all)                : ' + str(len(scryfall_cards_list)))
        print('Manabox (all)                 : ' + str(len(manabox_cards_list)))

        scryfall_cards = set(scryfall_cards_list)
        manabox_cards = set(manabox_cards_list)

        print('Scryfall (set)                : ' + str(len(scryfall_cards)))
        print('ManaBox (set)                 : ' + str(len(manabox_cards)))

        intersection = manabox_cards.intersection(scryfall_cards)
        print('Intersection                  : ' + str(len(intersection)))

        scryfall_minus_manabox = scryfall_cards.difference(manabox_cards)
        print('Scryfall minus ManaBox (SF-MB): ' + str(len(scryfall_minus_manabox)))

        manabox_minus_scryfall = manabox_cards.difference(scryfall_cards)
        print('ManaBox minus Scryfall (MB-SF): ' + str(len(manabox_minus_scryfall)))

        of = outfile_base + '__Intersection.tsv'
        if os.path.isfile(of):
            print('\n' + of + ' already exists' + '\n')
            sys.exit()

        with open(of, 'w') as f:
            for c in sorted(intersection):
                f.write(
                    c.get_name() + '\t(' + c.get_set() + ')' + '\t' + c.get_set_code() +
                    '\t' + c.get_card_number() + ' \t' + c.get_artist() +
                    '\t' + c.get_finish())
                f.write('\n')

        of = outfile_base + '__SF-MB.tsv'
        if os.path.isfile(of):
            print('\n' + of + ' already exists' + '\n')
            sys.exit()

        with open(of, 'w') as f:
            for c in sorted(scryfall_minus_manabox):
                f.write(
                    c.get_name() + '\t(' + c.get_set() + ')' + '\t' + c.get_set_code() +
                    '\t' + c.get_card_number() + ' \t' + c.get_artist() +
                    '\t' + c.get_finish())
                f.write('\n')

        of = outfile_base + '__MB-SF.tsv'
        if os.path.isfile(of):
            print('\n' + of + ' already exists' + '\n')
            sys.exit()

        with open(of, 'w') as f:
            for c in sorted(manabox_minus_scryfall):
                f.write(
                    c.get_name() + '\t(' + c.get_set() + ')' + '\t' + c.get_set_code() +
                    '\t' + c.get_card_number() + ' \t' + c.get_artist() +
                    '\t' + c.get_finish())
                f.write('\n')


if __name__ == '__main__':
    argument_parser = ap.ArgumentParser(prog='mtg_compare')

    argument_parser.add_argument(dest='scryfall_infile', help='scryfall export (example \'scryfall.csv\')',
                                 type=str)

    argument_parser.add_argument(dest='manabox_infile',
                                 help='manabox file (example \'ManaBox_Collection.tsv\')',
                                 type=str)
    argument_parser.add_argument(dest='outfile_base',
                                 help='outfile base (example \'A\')',
                                 type=str)

    argument_parser.add_argument('--version', action='version', version='%(prog)s ' + MtgCompare.VERSION)

    args = argument_parser.parse_args()

    scryfall = args.scryfall_infile
    manabox = args.manabox_infile
    outb = args.outfile_base

    MtgCompare.compare(scryfall, manabox, outb)
