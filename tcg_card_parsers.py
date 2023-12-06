import csv
import tcg_card as card
import util as util


class CardParsers(object):

    @staticmethod
    def parse_scryfall_spreadsheet(csv_infile):
        TYPE_SUPTYPE_SEP = ' — '
        GAME = 'Magic'
        cards = []
        with open(csv_infile, newline='\n') as infile:
            reader = csv.reader(infile, delimiter=',')
            for r in reader:
                if r[0] == 'section':
                    continue

                name = r[2]
                card_type_subtype = util.clean_str(r[4]).strip()
                card_subtype = ''
                game_set = r[5]
                set_code = r[6]
                card_number = r[7]
                lang = r[8]
                rarity = r[9]
                artist = r[10]
                finish = r[11].strip()
                art_type = ''
                uri = r[15]
                scryfall_id = r[16]
                release_date = ''
                illustration_id = ''

                if finish.lower() == 'nonfoil':
                    finish = ''
                if TYPE_SUPTYPE_SEP in card_type_subtype:
                    x = card_type_subtype.split(TYPE_SUPTYPE_SEP)
                    card_type = x[0]
                    card_subtype = x[1]
                else:
                    card_type = card_type_subtype

                c = card.Card(GAME, name, card_type, card_subtype, game_set, set_code, card_number, lang, rarity,
                              artist,
                              finish, art_type, uri, scryfall_id, release_date, illustration_id)
                cards.append(c)

            return cards

    def parse_manabox_export(csv_infile):
        GAME = 'Magic'
        cards = []
        with open(csv_infile, newline='\n') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for r in reader:
                if r[0] == 'BinderName':
                    continue

                name = r[2]
                card_type = ''
                card_subtype = ''
                set_code = r[3]
                set_name = r[4]
                card_number = r[5]
                foil = r[6].strip()
                rarity = r[7]
                scryfall_id = r[10]
                lang = r[15]
                artist = ''
                finish = ''
                art_type = ''
                uri = ''
                card_id = ''
                release_date = ''
                illustration_id = ''

                if foil.lower() == 'normal':
                    finish = ''


                c = card.Card(GAME, name, card_type, card_subtype, set_name, set_code, card_number, lang, rarity,
                              artist,
                              finish, art_type, uri, scryfall_id, release_date, illustration_id)
                cards.append(c)

            return cards
    @staticmethod
    def parse_tcg_player_collection(tcg_infile):
        GAME = 'Magic'
        cards = []
        with open(tcg_infile, newline='\n') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for r in reader:
                name_art_type_finish = util.clean_str(r[1]).strip()
                card_subtype = ''
                game_set = r[3].strip()
                card_type = ''
                set_code = ''
                card_number = ''
                lang = ''
                rarity = ''
                artist = ''
                finish = ''
                art_type = ''
                uri = ''
                scryfall_id = ''
                release_date = ''
                illustration_id = ''

                name_art_type_finish_lower = name_art_type_finish.lower()
                if '(halo foil)' in name_art_type_finish_lower:
                    art_type = 'Halo Foil'
                if '(gilded foil)' in name_art_type_finish_lower:
                    art_type = 'Gilded Foil'
                if '(foil etched)' in name_art_type_finish_lower:
                    art_type = 'Foil Etched'
                if '(extended art)' in name_art_type_finish_lower:
                    art_type = 'Extended Art'
                if '(borderless)' in name_art_type_finish_lower:
                    art_type = 'Borderless'
                if '(jp alternate art)' in name_art_type_finish_lower:
                    art_type = 'JP Alternate Art'
                if '[foil]' in name_art_type_finish_lower:
                    finish = 'foil'

                name = name_art_type_finish
                if finish:
                    name = name[0: name.find('- [')].strip()
                if art_type:
                    name = name[0: name.find('(')].strip()

                c = card.Card(GAME, name, card_type, card_subtype, game_set, set_code, card_number, lang, rarity,
                              artist,
                              finish, art_type, uri, scryfall_id, release_date, illustration_id)
                cards.append(c)

            return cards



