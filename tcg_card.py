import util as util


def set_comparator(comparator):
    __comparator = comparator


class Card(object):
    COMPARATOR_NAME = 0
    COMPARATOR_NAME_SET = 1
    COMPARATOR_NAME_SET_FINISH = 2
    COMPARATOR_NAME_SET_FINISH_ART = 3
    COMPARATOR_NAME_SET_ART = 4
    COMPARATOR_SCRYFALL_ID = 5
    COMPARATOR_NAME_SETCODE_NUMBER = 6

    __comparator = COMPARATOR_NAME_SETCODE_NUMBER

    def __init__(self, game, name, card_type, card_subtype, game_set, set_code, card_number, lang, rarity, artist,
                 finish, art_type, uri, scryfall_id, release_date, illustration_id):
        self.__game = util.clean_str(str(game))
        self.__name = util.clean_str(str(name))
        self.__card_type = util.clean_str(str(card_type))
        self.__card_subtype = util.clean_str(str(card_subtype))
        self.__game_set = util.clean_str(str(game_set))
        self.__set_code = util.clean_str(str(set_code)).upper()
        self.__card_number = util.clean_str(str(card_number))
        self.__lang = util.clean_str(str(lang)).lower()
        self.__rarity = util.clean_str(str(rarity)).lower()
        self.__artist = util.clean_str(str(artist))
        self.__finish = util.clean_str(str(finish)).lower()
        self.__art_type = util.clean_str(str(art_type)).lower()
        self.__uri = util.clean_str(str(uri)).lower()
        self.__scryfall_id = util.clean_str(str(scryfall_id)).lower()
        self.__release_date = util.clean_str(str(release_date)).lower()
        self.__illustration_id = util.clean_str(str(illustration_id))

    def get_game(self):
        return self.__game

    def get_name(self):
        return self.__name

    def get_type(self):
        return self.__card_type

    def get_subtype(self):
        return self.__card_subtype

    def get_set(self):
        return self.__game_set

    def get_set_code(self):
        return self.__set_code

    def get_card_number(self):
        return self.__card_number

    def get_lang(self):
        return self.__lang

    def get_rarity(self):
        return self.__rarity

    def get_artist(self):
        return self.__artist

    def get_finish(self):
        return self.__finish

    def get_art_type(self):
        return self.__art_type

    def get_uri(self):
        return self.__uri

    def get_scryfall_id(self):
        return self.__scryfall_id

    def get_release_date(self):
        return self.__release_date

    def get_illustration_id(self):
        return self.__illustration_id

    def __str__(self):
        s = '{} ({}) ({})'.format(self.get_name(), self.get_set(), self.get_card_number())
        if self.get_finish():
            s += ' [' + self.get_finish() + ']'
        if self.get_art_type():
            s += ' [' + self.get_art_type() + ']'
        return s

    def __repr__(self):
        if Card.__comparator == Card.COMPARATOR_NAME:
            return '{}::{}::{}'.format(self.__class__.__name__, self.get_game().lower(),
                                       self.get_name().lower())
        elif Card.__comparator == Card.COMPARATOR_NAME_SET:
            return '{}::{}::{}::{}'.format(self.__class__.__name__, self.get_game().lower(),
                                           self.get_name().lower(),
                                           self.get_set().lower())
        elif Card.__comparator == Card.COMPARATOR_NAME_SET_ART:
            return '{}::{}::{}::{}::{}'.format(self.__class__.__name__, self.get_game().lower(),
                                               self.get_name().lower(),
                                               self.get_set().lower(), self.get_art_type().lower())
        elif Card.__comparator == Card.COMPARATOR_NAME_SET_FINISH:
            return '{}::{}::{}::{}::{}'.format(self.__class__.__name__, self.get_game().lower(),
                                               self.get_name().lower(),
                                               self.get_set().lower(), self.get_finish())
        elif Card.__comparator == Card.COMPARATOR_NAME_SET_FINISH_ART:
            return '{}::{}::{}::{}::{}::{}'.format(self.__class__.__name__, self.get_game().lower(),
                                                   self.get_name().lower(),
                                                   self.get_set().lower(), self.get_finish().lower(),
                                                   self.get_art_type().lower())
        elif Card.__comparator == Card.COMPARATOR_SCRYFALL_ID:
            return '{}::{}::{}'.format(self.__class__.__name__, self.get_game().lower(),
                                       self.get_scryfall_id().lower())
        elif Card.__comparator == Card.COMPARATOR_NAME_SETCODE_NUMBER:
            return '{}::{}::{}::{}::{}'.format(self.__class__.__name__, self.get_game().lower(),
                                               self.get_name().lower(), self.get_set_code().upper(),
                                               self.get_card_number().lower())

    def __eq__(self, other):

        if isinstance(other, Card):
            return self.__repr__() == other.__repr__()
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.__repr__())

    def __lt__(self, other):
        return self.get_name() + '::' + self.get_set() < other.get_name() + '::' + other.get_set()


if __name__ == "__main__":
    card_0 = Card('Magic', 'Elspeth, Knight-Errant', 'Legendary Planeswalker', 'Elspeth',
                  'Mythic Edition: Guilds of Ravnica', 'MED', 'GR1', 'en', 'm', 'Zack Stella', 'foil', 'borderless', '',
                  '', '2022-06-13', '')

    card_1 = Card('Magic ', ' Elspeth,  Knight-Errant', 'Legendary  Planeswalker ', 'Elspeth ',
                  ' Mythic Edition: Guilds of Ravnica ', 'MED', 'GR1', 'En', 'm', 'Zack     Stella', 'foil',
                  'borderless', '',
                  '', '2022-06-13', '')

    print(repr(card_0))
    print(repr(card_1))

    print(card_0)
    print(card_1)
