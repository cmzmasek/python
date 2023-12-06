import re

def clean_str(s):
    return re.sub(r'\s+', ' ', s.strip())
