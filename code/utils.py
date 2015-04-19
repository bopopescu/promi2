def lmean (inlist):
    """returns the aritmetic mean of the value in the passed list"""
    sum = 0
    for item in inlist:
        sum = sum + item
    return sum/float(len(inlist))

def random_string(l):
    import string
    import random
    return ''.join(random.choice(string.ascii_uppercase+string.digits) for i in range(l))

def line_count(f):
    num_lines = sum(1 for line in open(f))
    return num_lines

def all_same(items):
    return all(x == items[0] for x in items)

def get_value_from_keycolonvalue_list(string, kcvlist, sep=':'):
    ## get value from a keycolonvalue list: ['a:x', 'b:y', 'c:z']
    ## given a string before colon, returns the string after the colon
    try:
        return [i for i in kcvlist if i.startswith(string+sep)][0].split(sep)[1]
    except IndexError:
        return ''

def ensure_dir(f, is_file=True):
    import os
    if is_file:
        d = os.path.dirname(f)
    else:
        d = f
    if not os.path.exists(d):
        os.makedirs(d)
