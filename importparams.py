import re
"""
1) Match all lines with "const" at beginning
2) Strip them.
3) Split them by "="
4) Dictionary assign them.
"""

def num(s):
    try:
        return int(s)
    except ValueError:
        return eval(s)

def parse_params():
    # define regular expression pattern
    pattern = "[\w]+[ ]+=[\w\W]+$"

    name = "[A-Z\d_]+[ ]+"
    eq = "="
    value = "[\w\W ]+;$"

    # import file
    filename = "armparams.h"
    with open(filename, "r") as f:
        raw_lines = f.readlines()

    # find lines with "const" at beginning
    pattern = "^const"
    consts = []
    for raw_line in raw_lines:
        found = re.search(pattern, raw_line)
        if found:
            consts.append(raw_line)

    # strip lines
    pattern = "[\w]+[\s]*=[\s]*[^;]+"
    stripped = []
    i = 0
    for raw_line in consts:
        found = re.search(pattern, raw_line)
        if found:
            stripped.append(found.group())

    # split lines
    names = []
    values = []
    params = {}
    for raw_line in stripped:
        new_line = re.split("=", raw_line) # split by =
        new_line[0] = re.sub("\s", "", new_line[0]) # strip whitespace
        new_line[1] = re.sub("[()]", "", new_line[1]) # remove parentheses
        params[new_line[0]] = num(new_line[1])
        #i += 1
    #print i
    #for p in params:
    #    print p, "=", params[p]
    
    return params