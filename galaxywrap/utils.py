def isiterable(obj):
    ''' apply iter on obj to see if it works and obj is iterable or not'''

    iterable = True
    try:
        iter(obj)
    except TypeError:
        iterable = False

    return iterable


def change_tuple_unit(tpl, unit):
    tpl = list(tpl)
    if unit is not None:
        for i, value in enumerate(tpl):
            if value is not None:
                tpl[i] = value.to(unit)

    return tpl


def translate_to_constraints_names(name):
    if name == 'r':
        name = 're'
    return name


def check_all_or_no_None(tpl):
    if None in tpl:
        for entry in tpl:
            assert entry is None
