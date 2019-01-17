def isiterable(obj):
    ''' apply iter on obj to see if it works and obj is iterable or not'''

    iterable = True
    try:
        iter(obj)
    except TypeError:
        iterable = False

    return iterable
