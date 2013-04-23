def bsearch_edge(items, key, f=cmp):
    
    """Returns the first i for which f(items[i - 1]) < key and f(items[i])
    >= key. An alternate comparison function can be given as f.

    >>> bsearch_edge('abcdefg', 'c')
    2

    >>> bsearch_edge('cdefg', 'a')
    0
    
    >>> bsearch_edge('abcde', 'g')
    5

    >>> bsearch_edge('abfg', 'e')
    2

    """

    p = 0
    q = len(items)


    while True:

        r = p + ((q - p) // 2)

        # Find the smallest i where item i-1's end is less than my end
        # and item i's end is greater than or equal to than my end.

        if r == 0 or f(items[r - 1], key) < 0:

            #            key
            # [ ... prev this ...]
            if f(items[r], key) >= 0:
                return r

            #                   key
            # [ ... prev this ] 
            elif r == len(items) - 1:
                return len(items)

            #                  key
            # [ ... prev this ...  ]
            else:
                p = r

        # key
        #    last
        else:
            q = r

def bsearch_overlap(items, key, key_left, key_right):

    cmp_left  = lambda a, b: cmp(key_right(a), key_left(b))

    i = bsearch_edge(items, key, cmp_left)

    while i < len(items) and key_left(items[i]) <= key_right(key):
        yield items[i]
        i += 1



