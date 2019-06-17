
def create_keys(ids):
    d = len(str(len(ids)))
    new_keys = [f"%0{d}d" % x for x in  list(range(len(ids)))]
#    return dict(zip(ids, new_keys))
    return [new_keys, ids]