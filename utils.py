""" Utilities for quarantine repo:
    - Generic helpers for iterables
    - Helpers for plotting things
"""

def selector(els, idxs):
    for idx in idxs:
        yield els[idx]

def tuple_filter(tup_iter, idx):
    return [tup[idx] for tup in tup_iter]

def invert_dict(d):
    new_dict = {}
    for k,v in d.items():
        if v in new_dict:
            new_dict[v].append(k)
        else:
            new_dict[v] = [k]
    return new_dict

def argmax(iterable):
    # Returns max_idx, max_value
    max_val, max_idx = -float('inf'), None
    for i, el in enumerate(iterable):
        if el > max_val:
            max_val = el
            max_idx = i
    return max_idx, max_val

def mean(iterable, lambda_=None):
    count, runsum = 0, 0.0 
    if lambda_ is None:
        lambda_ = lambda x: x 
    for el in iterable:
        count += 1
        runsum += lambda_(el)
    return runsum / float(count)
    

def mergesum(dicts):
    """ Given a list/iterable of (nested) dicts, will merge them together 
        where merge at the base level means summing values for shared keys 
    """
    def looper(iter_input, output=None):
        if all(isinstance(_, (float, int)) for _ in iter_input):
            return sum(iter_input)
        # Collect by shared keys:
        shared_keys = {}
        for el in iter_input:
            for k,v in el.items():
                if k not in shared_keys:
                    shared_keys[k] = []
                shared_keys[k].append(v)

        return {k: looper(v) for k, v in shared_keys.items()}
    return looper(dicts)



### Plotting helpers

def c(i):
    return 'bgrcmyk'[i]

