def gen_partitions_ms(n):
    """Generate all partitions of integer n (>= 0).
    by Tim Peter's code in:
    http://code.activestate.com/recipes/218332/


    Each partition is represented as a multiset, i.e. a dictionary
    mapping an integer to the number of copies of that integer in
    the partition.  For example, the partitions of 4 are {4: 1},
    {3: 1, 1: 1}, {2: 2}, {2: 1, 1: 2}, and {1: 4}.  In general,
    sum(k * v for k, v in a_partition.iteritems()) == n, and
    len(a_partition) is never larger than about sqrt(2*n).

    Note that the _same_ dictionary object is returned each time.
    This is for speed:  generating each partition goes quickly,
    taking constant time independent of n.
    """

    if n<0:
        raise ValueError("n must be >= 0")

    if n == 0:
        yield {}
        return

    ms = {n: 1}
    keys = [n]  # ms.keys(), from largest to smallest
    yield ms

    while keys != [1]:
        # Reuse any 1's.
        if keys[-1] == 1:
            del keys[-1]
            reuse = ms.pop(1)
        else:
            reuse = 0

        # Let i be the smallest key larger than 1.
        #  Reuse one instance of i.
        i = keys[-1]
        newcount = ms[i] = ms[i] - 1
        reuse += i
        if newcount == 0:
            del keys[-1], ms[i]

        # Break the remainder into pieces of size i-1.
        i -= 1
        q, r = divmod(reuse, i)
        ms[i] = q
        keys.append(i)
        if r:
            ms[r] = 1
            keys.append(r)
        yield ms

def partitions_rev(n):
    # reverse order
    if n == 0:
        yield []
        return
    for p in partitions_rev(n-1):
        yield p + [1]
        if p and (len(p) < 2 or p[-2] > p[-1]):
            yield p[:-1] + [p[-1] + 1]
            
# tested both performance by doing
#len(list(ruleAsc(60)))
#len(list(partitions_rev(60)))

def ruleAsc(n):
    """     Iterative Algorithm from Jerome Kelleher
    http://homepages.ed.ac.uk/jkellehe/partitions.php
    Algorithm to generate all ascending compositions.
    This algorithm is written as a Python generator.
    
    Although this algorithm is very simple, it is also very efficient.
    It is Constant Amortised Time.    
    """
    a = [0 for i in range(n + 1)]
    k = 1
    a[0] = 0
    a[1] = n
    while k != 0:
        x = a[k - 1] + 1
        y = a[k] - 1
        k -= 1
        while x <= y:
            a[k] = x
            y -= x
            k += 1
        a[k] = x + y
        yield a[:k + 1]