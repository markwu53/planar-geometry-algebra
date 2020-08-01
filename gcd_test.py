
def gcd1(vec):
    coll = [abs(e) for e in vec]
    #assume all nonnegative integers
    length = len(coll)
    am = [[0]*i+[1]+[0]*(length-i-1) for i in range(length)]
    coll2 = [(v,i) for i,v in enumerate(coll) if v != 0]
    while len(coll2) != 1:
        print(coll2)
        mv,mi = min(coll2)
        coll3 = [(v%mv,v//mv,i) if i != mi else (v,0,i) for (v,i) in coll2]
        for v,n,i in coll3:
            am[i] = [v-n*u for v,u in zip(am[i], am[mi])]
        coll2 = [(v,i) for v,n,i in coll3 if v != 0]
        print(coll3)
    print(coll2)
    for i,v in enumerate(vec):
        if v < 0:
            for e in am:
                e[i] = -1 * e[i]
    for e in am: print(e)

def gcd(seq):
    aseq = [abs(e) for e in seq]
    amat = [[0]*i+[1]+[0]*(len(seq)-i-1) for i in range(len(seq))]
    seq2 = [(v,i) for i,v in enumerate(aseq) if v != 0]
    while len(seq2) != 1:
        mv,mi = min(seq2)
        seq3 = [(v%mv, v//mv, i) if i != mi else (v,0,i) for v,i in seq2]
        amat = [[v-n*u for v,u in zip(amat[i], amat[mi])] for v,n,i in seq3]
        seq2 = [(v,i) for v,n,i in seq3 if v != 0]
    mat = [[-a if v < 0 else a for a,v in zip(ax,seq)] for ax in amat]
    return seq2, mat

def run():
    coll = [-22, 4, 26, 94]
    result = gcd(coll)
    print(result)

if __name__ == "__main__":
    run()
