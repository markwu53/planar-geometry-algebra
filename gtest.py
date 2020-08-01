
def gcd(seq):
    length = len(seq)
    aseq = [abs(e) for e in seq]
    amat = [[0]*i+[1]+[0]*(length-i-1) for i in range(length)]
    seq2 = [(v,i) for i,v in enumerate(aseq)]
    nz = [(v,i) for v,i in seq2 if v != 0]
    while len(nz) != 1:
        mv,mi = min(nz)
        seq3 = [(v%mv, v//mv, i) if i != mi else (v,0,i) for v,i in seq2]
        amat = [[v-n*u for v,u in zip(amat[i], amat[mi])] for v,n,i in seq3]
        seq2 = [(v,i) for v,n,i in seq3]
        nz = [(v,i) for v,i in seq2 if v != 0]
    mat = [[-a if v < 0 else a for a,v in zip(ax,seq)] for ax in amat]
    vert_space = [e for i,e in enumerate(mat) if i != mi]
    return mv, mat[mi], vert_space

if __name__ == "__main__":
    seq = [4, 26, -20, 94]
    result = gcd(seq)
    print(result)
