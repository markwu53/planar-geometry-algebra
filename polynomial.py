
def monic(d):
    return d+[1]

def shift(p, pos):
    return  [0]*pos + p

def mult(p, n):
    return [i*n for i in p]

def poly(p):
    pos = 1
    while True:
        if pos > len(p): return [0]
        if p[-pos] != 0: break
        pos += 1
    if pos == 1: return p
    return p[:1-pos]

def minus(p, q):
    q0 = q+[0]*(len(p)-len(q))
    return poly([a-b for a,b in zip(p, q0)])

def mr(q, d):
    m = monic(d)
    def remainder(p):
        if len(p) < len(m):
            return p
        base = mult(shift(m, len(p)-len(m)), p[-1])
        return remainder(minus(p, base))
    return remainder(q)

def mult2(p, q):
    result = {}
    for ((i,j),v) in p.items():
        for ((a,b),u) in q.items():
            if (i+a,j+b) not in result:
                result[(i+a,j+b)] = []
            result[(i+a,j+b)].append(v*u)
    result = {k:sum(result[k]) for k in result}
    return result

def x_plus_y_n(n):
    p={}
    p[(0,1)] = 1
    p[(1,0)] = 1
    r = p
    for i in range(n-1):
        r = mult2(r, p)
    return r

def x2x(p):
    result = {}
    for i,v in enumerate(p):
        if v == 0: continue
        result[(i,0)] = v
    return result

def x2y(p):
    result = {}
    for i,v in enumerate(p):
        if v == 0: continue
        result[(0,i)] = v
    return result

def add2(p, q):
    result = {k:p[k] for k in p}
    for k in q:
        if k in result:
            result[k] += q[k]
        else:
            result[k] = q[k]
    return result

def q2(p, mx, my):
    result = {}
    for i,j in p:
        v = p[(i,j)]
        px = mr([0]*i+[v], mx)
        py = mr([0]*j+[1], my)
        result = add2(result, mult2(x2x(px), x2y(py)))
    return result

def pd1(p2, X, Y):
    result = []
    for i in range(X):
        for j in range(Y):
            if (i,j) in p2:
                result.append(p2[(i,j)])
            else:
                result.append(0)
    return result

def x_plus_y_base():
    mx = [4,4]
    my = [4,4,0]
    for n in range(1, 20):
        r = x_plus_y_n(n)
        result = q2(r, mx, my)
        result = pd1(result, len(mx), len(my))
        print(result)

def x_times_y_base():
    mx = [4,4]
    my = [4,4,0]
    coll = [[1]+[0]*(len(mx)*len(my)-1)]
    for n in range(1, 10):
        r = {(n,n):1}
        result = q2(r, mx, my)
        result = pd1(result, len(mx), len(my))
        coll.append(result)
    return coll

def linear_dep(coll):
    buff = [coll]
    used = []
    length = len(coll[0])
    height = len(coll)
    for ind in range(length):
        col = [e[ind] for e in coll]
        col_list = [(abs(v),i,v) for i,v in enumerate(col)]
        m = min([e for e in col_list if e[0] != 0 and e[1] not in used])
        used.append(m[1])
        col2 = [e//m[2] for e in col]
        coll = [[x-y*col2[i] for x,y in zip(coll[i], coll[m[1]])] if i!= m[1] else coll[i] for i in range(height)]
        #coll = [[e-m[2]*k for e in v] if i!=m[1] else v for i,(v,k) in enumerate(zip(coll, col2))]
        buff.append(coll)
        print(m, col2)
    for coll in buff:
        print("one")
        for e in coll: print(e)

def gcd(coll):
    coll = [abs(e) for e in coll if e != 0]
    if len(coll) == 1:
        return coll[0]
    coll = sorted(coll)
    m = coll[0]
    coll = [m]+[e%m for e in coll[1:]]
    return gcd(coll)

def gcd2(coll):
    path = []
    vi = [(v,i) for i,v in enumerate(coll)]
    while True:
        mv,mi = min([e for e in vi if e[0] != 0])
        vni = [(v%mv,v//mv,i) if i != mi else (v,0,i) for v,i in vi]
        path.append([(mv,mi), vni])
        nz = [e for e in vni if e[0] != 0]
        if len(nz) == 1: break
        vi = [(v,i) for v,n,i in vni]
    for e in path: print(e)

def run():
    coll = [124, 123456, 236]
    gcd2(coll)
    a = [-2983, 3, -2]
    r = sum([x*y for x,y in zip(coll, a)])
    print(r)

if __name__ == "__main__":
    run()
