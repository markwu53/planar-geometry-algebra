import sympy


def simplify(x):
    return sympy.cancel(x)


class Point:

    def __init__(self, x, y):
        self.x = simplify(x)
        self.y = simplify(y)

    def subs(self, subs):
        return Point(self.x.subs(subs), self.y.subs(subs))

    def __str__(self):
        return str([self.x, self.y])

    def __add__(self, q):
        return Point(self.x + q.x, self.y + q.y)

    def __sub__(self, q):
        return Point(self.x - q.x, self.y - q.y)

    def __rmul__(self, s):
        return Point(s * self.x, s * self.y)

    def __iter__(self):
        self.pos = 0
        return self

    def __next__(self):
        self.pos += 1
        if self.pos == 1:
            return self.x
        if self.pos == 2:
            return self.y
        raise StopIteration

# all rational points can either be resolved by triangulation or line intersection.

# triangulation.
def findC(pA, pB, tA, tB):
    x = tB / (tA + tB)
    y = x * tA
    p1 = pB - pA
    p2 = Point(-p1.y, p1.x)
    return pA + x * p1 + y * p2

# line intersection
def x2line(l1, l2):
    m1, b1 = l1
    m2, b2 = l2
    x = (b1 - b2) / (m2 - m1)
    y = m1 * x + b1
    return Point(x, y)

# some line utility functions
def slope(pA, pB): return (pB.y - pA.y) / (pB.x - pA.x)
def linemp(m, p): return Point(m, p.y - m * p.x)
def line2p(pA, pB): return linemp(slope(pA, pB), pA)
def dropfoot(p, l): return x2line(l, linemp(-1 / l.x, p))

# some utility funcitons
def dot(v1, v2): return v1.x * v2.x + v1.y * v2.y
def dsqr(v): return dot(v, v)
def tanAplusB(tA, tB): return (tA + tB) / (1 - tA * tB)

# given a triangle, the in-circle center is not a rational point

# orthocenter resolved by line intersection
# it can also be resolved by triangulation
def orthocenter(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1 / m1, pB)
    l2 = linemp(-1 / m2, pA)
    return x2line(l1, l2)

# circumcenter resolved by line intersection
# it can also be resolved by triangulation
def circumcenter(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1 / m1, sympy.Rational(1, 2) * (pA + pC))
    l2 = linemp(-1 / m2, sympy.Rational(1, 2) * (pB + pC))
    return x2line(l1, l2)



#################

def run():
    #Feuerbach theorem
    x, y, z = sympy.symbols("x y z")
    pA = Point(0,0)
    pB = Point(1,0)
    tA = tanAplusB(x, x)
    tB = tanAplusB(y, y)
    pC = findC(pA, pB, tA, tB)
    pI = findC(pA, pB, x, y)
    pD = sympy.Rational(1,2) * (pA+pB)
    pE = sympy.Rational(1,2) * (pC+pB)
    pF = sympy.Rational(1,2) * (pA+pC)
    pO2 = circumcenter(pE, pF, pD)

    def c2tangent(pP,pA, pQ,pB):
        c = dsqr(pP-pQ)
        a = dsqr(pP-pA)
        b = dsqr(pQ-pB)
        return simplify((c-a-b)**2-4*a*b)

    # 3 ex-circle tangent to Euler circle
    pABo = findC(pB, pA, 1/y, 1/x)
    pAC = dropfoot(pABo, line2p(pA, pC))
    print(c2tangent(pABo, pAC, pO2, pD))

    pBCo = findC(pC, pB, tanAplusB(x, y), 1/y)
    pBC = dropfoot(pBCo, line2p(pB, pC))
    print(c2tangent(pBCo, pBC, pO2, pD))

    pACo = findC(pA, pC, 1/x, tanAplusB(x, y))
    pAC = dropfoot(pACo, line2p(pA, pC))
    print(c2tangent(pACo, pAC, pO2, pD))

    # in-circle tangent to Euler circle
    pIBC = dropfoot(pI, line2p(pC, pB))
    print(c2tangent(pO2, pD, pI, pIBC))

def run():
    #twomiles problem 1/2
    x, y, z = sympy.symbols("x y z")
    pB = Point(0,0)
    pC = Point(1,0)
    tB = tanAplusB(x, x)
    tC = tanAplusB(y, y)
    pA = findC(pB, pC, tB, tC)
    pI = findC(pB, pC, x, y)
    pD = x2line(line2p(pB, pI), line2p(pA, pC))
    pE = x2line(line2p(pC, pI), line2p(pA, pB))
    pO = circumcenter(pB, pC, pA)
    pBh = sympy.Rational(1,2) * (pB+pI)
    pCh = sympy.Rational(1,2) * (pC+pI)
    pG2 = findC(pBh, pI, -x, z)
    pG2o = pG2 + 2*(pO-pG2)
    pF2 = dropfoot(pG2o, line2p(pG2, pI))
    pG = Point(pG2o.x, pBh.y)
    pGo = pG + 2*(pO-pG)
    pF = dropfoot(pGo, line2p(pG, pI))
    a1,b1 = sympy.fraction(simplify(slope(pD, pE)))
    a2,b2 = sympy.fraction(simplify(slope(pD, pF)))
    a3,b3 = sympy.fraction(simplify(slope(pD, pF2)))
    t1 = simplify(a1*b2-a2*b1)
    t2 = simplify(a1*b3-a3*b1)
    constraint,_ = sympy.fraction(simplify(dsqr(pO-pG2)-dsqr(pO-pB)))
    print(sympy.factor(constraint))
    print(sympy.factor(t1))
    print(sympy.factor(t2))

def run():
    #twomiles problem t
    x, y, z, t = sympy.symbols("x y z t")
    pB = Point(0,0)
    pC = Point(1,0)
    tB = tanAplusB(x, x)
    tC = tanAplusB(y, y)
    pA = findC(pB, pC, tB, tC)
    pI = findC(pB, pC, x, y)
    pD = x2line(line2p(pB, pI), line2p(pA, pC))
    pE = x2line(line2p(pC, pI), line2p(pA, pB))
    pO = circumcenter(pB, pC, pA)
    pBh = pI+t*(pB-pI)
    pCh = pI+t*(pC-pI)
    pG2 = findC(pBh, pI, -x, z)
    constraint,_ = sympy.fraction(simplify(dsqr(pO-pG2)-dsqr(pO-pB)))
    print(sympy.factor(constraint))
    pG2o = pG2 + 2*(pO-pG2)
    pF2 = dropfoot(pG2o, line2p(pG2, pI))
    pG = Point(pG2o.x, pBh.y)
    pGo = pG + 2*(pO-pG)
    pF = dropfoot(pGo, line2p(pG, pI))
    m,b = line2p(pD, pE)
    pX = Point((pI.y-b)/m, pI.y)
    print(pX)
    mt,bt = line2p(pF, pF2)
    pXt = Point((pI.y-bt)/mt, pI.y)
    print(pXt)
    target,_ = sympy.fraction(simplify(pX.x-pXt.x))
    print(sympy.factor(target))

if __name__ == "__main__":
    run()
