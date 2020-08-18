import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import sympy
import sympy.solvers

def simplify(x):
    return sympy.cancel(x)

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def simplify(self):
        return Point(simplify(self.x),simplify(self.y))

    def __str__(self):
        return str([self.x, self.y])

    def __add__(self, q):
        return Point(self.x+q.x, self.y+q.y)

    def __sub__(self, q):
        return Point(self.x-q.x, self.y-q.y)

    def __rmul__(self, s):
        return Point(s*self.x, s*self.y)

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

    def evalf(self, subs):
        return np.array([self.x.evalf(subs=subs), self.y.evalf(subs=subs)])

fig, ax = plt.subplots()
#plt.xlim(-.5,2.)
#plt.ylim(-.6,1.25)
ax.set_aspect(1)

def lseg(pA, pB, **karg):
    a = np.array([pA, pB])
    plt.plot(*a.T, **karg)

def polyline(ps, **karg):
    for pA,pB in zip(ps[:-1], ps[1:]):
        lseg(pA, pB, **karg)

def polygon(ps, **karg):
    polyline(ps+ps[:1], **karg)

def mtanAB(tA, tB):
    return (tA+tB)/(1-tA*tB)

def dsqr(v):
    x,y = v
    return x**2+y**2

def vlen(v):
    return sympy.sqrt(dsqr(v))

def vunit(v):
    vl = vlen(v)
    x,y = v
    return Point(x/vl, y/vl)

def urotate(v, u):
    x,y = v
    c,s = u
    return Point(x*c-y*s, x*s+y*c)

def dot(v1, v2):
    x1,y1 = v1
    x2,y2 = v2
    return x1*x2+y1*y2

def ca2(v1, v2):
    return dot(v1,v2)**2/(dsqr(v1)*dsqr(v2))

def linemp(m, p):
    x,y = p
    b = y-m*x
    return Point(m,b)

def slope(pA, pB):
    xA,yA = pA
    xB,yB = pB
    m = (yB-yA)/(xB-xA)
    return m

def line2p(pA, pB):
    return linemp(slope(pA,pB),pA)

def lineabc(a, b, c):
    return Point(-a/b, -c/b)

def x2line(l1, l2):
    m1,b1 = l1
    m2,b2 = l2
    x = (b1-b2)/(m2-m1)
    y = m1*x+b1
    return Point(x,y)

def findC(pA, pB, tA, tB):
    x = tB/(tA+tB)
    y = x*tA
    vAB = pB-pA
    v = vlen(vAB)*Point(x,y)
    return pA+urotate(v, vunit(vAB))

def orthocenter(pA, pB, tA, tB):
    x = tB/(tA+tB)
    y = x/tB
    vAB = pB-pA
    v = vlen(vAB)*Point(x,y)
    return pA+urotate(v, vunit(vAB))

def orthocenter2(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1/m1, pB)
    l2 = linemp(-1/m2, pA)
    return x2line(l1,l2)

def circumcircle(pA, pB, tA, tB):
    x = sympy.Rational(1,2)
    y = -x/mtanAB(tA,tB)
    vAB = pB-pA
    v = vlen(vAB)*Point(x,y)
    return pA+urotate(v, vunit(vAB)), vlen(v)

def circumcenter2(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1/m1, sympy.Rational(1,2)*(pA+pC))
    l2 = linemp(-1/m2, sympy.Rational(1,2)*(pB+pC))
    return x2line(l1,l2)

def xlinecircle(l, pO, rO):
    m,b = l
    xO,yO = pO
    a = 1+m**2
    b = -2*xO+2*m*(b-yO)
    c = xO**2+(b-yO)**2-rO**2
    d = b**2-4*a*c
    x1 = (-b+sympy.sqrt(d))/(2*a)
    x2 = (-b-sympy.sqrt(d))/(2*a)
    y1 = m*x1+b
    y2 = m*x2+b
    return Point(x1,y1), Point(x2,y2)

def dropfoot(p, l):
    m,b = l
    hl = linemp(-1/m, p)
    return x2line(l, hl)

def run1():
    x,y = sympy.symbols("x y")
    pA = Point(0,0)
    pB = Point(1,0)
    taddx = mtanAB(x, y)
    tadd2x = simplify(taddx.subs(y,taddx))
    tadd3x = simplify(tadd2x.subs(y,taddx))
    tadd4x = simplify(tadd3x.subs(y,taddx))
    aA = simplify(tadd4x.subs(y,x))
    aX = simplify(tadd3x.subs(y,x))
    aY = simplify(tadd2x.subs(y,x))
    pC = findC(pA, pB, aA, aA).simplify()
    pD = findC(pA, pB, aX, aY).simplify()
    pH = dropfoot(pD, line2p(pA,pC)).simplify()
    dCH2 = simplify(dsqr(pC-pH))
    dDH2 = simplify(dsqr(pD-pH))
    t2 = simplify(dDH2/dCH2)
    target = simplify(taddx.subs(y,x)**2)
    constraint = simplify(tadd4x.subs(y,x)*tadd3x.subs(y,x))
    a,b = sympy.fraction(constraint)
    p = sympy.poly(sympy.expand(b-a))
    print(p)
    a,b = sympy.fraction(tadd2x.subs(y,x)**2)
    q = sympy.poly(sympy.expand(3*a-b))
    print(q)
    s,t,h = sympy.gcdex(p,q)
    print(s)
    print(t)
    print(h)
    print(simplify(p/q))
    result = sympy.solvers.solve(3*x**3 - 27*x**2 + 33*x - 1, x)
    for e in result:
        print(e, e.evalf(subs={x:1}))
    result = sympy.solvers.nsolve(3*x**3 - 27*x**2 + 33*x - 1, x, 0)

    t2_num,t2_den = sympy.fraction(t2)
    target_num,target_den = sympy.fraction(target)
    tx = sympy.expand(t2_num*target_den-t2_den*target_num)
    qotient = sympy.cancel(tx/constraint)

def test():
    x = np.linspace(-2,8,100)
    y = 3*x**3 - 27*x**2 + 33*x - 1
    plt.plot(x,y)
    x = np.array([x[0],x[-1]])
    y = np.array([0,0])
    print(math.tan(10/180*math.pi)**2)
    plt.plot(x,y)
    plt.show()

def run2():
    x,y,z = sympy.symbols("x y z")
    pA = Point(0,0)
    pB = Point(1,0)
    pC = findC(pA, pB, x, y).simplify()
    pO = circumcenter2(pA, pB, pC)
    mO = simplify(pO.y/pO.x)
    pP = Point(z,mO*z)
    consa,consb = sympy.fraction(simplify(ca2(pC-pP,pA-pP)/ca2(pB-pP,pA-pP)))
    qa = consa-consb
    c = qa.subs(z,0)
    qa = simplify((qa-c)/z)
    b = qa.subs(z,0)
    a = simplify((qa-b)/z)
    print(a)
    print(b)
    print(c)
    d = simplify(b**2-4*a*c)
    print(d)
    print(simplify(b/a/2))
    print(simplify(c/a))
    cf = simplify(a/(2*x*y**3 + 2*x*y + y**4 - 1))
    print(simplify(a/cf))
    print(simplify(b/cf))
    print(simplify(c/cf))
    #print(consb)
    pD = dropfoot(pA, line2p(pB,pC))
    pE = dropfoot(pP, line2p(pA,pC))
    pF = Point(pP.x,0)
    d2DE = dsqr(pE-pD)
    d2PF = pP.y**2
    tga,tgb = sympy.fraction(simplify(d2DE/d2PF))
    qa = tga-tgb
    c = qa.subs(z,0)
    qa = simplify((qa-c)/z)
    b = qa.subs(z,0)
    a = simplify((qa-b)/z)
    print(a)
    print(b)
    print(c)
    
    d2PE = dsqr(pE-pP)
    d2DF = dsqr(pF-pD)
    tga,tgb = sympy.fraction(simplify(d2PE/d2DF))
    #print(tga-tgb)

def run3():
    x,y,z = sympy.symbols("x y z")
    pA = Point(0,0)
    pB = Point(1,0)
    pC = findC(pA, pB, x, y).simplify()
    pH = orthocenter2(pA, pB, pC).simplify()
    pF = Point(pC.x,0).simplify()
    pD = dropfoot(pA, line2p(pB,pC)).simplify()
    pE = dropfoot(pB, line2p(pA,pC)).simplify()
    pP = x2line(line2p(pD,pF), line2p(pA,pC)).simplify()
    lBE = line2p(pB, pE).simplify()
    lHQ = linemp(-1/lBE.x, pH).simplify()
    pQ = x2line(lHQ, line2p(pB,pP)).simplify()
    pR = x2line(lHQ, line2p(pD,pP)).simplify()
    d2QR = simplify(dsqr(pQ-pR))
    d2HR = simplify(dsqr(pH-pR))
    print("C:",pC)
    print("H:",pH)
    print("D:",pD)
    print("E:",pE)
    print("F:",pF)
    print("P:",pP)
    print("Q:",pQ)
    print("R:",pR)
    print("QR^2:",d2QR)
    print("HR^2:",d2HR)
    #aq,bq = sympy.fraction(simplify(d2QR/d2HR))
    subs = {x: math.tan(60/180*math.pi), y: math.tan(50/180*math.pi)}
    print("QR numeric:",sympy.sqrt(d2QR).evalf(subs=subs))
    print("HR numeric:",sympy.sqrt(d2HR).evalf(subs=subs))

    fpA = np.array([0,0])
    fpB = np.array([1,0])
    fpC = pC.evalf(subs)
    fpF = np.array([pF.x.evalf(subs=subs),0])
    fpE = pE.evalf(subs)
    fpD = pD.evalf(subs)
    fpP = pP.evalf(subs)
    fpH = pH.evalf(subs)
    fpQ = pQ.evalf(subs)
    fpR = pR.evalf(subs)
    polygon([fpA,fpB,fpC], color="k")
    lseg(fpA,fpD, color="k")
    lseg(fpB,fpE, color="k")
    lseg(fpC,fpF, color="k")
    lseg(fpA,fpP, color="k")
    lseg(fpP,fpB, color="k")
    lseg(fpP,fpD, color="k")
    lseg(fpR,fpQ, color="r")
    lseg(fpH,fpR, color="g")
    #ax.add_patch(matplotlib.patches.Circle(fpD, .1, color="g", fill=False))
    plt.show()

def run4():
    x,y,u,v = sympy.symbols("x y u v")
    pA = Point(0,0)
    pB = Point(1,0)
    pC = findC(pA, pB, x, y).simplify()
    pF = Point(u,0).simplify()
    pD = findC(pA, pB, v, y).simplify()
    #pH = orthocenter2(pA, pB, pC).simplify()
    pH = x2line(line2p(pA,pD), line2p(pC,pF)).simplify()
    pE = x2line(line2p(pA,pC), line2p(pB,pH)).simplify()
    #pD = dropfoot(pA, line2p(pB,pC)).simplify()
    #pE = dropfoot(pB, line2p(pA,pC)).simplify()
    pP = x2line(line2p(pD,pF), line2p(pA,pC)).simplify()
    #lBE = line2p(pB, pE).simplify()
    #lHQ = linemp(-1/lBE.x, pH).simplify()
    lAC = line2p(pA,pC).simplify()
    lHQ = linemp(lAC.x, pH).simplify()
    pQ = x2line(lHQ, line2p(pB,pP)).simplify()
    pR = x2line(lHQ, line2p(pD,pP)).simplify()
    d2QR = simplify(dsqr(pQ-pR))
    d2HR = simplify(dsqr(pH-pR))
    print("C:",pC)
    print("H:",pH)
    print("D:",pD)
    print("E:",pE)
    print("F:",pF)
    print("P:",pP)
    print("Q:",pQ)
    print("R:",pR)
    print("QR^2:",d2QR)
    print("HR^2:",d2HR)
    #aq,bq = sympy.fraction(simplify(d2QR/d2HR))
    subs = {x: math.tan(60/180*math.pi), y: math.tan(50/180*math.pi), u: 0.75, v: math.tan(40/180*math.pi)}
    print("QR numeric:",sympy.sqrt(d2QR).evalf(subs=subs))
    print("HR numeric:",sympy.sqrt(d2HR).evalf(subs=subs))

    fpA = np.array([0,0])
    fpB = np.array([1,0])
    fpC = pC.evalf(subs)
    fpF = np.array([pF.x.evalf(subs=subs),0])
    fpE = pE.evalf(subs)
    fpD = pD.evalf(subs)
    fpP = pP.evalf(subs)
    fpH = pH.evalf(subs)
    fpQ = pQ.evalf(subs)
    fpR = pR.evalf(subs)
    polygon([fpA,fpB,fpC], color="k")
    lseg(fpA,fpD, color="k")
    lseg(fpB,fpE, color="k")
    lseg(fpC,fpF, color="k")
    lseg(fpA,fpP, color="k")
    lseg(fpP,fpB, color="k")
    lseg(fpP,fpD, color="k")
    lseg(fpF,fpD, color="k")
    lseg(fpR,fpQ, color="r")
    lseg(fpH,fpR, color="g")
    #ax.add_patch(matplotlib.patches.Circle(fpD, .1, color="g", fill=False))
    plt.show()

def run5():
    x,y,z = sympy.symbols("x y z")
    pA = Point(0,0)
    pB = Point(1,0)
    pC = findC(pA, pB, x, y).simplify()
    pP = findC(pA, pB, (x-z)/(1+x*z), (y+z)/(1-y*z)).simplify()
    pD = dropfoot(pP, line2p(pB,pC)).simplify()
    pE = dropfoot(pP, line2p(pA,pC)).simplify()
    pF = Point(pP.x,0)
    pH = orthocenter2(pA, pB, pC).simplify()
    lDE = line2p(pD,pE).simplify()
    lDF = line2p(pD,pF).simplify()
    pX = x2line(lDE, line2p(pP, pH)).simplify()
    m1 = simplify(lDE.x)
    m2 = simplify(lDF.x)
    print(m1)
    print(m2)
    d2XH = simplify(dsqr(pX-pH))
    d2XP = simplify(dsqr(pX-pP))
    print(d2XH)
    print(d2XP)
    subs = {x: math.tan(60/180*math.pi), y: math.tan(50/180*math.pi), z: math.tan(15/180*math.pi)}

    fpA = np.array([0,0])
    fpB = np.array([1,0])
    fpC = pC.evalf(subs)
    fpD = pD.evalf(subs)
    fpE = pE.evalf(subs)
    fpF = np.array([pF.x.evalf(subs=subs),0])
    fpP = pP.evalf(subs)
    fpH = pH.evalf(subs)
    fpX = pX.evalf(subs)
    polygon([fpA,fpB,fpC], color="k")
    lseg(fpD,fpP, color="k")
    lseg(fpE,fpP, color="k")
    lseg(fpF,fpP, color="k")
    lseg(fpC,fpE, color="k")
    lseg(fpF,fpD)
    lseg(fpE,fpD)
    lseg(fpX,fpH)
    lseg(fpX,fpP)
    #ax.add_patch(matplotlib.patches.Circle(fpD, .1, color="g", fill=False))
    plt.show()

def run():
    """
    Feurbach Theorem
    """
    x,y = sympy.symbols("x y")
    #sin2x/cos2x=2sinxcosx/(cosxcosx-sinxsinx)=2tanx/(1-tanxtanx)
    tan2x = 2*x/(1-x**2)
    tan2y = 2*y/(1-y**2)
    pA = Point(0,0)
    pB = Point(1,0)
    pC = findC(pA, pB, tan2x, tan2y).simplify()
    phAB = sympy.Rational(1,2) *(pA+pB)
    phBC = sympy.Rational(1,2) * (pC+pB)
    phAC =  sympy.Rational(1,2) * (pA+pC)
    pnO = circumcenter2(phBC, phAC, phAB).simplify()
    rnO2 = simplify(dsqr(phAB-pnO))
    p2O = findC(pA, pB, x, y).simplify()
    p2B = dropfoot(p2O, line2p(pA,pC)).simplify()
    r2O2 = simplify(dsqr(p2O-p2B))
    pxA = findC(pA, pB, x, -1/y).simplify()
    rxA = pxA.y
    d2xA = simplify(dsqr(pnO-pxA))
    lhs = simplify((d2xA-rxA**2-rnO2)**2)
    rhs = simplify(4*rxA**2*rnO2)
    print(lhs)
    print(rhs)
    d2OO = simplify(dsqr(pnO-p2O))
    lhs2 = simplify((r2O2+rnO2-d2OO)**2)
    rhs2 = simplify(4*r2O2*rnO2)
    print(lhs2)
    print(rhs2)
    subs = {x: math.tan(60/2/180*math.pi), y: math.tan(40/2/180*math.pi)}

    fpA = np.array([0,0])
    fpB = np.array([1,0])
    fpC = pC.evalf(subs)
    fphAB = phAB.evalf(subs)
    fphAC = phAC.evalf(subs)
    fphBC = phBC.evalf(subs)
    fpnO = pnO.evalf(subs)
    fp2O = p2O.evalf(subs)
    fp2B = p2B.evalf(subs)
    fpxA = pxA.evalf(subs)
    frnO2 = math.sqrt(rnO2.evalf(subs=subs))
    fr2O2 = math.sqrt(r2O2.evalf(subs=subs))
    frxA = rxA.evalf(subs=subs)
    polygon([fpA,fpB,fpC], color="k")
    ax.add_patch(matplotlib.patches.Circle(fpnO, frnO2, color="g", fill=False))
    ax.add_patch(matplotlib.patches.Circle(fp2O, fr2O2, color="r", fill=False))
    ax.add_patch(matplotlib.patches.Circle(fpxA, frxA, color="g", fill=False))
    plt.show()

if __name__ == "__main__":
    run()