import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import sympy

fig, ax = plt.subplots()
plt.xlim(-.5,2.)
plt.ylim(-.6,1.25)
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

def plus(pA, pB):
    xA,yA = pA
    xB,yB = pB
    return (xA+xB, yA+yB)

def minus(pA, pB):
    xA,yA = pA
    xB,yB = pB
    return (xA-xB, yA-yB)

def smul(s, v):
    x,y = v
    return (s*x, s*y)

def dsqr(v):
    x,y = v
    return x**2+y**2

def vlen(v):
    return sympy.sqrt(dsqr(v))

def vunit(v):
    vl = vlen(v)
    x,y = v
    return (x/vl, y/vl)

def urotate(v, u):
    x,y = v
    c,s = u
    return (x*c-y*s, x*s+y*c)

def linemp(m, p):
    x,y = p
    b = y-m*x
    return (m,b)

def slope(pA, pB):
    xA,yA = pA
    xB,yB = pB
    m = (yB-yA)/(xB-xA)
    return m

def line2p(pA, pB):
    return linemp(slope(pA,pB),pA)

def lineabc(a, b, c):
    return (-a/b, -c/b)

def x2line(l1, l2):
    m1,b1 = l1
    m2,b2 = l2
    x = (b1-b2)/(m2-m1)
    y = m1*x+b1
    return (x,y)

def findC(pA, pB, tA, tB):
    x = tB/(tA+tB)
    y = x*tA
    vAB = minus(pB,pA)
    v = smul(vlen(vAB), (x,y))
    return plus(pA,urotate(v, vunit(vAB)))

def orthocenter(pA, pB, tA, tB):
    x = tB/(tA+tB)
    y = x/tB
    vAB = minus(pB,pA)
    v = smul(vlen(vAB), (x,y))
    return plus(pA,urotate(v, vunit(vAB)))

def orthocenter2(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1/m1, pB)
    l2 = linemp(-1/m2, pA)
    return x2line(l1,l2)

def circumcircle(pA, pB, tA, tB):
    x = sympy.Rational(1,2)
    y = -x/mtanAB(tA,tB)
    vAB = minus(pB,pA)
    v = smul(vlen(vAB), (x,y))
    return plus(pA,urotate(v, vunit(vAB))), vlen(v)

def circumcenter2(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1/m1, smul(sympy.Rational(1,2),plus(pA,pC)))
    l2 = linemp(-1/m2, smul(sympy.Rational(1,2),plus(pB,pC)))
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
    return ((x1,y1),(x2,y2))

def dropfoot(p, l):
    m,b = l
    hl = linemp(-1/m, p)
    return x2line(l, hl)

def eqgrav(pA, rA, pB, rB):
    return plus(pA,smul(rA/(rA+rB),minus(pB,pA)))
def mcancel(p):
    return (sympy.cancel(p[0]),sympy.cancel(p[1]))

def m1cancel(a):
    return sympy.cancel(a)

def m2eval(p, subs):
    return np.array([p[0].evalf(subs=subs),p[1].evalf(subs=subs)])

def m1eval(a, subs):
    return a.evalf(subs=subs)

def run():
    #tangent of angle A and B in a triangle
    tA,tB = sympy.symbols("x y")

    pA = 0,0
    pB = 1,0

    pC = mcancel(findC(pA, pB, tA, tB))
    pO,rO = circumcircle(pA, pB, tA, tB)
    pO = mcancel(pO)
    rO = m1cancel(rO)
    pH = mcancel(orthocenter(pA, pB, tA, tB))
    pD = mcancel(smul(sympy.Rational(1,2),plus(pA,pC)))
    pE = mcancel(smul(sympy.Rational(1,2),plus(pB,pC)))
    pF = mcancel(dropfoot(pB, line2p(pD, pH)))
    pG = mcancel(dropfoot(pA, line2p(pE, pH)))
    pI = mcancel(x2line(line2p(pD, pE), line2p(pF, pG)))
    pX = mcancel(dropfoot(pI, line2p(pB, pC)))
    tu = sympy.cancel(dsqr(minus(pI,pX)))
    td = sympy.cancel(dsqr(minus(pC,pX)))
    target = sympy.cancel(tu/td)
    #for e in [pC,pO,rO,pH,pD,pE,pF,pG,pI,pX]: print(e)
    print("pA:", pC)
    print("pO:", pO)
    print("rO:", rO)
    print("pH:", pH)
    print("pD:", pD)
    print("pE:", pE)
    print("pF:", pF)
    print("pG:", pG)
    print("pI:", pI)
    print("pX:", pX)
    print(target)

    subs={tA:math.tan(math.pi*40/180), tB:math.tan(math.pi*75/180)}
    pA=np.array(pA)
    pB=np.array(pB)
    pC=m2eval(pC, subs)
    pO = m2eval(pO, subs)
    rO = m1eval(rO, subs)
    pD=m2eval(pD, subs)
    pE=m2eval(pE, subs)
    pH=m2eval(pH, subs)
    pF=m2eval(pF, subs)
    pG=m2eval(pG, subs)
    pI=m2eval(pI, subs)

    polygon([pA,pB,pC], color="k")
    #ax.add_patch(matplotlib.patches.Circle(pO, rO, color="g", fill=False))
    lseg(pD, pI)
    lseg(pD, pF)
    lseg(pE, pG)
    lseg(pA, pG)
    lseg(pG, pI)
    lseg(pC, pI)
    plt.show()

def test():
    x,y,z = sympy.symbols("x y z")
    pA = 0,0
    pB = 1,0
    pC = findC(pA,pB,1/z,1/z)
    pD = findC(pA,pB,x,y)
    d2AD = dsqr(minus(pD, pA))
    d2AC = dsqr(minus(pC, pA))
    d2CD = dsqr(minus(pC, pD))
    pX = dropfoot(pD, line2p(pA, pC))
    d2DX = dsqr(minus(pD, pX))
    d2CX = dsqr(minus(pC, pX))
    t2 = m1cancel(d2DX/d2CX)
    print(pC,pD)
    print(d2AD)
    print(d2AC)
    print(d2CD)
    print(pX)
    print(t2)

def test2():
    x,y,r,x1,y1,r1,x2,y2,r2,x3,y3,r3 = sympy.symbols("x y r x1 y1 r1 x2 y2 r2 x3 y3 r3")
    R2 = sympy.Rational(2)
    l12 = lineabc(R2*(x2-x1), R2*(y2-y1), -R2*(r1-r2)*r-r1**2+r2**2+x1**2-x2**2+y1**2-y2**2)
    l32 = lineabc(R2*(x2-x3), R2*(y2-y3), -R2*(r3-r2)*r-r3**2+r2**2+x3**2-x2**2+y3**2-y2**2)
    #print(l12[1].subs(r,0))
    xr,yr = x2line(l12, l32)
    xr = m1cancel(xr)
    yr = m1cancel(yr)
    xr = m1cancel(xr.subs([(x1,0),(y1,0),(x2,1),(y2,0)]))
    yr = m1cancel(yr.subs([(x1,0),(y1,0),(x2,1),(y2,0)]))
    #print(xr)
    #print(yr)
    #e = (xr-x1)**2+(yr-y1)**2-(r+r1)**2
    e = xr**2+yr**2-(r+r1)**2
    e = m1cancel(e)
    #print(e)
    e = m1cancel(e.subs([(x1,0),(y1,0),(x2,1),(y2,0)]))
    p = e
    q = m1cancel(p.subs(r,0))
    c = q
    p = m1cancel((p-q)/r)
    q = m1cancel(p.subs(r,0))
    b = q
    p = m1cancel((p-q)/r)
    q = m1cancel(p.subs(r,0))
    a = q
    #print(a)
    #print(b)
    #print(c)
    d = m1cancel(b**2-R2*R2*a*c)
    #print(d)
    u = -b/(R2*a)
    v = sympy.sqrt(d)/(R2*a)
    rr1,rr2 = u+v,u-v
    #print(rr1)

    pA = 0,0
    pB = 1,0
    tA,tB = sympy.symbols("tA tB")
    pC = findC(pA, pB, tA, tB)
    subs={tA:math.tan(math.pi*40/180), tB:math.tan(math.pi*75/180)}
    pC = m2eval(pC, subs)
    xC,yC = pC
    r1v = .3
    r2v = .2
    r3v= .15
    subs2 = {
        x3:xC,
        y3:yC,
        r1:r1v,
        r2:r2v,
        r3:r3v,
        }

    pA = np.array(pA)
    pB = np.array(pB)
    polygon([pA,pB,pC], color="k")
    ax.add_patch(matplotlib.patches.Circle(pA, r1v, color="k", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pB, r2v, color="k", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pC, r3v, color="k", fill=False))

    rr = rr1
    pO = xr.subs(r,rr),yr.subs(r,rr)
    pO = m2eval(pO, subs2)
    rO = m1eval(rr, subs2)
    print(pO,rO)
    ax.add_patch(matplotlib.patches.Circle(pO, rO, color="g", fill=False))

    rr = rr2
    pO = xr.subs(r,rr),yr.subs(r,rr)
    pO = m2eval(pO, subs2)
    rO = m1eval(rr, subs2)
    print(pO,rO)
    ax.add_patch(matplotlib.patches.Circle(pO, rO, color="r", fill=False))

    plt.show()


if __name__ == "__main__":
    test2()