import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
#from sympy import *
import sympy

def fadd(a, b):
    a1, a2 = a
    b1, b2 = b
    return (a1*b2+a2*b1,a2*b2)

def fmin(a, b):
    a1, a2 = a
    b1, b2 = b
    return (a1*b2-a2*b1,a2*b2)

def fmul(a, b):
    a1, a2 = a
    b1, b2 = b
    return (a1*b1, a2*b2)

def fdiv(a, b):
    a1, a2 = a
    b1, b2 = b
    return (a1*b2, a2*b1)

def finv(a):
    a1, a2 = a
    return (a2,a1)

def fneg(a):
    a1, a2 = a
    return (-a1,a2)

def fsqrt(a):
    a1, a2 = a
    return (sympy.sqrt(a1),sympy.sqrt(a2))

def mtanAB(tA, tB):
    return fdiv(fadd(tA,tB),fmin((1,1),fmul(tA,tB)))

def plus(pA, pB):
    xA,yA = pA
    xB,yB = pB
    return (fadd(xA,xB), fadd(yA,yB))

def minus(pA, pB):
    xA,yA = pA
    xB,yB = pB
    return (fmin(xA,xB), fmin(yA,yB))

def smul(s, v):
    x,y = v
    return (fmul(s,x), fmul(s,y))

def dsqr(v):
    x,y = v
    return fadd(fmul(x,x),fmul(y,y))

def vlen(v):
    return fsqrt(dsqr(v))

def vunit(v):
    vl = vlen(v)
    x,y = v
    return (fdiv(x,vl), fdiv(y,vl))

def urotate(v, u):
    x,y = v
    c,s = u
    return (fmin(fmul(x,c),fmul(y,s)), fadd(fmul(x,s),fmul(y,c)))

def linemp(m, p):
    x,y = p
    b = fmin(y,fmul(m,x))
    return (m,b)

def slope(pA, pB):
    xA,yA = pA
    xB,yB = pB
    return fdiv(fmin(yB,yA),fmin(xB,xA))

def line2p(pA, pB):
    return linemp(slope(pA,pB),pA)

def x2line(l1, l2):
    m1,b1 = l1
    m2,b2 = l2
    x = fdiv(fmin(b1,b2),fmin(m2,m1))
    y = fadd(fmul(m1,x),b1)
    return (x,y)

def findC(pA, pB, tA, tB):
    x = fdiv(tB,fadd(tA,tB))
    y = fmul(x,tA)
    vAB = minus(pB,pA)
    v = smul(vlen(vAB), (x,y))
    return plus(pA,urotate(v, vunit(vAB)))

def orthocenter(pA, pB, tA, tB):
    x = fdiv(tB,fadd(tA,tB))
    y = fdiv(x,tA)
    vAB = minus(pB,pA)
    v = smul(vlen(vAB), (x,y))
    return plus(pA,urotate(v, vunit(vAB)))

def orthocenter2(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(finv(m1), pB)
    l2 = linemp(finv(m2), pA)
    return x2line(l1,l2)

def circumcircle(pA, pB, tA, tB):
    x = (1,2)
    y = fneg(fdiv(x,mtanAB(tA, tB)))
    vAB = minus(pB,pA)
    v = smul(vlen(vAB), (x,y))
    return plus(pA,urotate(v, vunit(vAB))), vlen(v)

def circumcenter2(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(finv(m1), smul((1,2),plus(pA,pC)))
    l2 = linemp(finv(m2), smul((1,2),plus(pB,pC)))
    return x2line(l1,l2)

def xlinecircle(l, pO, rO):
    m,b = l
    xO,yO = pO
    a = fadd((1,1),fmul(m,m))
    b = fadd(fneg((2,1),xO),fmul(fmul((2,1),m),fmin(b,yO)))
    c = fmin(fadd(fmul(xO,xO),fmul(fmin(b,yO),fmin(b,yO))),fmul(rO,rO))
    d = fmin(fmul(b,b),fmul(fmul((4,1),a),c))
    x1 = fdiv(fadd(fneg(b),fsqrt(d)),fmul((2,1),a))
    x2 = fdiv(fmin(fneg(b),fsqrt(d)),fmul((2,1),a))
    y1 = fadd(fmul(m,x1),b)
    y2 = fadd(fmul(m,x2),b)
    return ((x1,y1),(x2,y2))

def dropfoot(p, l):
    m,b = l
    hl = linemp(finv(m), p)
    return x2line(l, hl)

def run():
    #tangent of angle A and B in a triangle
    #tA,tB = sympy.symbols("tA tB")
    x,y = sympy.symbols("x y")
    tA = (x,1)
    tB = (y,1)

    pA = ((0,1),(0,1))
    pB = ((1,1),(0,1))
    #xA,yA,xB,yB = sympy.symbols("xA yA xB yB")
    #pA = (xA, yA)
    #pB = (xB, yB)

    pC = findC(pA, pB, tA, tB)
    print(pC)
    print(sympy.cancel(pC[0][0]/pC[0][1]))
    pO,rO = circumcircle(pA, pB, tA, tB)
    #print(pO)
    #*print(rO)
    pH = orthocenter(pA, pB, tA, tB)
    #print(pH)
    pD = smul((1,2),plus(pA,pC))
    #print(pD)
    pE = smul((1,2),plus(pB,pC))
    #pF,pF2 = xlinecircle(line2p(pD, pH), pO, rO)
    pF = dropfoot(pB, line2p(pD, pH))
    #print(pF)
    #pG,pG2 = xlinecircle(line2p(pE, pH), pO, rO)
    pG = dropfoot(pA, line2p(pE, pH))
    pI = x2line(line2p(pD, pE), line2p(pF, pG))
    #print(pI)
    pX = dropfoot(pI, line2p(pB, pC))
    #print(pX)
    #tsqr = fdiv(dsqr(minus(pI,pX)),dsqr(minus(pC,pX)))
    #t1,t2 = tsqr
    #print(sympy.expand(t1))
    #target = fdiv(vlen(minus(pI,pX)),vlen(minus(pC,pX)))
    #print(sympy.expand(target))
    #print(target)

def test():
    x = (1,1)
    y = (0,1)
    v = (x,y)
    print(v)
    #r = vlen(v)
    x,y = v
    print(fsqrt(fadd(fmul(x,x),fmul(y,y))))
    #print(fmul(x,x))
    #print(fmul(y,y))
    #print(fadd(fmul(x,x),fmul(y,y)))

if __name__ == "__main__":
    run()