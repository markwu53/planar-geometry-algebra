import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import sympy
import numpy as np

fig, ax = plt.subplots()
plt.xlim(-2.,2.5)
plt.ylim(-2.5,2.5)
ax.set_aspect(1)

def lseg(pA, pB, **karg):
    a = np.array([pA, pB])
    plt.plot(*a.T, **karg)

def polyline(ps, **karg):
    for pA,pB in zip(ps[:-1], ps[1:]):
        lseg(pA, pB, **karg)

def polygon(ps, **karg):
    polyline(ps+ps[:1], **karg)

def rad(degree):
    return math.pi*degree/180

def dsqr(v):
    x,y = v
    return x**2+y**2

def vlen(v):
    return math.sqrt(dsqr(v))

def vunit(v):
    return v/vlen(v)

def urotate(v, u):
    x,y = v
    c,s = u
    return np.array([x*c-y*s, x*s+y*c])

def rotate(v, a):
    x,y = v
    c = math.cos(a)
    s = math.sin(a)
    return np.array([x*c-y*s, x*s+y*c])

def dot(vA, vB):
    xA,yA = vA
    xB,yB = vB
    return xA*xB+yA*yB

def vangle(vA, vB):
    va = vlen(vA)
    vb = vlen(vB)
    c = dot(vA, vB)/(va*vb)
    return math.acos(c)

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
    return np.array([x,y])

def vfoot(pX, pA, pB):
    AB = vlen(pB-pA)
    BC = vlen(pX-pB)
    CA = vlen(pA-pX)
    x = (CA**2-BC**2)/AB**2
    qA = .5*(1+x)
    qB = .5*(1-x)
    return qB*pA+qA*pB

def vfoot3(pA, pB, pC):
    phA = vfoot(pA, pB, pC)
    phB = vfoot(pB, pC, pA)
    phC = vfoot(pC, pA, pB)
    return phA,phB,phC

def pldist(pX, pA, pB):
    pH = vfoot(pX, pA, pB)
    return vlen(pX-pH)

def incircle(pA, pB, pC):
    l1 = bisect(pA, pB-pA, pC-pA)
    l2 = bisect(pB, pC-pB, pA-pB)
    pI = x2line(l1, l2)
    rI = pldist(pI, pA, pB)
    return pI,rI

def circumcircle(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1/m1, 1/2*(pA+pC))
    l2 = linemp(-1/m2, 1/2*(pB+pC))
    pO = x2line(l1,l2)
    rO = vlen(pO-pA)
    return pO,rO

def findC(pA, pB, tA, tB):
    x = tB/(tA+tB)
    y = x*tA
    vAB = pB-pA
    v = vlen(vAB)*np.array([x,y])
    return pA+urotate(v, vunit(vAB))

def bisect(pA, v1, v2):
    a = vangle(v1, v2)
    pM = pA+rotate(v1, a/2)
    return line2p(pA,pM)

def excircle(pA, pB, rA, rB):
    a = math.acos((rA+rB)/vlen(pB-pA))
    uvAB = vunit(pB-pA)

    ptA = pA+rotate(rA*uvAB,a)
    ptB = pB+rotate(-rB*uvAB,a)
    l1 = line2p(ptA,ptB)
    ptA = pA+rotate(rA*uvAB,-a)
    ptB = pB+rotate(-rB*uvAB,-a)
    l2 = line2p(ptA,ptB)

    b = math.acos((rA-rB)/vlen(pB-pA))
    ptA = pA+rotate(rA*uvAB,b)
    ptB = pB+rotate(rB*uvAB,b)
    l3 = line2p(ptA,ptB)

    pX = x2line(l1, l2)
    pY = x2line(l2, l3)
    pZ = x2line(l1, l3)
    bl1 = bisect(pY, pY-pX, pZ-pY)
    bl2 = bisect(pZ, pY-pZ, pZ-pX)
    pO = x2line(bl1,bl2)
    pH = vfoot(pO, pY, pZ)
    rO = vlen(pO-pH)
    return pX,pY,pZ,pO,rO

def t3circle(pA, rA, pB, rB, pC, rC):
    xA,yA = pA
    xB,yB = pB
    xC,yC = pC
    #(x-xA)^2+(y-yA)^2=(r+rA)^2
    #(x-xB)^2+(y-yB)^2=(r+rB)^2
    #(x-xC)^2+(y-yC)^2=(r+rC)^2
    #A-B
    #ax+by+cr+d=0
    a1 = 2*(xB-xA)
    b1 = 2*(yB-yA)
    c1 = -2*(rA-rB)
    d1 = -(rA**2-rB**2)+(xA**2-xB**2)+(yA**2-yB**2)
    a2 = 2*(xC-xA)
    b2 = 2*(yC-yA)
    c2 = -2*(rA-rC)
    d2 = -(rA**2-rC**2)+(xA**2-xC**2)+(yA**2-yC**2)
    #express x,y as a linear function of r
    #x=xm r +xb, y=ym r +yb
    xm,xb = lineabc(c1*b2-c2*b1,a1*b2-a2*b1,d1*b2-d2*b1)
    ym,yb = lineabc(c1*a2-c2*a1,-a1*b2+a2*b1,d1*a2-d2*a1)
    #plug back in eq 1
    #ar^2+br+c=0
    a = xm**2+ym**2-1
    b = 2*xm*(xb-xA)+2*ym*(yb-yA)-2*rA
    c = (xb-xA)**2+(yb-yA)**2-rA**2
    d = b**2-4*a*c
    u = -b/(2*a)
    v = math.sqrt(d)/(2*a)
    r1,r2 = u+v,u-v
    pO1 = np.array([xm*r1+xb, ym*r1+yb])
    pO2 = np.array([xm*r2+xb, ym*r2+yb])
    return pO1,r1,pO2,r2

def run():
    """
    pA = np.array([0,0])
    pB = np.array([.2,0])
    aA = rad(15)
    aB = rad(150)
    rA = .08
    rB = .06
    rC = .12
    """
    pA = np.array([0,0])
    pB = np.array([1,0])
    aA = rad(55)
    aB = rad(75)
    rA = .4
    rB = .3
    rC = .6

    pC = findC(pA, pB, math.tan(aA), math.tan(aB))

    polygon([pA,pB,pC], color="k")
    ax.add_patch(matplotlib.patches.Circle(pA, rA, color="k", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pB, rB, color="k", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pC, rC, color="k", fill=False))

    p1,p2,p3,p4,d1 = excircle(pB,pA,rB,rA)
    polygon([p1,p2,p3], color="b")
    polyline([pA,p4,pB], color="g")
    ax.add_patch(matplotlib.patches.Circle(p4, d1, color="g", fill=False))
    peAB = p4
    p4,d1 = circumcircle(1/2*(p1+p2), 1/2*(p2+p3), 1/2*(p3+p1))
    ax.add_patch(matplotlib.patches.Circle(p4, d1, color="b", fill=False))
    p1,p2,p3,p4,d1 = excircle(pA,pC,rA,rC)
    polygon([p1,p2,p3], color="b")
    polyline([pA,p4,pC], color="g")
    ax.add_patch(matplotlib.patches.Circle(p4, d1, color="g", fill=False))
    peAC = p4
    p4,d1 = circumcircle(1/2*(p1+p2), 1/2*(p2+p3), 1/2*(p3+p1))
    ax.add_patch(matplotlib.patches.Circle(p4, d1, color="b", fill=False))
    p1,p2,p3,p4,d1 = excircle(pC,pB,rC,rB)
    polygon([p1,p2,p3], color="b")
    polyline([pB,p4,pC], color="g")
    ax.add_patch(matplotlib.patches.Circle(p4, d1, color="g", fill=False))
    peBC = p4
    p4,d1 = circumcircle(1/2*(p1+p2), 1/2*(p2+p3), 1/2*(p3+p1))
    ax.add_patch(matplotlib.patches.Circle(p4, d1, color="b", fill=False))
    #polygon([peAB,peAC,peBC], color="b")
    p1,p2,p3 = vfoot3(peAB, peBC, peAC)
    #polygon([p1,p2,p3], color="b")
    p1,d1 = incircle(p1, p2, p3)
    #ax.add_patch(matplotlib.patches.Circle(p1, d1, color="b", fill=False))

    p1,d1,p2,d2 = t3circle(pA, rA, pB, rB, pC, rC)
    #ax.add_patch(matplotlib.patches.Circle(p1, d1, color="r", fill=False))
    ax.add_patch(matplotlib.patches.Circle(p2, d2, color="r", fill=False))
    l1=linemp(-1/slope(pA,p2),pA+rA*vunit(p2-pA))
    l2=linemp(-1/slope(pB,p2),pB+rB*vunit(p2-pB))
    l3=linemp(-1/slope(pC,p2),pC+rC*vunit(p2-pC))
    p1,p2,p3 = x2line(l1, l2),x2line(l2, l3),x2line(l3, l1)
    #polygon([p1,p2,p3], color="k")

    plt.show()

if __name__ == "__main__":
    run()