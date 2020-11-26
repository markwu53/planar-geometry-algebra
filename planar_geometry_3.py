import numpy as np
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
#plt.xlim(-1.,1.5)
#plt.ylim(-1.5,1.5)
ax.set_aspect(1)

def lseg(pA, pB, **karg):
    a = np.array([pA, pB])
    plt.plot(*a.T, **karg)

def polyline(ps, **karg):
    for pA,pB in zip(ps[:-1], ps[1:]):
        lseg(pA, pB, **karg)

def polygon(ps, **karg):
    polyline(ps+ps[:1], **karg)

def rad(degree): return np.pi*degree/180
def dot(u, v): return u[0]*v[0]+u[1]*v[1]
def dsqr(v): return dot(v, v)
def vlen(v): return np.sqrt(dsqr(v))
def linemp(m, p): return np.array([m, p[1]-m*p[0]])
def slope(p, q): return (q[1]-p[1])/(q[0]-p[0])
def line2p(p, q): return linemp(slope(p,q),p)

def x2line(l1, l2):
    m1,b1 = l1
    m2,b2 = l2
    x = (b1-b2)/(m2-m1)
    y = m1*x+b1
    return np.array([x,y])

def findC(pA, pB, tA, tB):
    x = tB / (tA + tB)
    y = x * tA
    p1 = pB - pA
    p2 = np.array([-p1[1], p1[0]])
    return pA + x * p1 + y * p2

def orthocenter(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    return x2line(linemp(-1/m1, pB), linemp(-1/m2, pA))

def circumcenter(pA, pB, pC):
    m1 = slope(pA, pC)
    m2 = slope(pB, pC)
    l1 = linemp(-1/m1, 1/2*(pA+pC))
    l2 = linemp(-1/m2, 1/2*(pB+pC))
    return x2line(l1,l2)

def shift(d):
    fontsize = .014
    sdir = [[0,0],[-1,0],[-2,0],[-2,-1],[-2,-2],[-1,-2],[0,-2],[0,-1]]
    return np.array(sdir[d])*fontsize + np.array([fontsize/2,fontsize/2])

def run():
    pA = np.array([0,0])
    pB = np.array([1,0])
    aA = rad(50)
    aB = rad(70)
    pC = findC(pA, pB, np.tan(aA), np.tan(aB))
    pH = orthocenter(pA, pB, pC)
    pD = x2line(line2p(pA, pH), line2p(pB, pC))
    pE = x2line(line2p(pB, pH), line2p(pA, pC))
    pF = x2line(line2p(pC, pH), line2p(pB, pA))

    polygon([pA, pB, pC])
    lseg(pA, pD)
    lseg(pB, pE)
    lseg(pC, pF)
    polygon([pD, pE, pF])
    ax.add_patch(matplotlib.patches.Circle(pH, .005, color="r", fill=False))
    ax.text(*(pA+shift(4)), "A", fontweight="normal", fontsize=14)
    ax.text(*(pB+shift(6)), "B", fontweight="normal", fontsize=14)
    ax.text(*(pC+shift(1)), "C", fontweight="normal", fontsize=14)
    ax.text(*(pD+shift(0)), "D", fontweight="normal", fontsize=14)
    ax.text(*(pE+shift(2)), "E", fontweight="normal", fontsize=14)
    ax.text(*(pF+shift(5)), "F", fontweight="normal", fontsize=14)
    ax.text(*(pH+shift(1)), "H", fontweight="normal", fontsize=14)

    plt.show()

if __name__ == "__main__":
    run()
