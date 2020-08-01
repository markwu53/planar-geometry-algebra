import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import sympy
import numpy as np

fig, ax = plt.subplots()
plt.xlim(-.5,1.5)
plt.ylim(-.5,1.5)
ax.set_aspect(1)

def lseg(pA, pB, **karg):
    a = np.array([pA, pB])
    plt.plot(*a.T, **karg)

def polyline(ps, **karg):
    for pA,pB in zip(ps[:-1], ps[1:]):
        lseg(pA, pB, **karg)

def vlen(p):
    x,y = p
    return math.sqrt(x**2+y**2)

def vfoot(pC, pA, pB):
    AB = vlen(pB-pA)
    BC = vlen(pC-pB)
    CA = vlen(pA-pC)
    x = (CA**2-BC**2)/AB**2
    qA = .5*(1+x)
    qB = .5*(1-x)
    return qB*pA+qA*pB

def pbetween(pA, pB, dA):
    AB = vlen(pB-pA)
    dB = AB-dA
    return dA/AB*pB+dB/AB*pA

def dot(vA, vB):
    xA,yA = vA
    xB,yB = vB
    return xA*xB+yA*yB

def xpv2(pA, vA, pB, vB):
    vAB = pB-pA
    vBA = pA-pB
    AB = vlen(vAB)
    va = vlen(vA)
    vb = vlen(vB)
    cosA = dot(vAB, vA)/(AB*va)
    sinA = math.sqrt(1-cosA**2)
    cosB = dot(vBA, vB)/(AB*vb)
    sinB = math.sqrt(1-cosB**2)
    if cosA == 0:
        return pB+AB/cosB*vB/vb
    if cosB == 0:
        return pA+AB/cosA*vA/va
    tanA = sinA/cosA
    tanB = sinB/cosB
    return pA+AB*tanB/(cosA*(tanA+tanB))*vA/va

def rotate(v, a):
    x,y = v
    c = math.cos(a)
    s = math.sin(a)
    return np.array([x*c-y*s, x*s+y*c])

def rotate2(v, cosa):
    x,y = v
    if cosa == 0:
        sina = 1.
    else:
        sina = math.sqrt(1-cosa**2)
    return np.array([x*cosa-y*sina, x*sina+y*cosa])

def vrotate(v):
    x,y = v
    return np.array([-y,x])

def vunit(v):
    return v/vlen(v)

def circumcircle(pA, pB, pC):
    vCA = pA-pC
    vCB = pB-pC
    CA = vlen(vCA)
    CB = vlen(vCB)
    AB = vlen(pB-pA)
    cosC = dot(vCA, vCB)/(CA*CB)
    phAB = (pA+pB)/2
    vhAB = phAB-pA
    if cosC == 0: return phAB, AB/2
    vnormal = vrotate(vunit(vhAB))
    sinC = math.sqrt(1-cosC**2)
    tanC = sinC/cosC
    vOh = vnormal*(AB/2/tanC)
    pO = phAB+vOh
    rO = vlen(vhAB+vOh)
    return pO, rO

def tang2(pA, pB, rA, rB, t):
    AB = vlen(pB-pA)
    OB = AB*math.sin(t)
    OA = AB*math.cos(t)
    OE = OA-rA
    a = 4.*(OE**2-rB**2)
    b = -4.*OE*(OE**2-OB**2-rB**2)
    c = (OE**2-OB**2-rB**2)**2-4.*rB**2*OB**2
    delta = b**2-4.*a*c
    x = (-b+math.sqrt(delta))/(2.*a)
    r = OA-x-rA
    return r

def vangle(vA, vB):
    va = vlen(vA)
    vb = vlen(vB)
    c = dot(vA, vB)/(va*vb)
    return math.acos(c)

pA = np.array([0.,0.])
pB = np.array([1.,0.])
#aAD = 30
#aAD = 50
aAD = 40
aBD = 70
aCD = 180-aAD-aBD
aA = math.pi*aAD/180
aB = math.pi*aBD/180
aC = math.pi*aCD/180

rA = 0.3
rB = 0.2
rC = 0.15

xC = math.tan(aB)/(math.tan(aA)+math.tan(aB))
yC = xC*math.tan(aA)
pC = np.array([xC,yC])

AB = vlen(pB-pA)
BC = vlen(pC-pB)
CA = vlen(pA-pC)

pcrC = rA/(rA+rB)*pB+rB/(rA+rB)*pA
pcrA = rB/(rB+rC)*pC+rC/(rB+rC)*pB
pcrB = rC/(rC+rA)*pA+rA/(rC+rA)*pC
pO, rO = circumcircle(pcrA, pcrB, pcrC)
phA = vfoot(pA, pcrC, pcrB)
phB = vfoot(pB, pcrC, pcrA)
phC = vfoot(pC, pcrA, pcrB)

"""
ptA = pbetween(pA, phA, rA)
ptB = pbetween(pB, phB, rB)
ptC = pbetween(pC, phC, rC)

ptvA = pcrB-pcrC
ptvB = pcrC-pcrA
ptvC = pcrA-pcrB
ptcA = xpv2(ptB, -ptvB, ptC, ptvC)
ptcB = xpv2(ptA, ptvA, ptC, -ptvC)
ptcC = xpv2(ptA, -ptvA, ptB, ptvB)
ptchA = (ptcB+ptcC)/2
ptchB = (ptcC+ptcA)/2
ptchC = (ptcA+ptcB)/2
pOx, rOx = circumcircle(ptchA, ptchB, ptchC)
"""

vutAB = vunit(rotate2(pB-pA, (rA-rB)/AB))
ptABA = pA+rA*vutAB
ptABB = pB+rB*vutAB
vutBC = vunit(rotate2(pC-pB, (rB-rC)/BC))
ptBCB = pB+rB*vutBC
ptBCC = pC+rC*vutBC
vutCA = vunit(rotate2(pA-pC, (rC-rA)/CA))
ptCAC = pC+rC*vutCA
ptCAA = pA+rA*vutCA

ta = math.asin(rA/vlen(pcrC-pA))
vua = vunit(pA-pcrC)
vA = rotate(vua, -ta)
vB = rotate(-vua, ta)
vC = ptABA-ptABB
pxABA = xpv2(pcrC, vA, ptABB, vC)
pxABB = xpv2(pcrC, vB, ptABA, -vC)
pxABC = pcrC

ta = math.asin(rB/vlen(pcrA-pB))
vua = vunit(pB-pcrA)
vA = rotate(vua, -ta)
vB = rotate(-vua, ta)
vC = ptBCB-ptBCC
pxBCA = xpv2(pcrA, vA, ptBCC, vC)
pxBCB = xpv2(pcrA, vB, ptBCB, -vC)
pxBCC = pcrA

ta = math.asin(rC/vlen(pcrB-pC))
vua = vunit(pC-pcrB)
vA = rotate(vua, -ta)
vB = rotate(-vua, ta)
vC = ptCAC-ptCAA
pxCAA = xpv2(pcrB, vA, ptCAA, vC)pxCAB = xpv2(pcrB, vB, ptCAC, -vC)
pxCAC = pcrB

def t3ABC():
    def po(t, r):
        return pA+(r+rA)*rotate(vunit(pB-pA), t)

    def rdiff(t, r):
        pO = po(t, r)
        rc = vlen(pC-pO)-rC
        return r-rc

    tolerance = 1e-10

    def solve(t0, t1):
        t = (t0+t1)/2
        r = tang2(pA, pB, rA, rB, t)
        rd = rdiff(t, r)
        if math.fabs(rd) < tolerance: return t, r
        if rd > 0: return solve(t0, t)
        return solve(t, t1)

    t0 = -math.pi/4
    t1 = vangle(pB-pA, pC-pA)
    ts = np.linspace(t0, t1, 20)
    arange = []
    for t in ts:
        r = tang2(pA, pB, rA, rB, t)
        arange.append((t, r, rdiff(t, r)))
    tr0 = [x for x in arange if x[2]<0][-1]
    tr1 = [x for x in arange if x[2]>0][0]

    t0, t1 = tr0[0], tr1[0]
    t, r = solve(t0, t1)
    pO = po(t, r)
    return pO, r

#tx = 0.516
#rx = tang2(pA, pB, rA, rB, tx)
#pOABC = pA+(rx+rA)*rotate(vunit(pB-pA), tx)
#rOABC = rx
pOABC, rOABC = t3ABC()

ptxA = xpv2(pB, pC-pB, pA, pOABC-pA)
ptxB = xpv2(pC, pA-pC, pB, pOABC-pB)
ptxC = xpv2(pA, pB-pA, pC, pOABC-pC)
vA = vrotate(pOABC-pA)
vB = vrotate(pOABC-pB)
vC = vrotate(pOABC-pC)
ptxxA = xpv2(pB, -vB, pC, vC)
ptxxB = xpv2(pC, -vC, pA, vA)
ptxxC = xpv2(pA, -vA, pB, vB)
#lseg(pB, ptxxA, color="b")
#lseg(pC, ptxxA, color="b")
#polyline([ptxxA,ptxxB,ptxxC,ptxxA], color="b")

lseg(pA, pB, color="b")
lseg(pB, pC, color="b")
lseg(pC, pA, color="b")
#lseg(pcrC, pcrA, color="g")
#lseg(pcrA, pcrB, color="g")
#lseg(pcrB, pcrC, color="g")
lseg(pA, pcrA, color="g")
lseg(pB, pcrB, color="g")
lseg(pC, pcrC, color="g")
#lseg(pA, phA, color="k")
#lseg(pB, phB, color="k")
#lseg(pC, phC, color="k")
#lseg(ptA, ptB, color="k")
#lseg(ptB, ptC, color="k")
#lseg(ptC, ptA, color="k")
#lseg(ptcA, ptcB, color="g")
#lseg(ptcB, ptcC, color="g")
#lseg(ptcC, ptcA, color="g")
#lseg(ptchA, ptchB, color="g")
#lseg(ptchB, ptchC, color="g")
#lseg(ptchC, ptchA, color="g")
#lseg(ptABA, ptABB, color="g")
#lseg(ptBCB, ptBCC, color="g")
#lseg(ptCAC, ptCAA, color="g")
#lseg(pA, ptxA, color="b")
#lseg(pB, ptxB, color="b")
#lseg(pC, ptxC, color="b")
#lseg(pxABA, pxABB, color="b")
#lseg(pxABB, pxABC, color="b")
#lseg(pxABC, pxABA, color="b")
#lseg(pxBCA, pxBCB, color="b")
#lseg(pxBCB, pxBCC, color="b")
#lseg(pxBCC, pxBCA, color="b")
#lseg(pxCAA, pxCAB, color="b")
#lseg(pxCAB, pxCAC, color="b")
#lseg(pxCAC, pxCAA, color="b")


ax.add_patch(matplotlib.patches.Circle(pA, rA, color="r", fill=False))
ax.add_patch(matplotlib.patches.Circle(pB, rB, color="r", fill=False))
ax.add_patch(matplotlib.patches.Circle(pC, rC, color="r", fill=False))
#ax.add_patch(matplotlib.patches.Circle(pO, rO, color="r", fill=False))
#ax.add_patch(matplotlib.patches.Circle(pOx, rOx, color="r", fill=False))
ax.add_patch(matplotlib.patches.Circle(pcrA, .01, color="r", fill=False))
ax.add_patch(matplotlib.patches.Circle(pcrB, .01, color="r", fill=False))
ax.add_patch(matplotlib.patches.Circle(pcrC, .01, color="r", fill=False))
ax.add_patch(matplotlib.patches.Circle(pOABC, .01, color="r", fill=False))
ax.add_patch(matplotlib.patches.Circle(pOABC, rOABC, color="g", fill=False))

plt.show()