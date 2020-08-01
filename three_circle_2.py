import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import sympy
import numpy as np

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

def rad(degree):
    return math.pi*degree/180

def vlen(v):
    x,y = v
    return math.sqrt(x**2+y**2)

def dot(vA, vB):
    xA,yA = vA
    xB,yB = vB
    return xA*xB+yA*yB

def rotate(v, a):
    x,y = v
    c = math.cos(a)
    s = math.sin(a)
    return np.array([x*c-y*s, x*s+y*c])

def urotate(v, u):
    x,y = v
    c,s = u
    return np.array([x*c-y*s,x*s+y*c])

def vrotate(v):
    x,y = v
    return np.array([-y,x])

def vunit(v):
    return v/vlen(v)

def vangle(vA, vB):
    va = vlen(vA)
    vb = vlen(vB)
    c = dot(vA, vB)/(va*vb)
    return math.acos(c)

def vcos(vA, vB):
    va = vlen(vA)
    vb = vlen(vB)
    c = dot(vA, vB)/(va*vb)
    return c

def vfoot(pX, pA, pB):
    AB = vlen(pB-pA)
    BC = vlen(pX-pB)
    CA = vlen(pA-pX)
    x = (CA**2-BC**2)/AB**2
    qA = .5*(1+x)
    qB = .5*(1-x)
    return qB*pA+qA*pB

def vfoot3(pA, pB, pC):
    pvA = vfoot(pA, pB, pC)
    pvB = vfoot(pB, pC, pA)
    pvC = vfoot(pC, pA, pB)
    return pvA, pvB, pvC

#given A,B and length of AC,BC find C
def ABfindCby2e(pA, pB, AC, BC):
    vAB = pB-pA
    AB = vlen(vAB)
    x = (AC**2+AB**2-BC**2)/(2*AB)
    y = math.sqrt(AC**2-x**2)
    vAC = urotate(np.array([x,y]), vunit(vAB))
    pC = pA+vAC
    return pC

#given A,B and angle of A,B find C
def ABfindCby2ang(pA, pB, aA, aB):
    x = math.tan(aB)/(math.tan(aA)+math.tan(aB))
    y = x*math.tan(aA)
    u = np.array([x,y])
    vAB = pB-pA
    AB = vlen(vAB)
    vAC = urotate(AB*u, vunit(vAB))
    pC = pA+vAC
    return pC

def p3circle(pA, pB, pC):
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

def xpv2(pA, vA, pB, vB):
    xA,yA = pA
    xB,yB = pB
    xvA,yvA = vA
    xvB,yvB = vB
    if xvA == 0 and xvB == 0: return
    if xvA ==0:
        mB = yvB/xvB
        x = xA
        y = yB+mB*(x-xB)
        return np.array([x,y])
    if xvB == 0:
        mA = yvA/xvA
        x = xB
        y = yA+mA*(x-xA)
        return np.array([x,y])
    mA = yvA/xvA
    mB = yvB/xvB
    if mA == mB: return
    x = (yB-yA+mA*xA-mB*xB)/(mA-mB)
    y = yA+mA*(x-xA)
    return np.array([x,y])

def p3_in_circle(pA, pB, pC):
    vAB = pB-pA
    vAC = pC-pA
    vBC = pC-pB
    vBA = pA-pB
    aA = vangle(vAB, vAC)
    aB = vangle(vBC, vBA)
    pI = ABfindCby2ang(pA, pB, aA/2, aB/2)
    rI = vlen(pI-pA)*math.sin(aA/2)
    return pI, rI

#different triangles may have the same vfoot3,
#so vfoot3_inv may not be "inverse" to vfoot3
def vfoot3_inv(pA, pB, pC):
    pI, rI = p3_in_circle(pA, pB, pC)
    vA = vrotate(pI-pA)
    vB = vrotate(pI-pB)
    vC = vrotate(pI-pC)
    pxA = xpv2(pB, vB, pC, vC)
    pxB = xpv2(pC, vC, pA, vA)
    pxC = xpv2(pA, vA, pB, vB)
    return pxA, pxB, pxC

def half_triangle(pA, pB, pC):
    phA = (pB+pC)/2
    phB = (pC+pA)/2
    phC = (pA+pB)/2
    return phA, phB, phC

def half_triangle_inv(pA, pB, pC):
    phiA = pB+pC-pA
    phiB = pC+pA-pB
    phiC = pA+pB-pC
    return phiA, phiB, phiC

def circle_tangent_3_circles():
    pA = np.array([0.,0.])
    pB = np.array([1.,0.])
    aA = rad(40)
    aB = rad(35)
    pC = ABfindCby2ang(pA, pB, aA, aB)
    rA = 0.3
    rB = 0.2
    rC = 0.25
    def solve(r0, r1):
        if r1-r0 < tolerance:
            return r0
        r = (r0+r1)/2
        AT = rA+r
        BT = rB+r
        pT = ABfindCby2e(pA, pB, AT, BT)
        r2C = vlen(pT-pC)-rC
        if r2C > r: return solve(r, r1)
        return solve(r0, r)
    def check(r):
        AT = rA+r
        BT = rB+r
        pT = ABfindCby2e(pA, pB, AT, BT)
        r2C = vlen(pT-pC)-rC
        #print(r, r2C)
        return r, r2C
    def checkr():
        AB = vlen(pB-pA)
        r = (AB-(rA+rB))/2
        return check(r)
    r, r2C = checkr()
    if r > r2C:
        pA,pB,pC = pB,pC,pA
        rA,rB,rC = rB,rC,rA
        r, r2C = checkr()
        if r > r2C:
            pA,pB,pC = pB,pC,pA
            rA,rB,rC = rB,rC,rA
            r, r2C = checkr()
    tolerance = 1e-10
    r0, r1 = r, r2C
    while math.fabs(r0-r1) > tolerance:
        r = (r0+r1)/2
        r, r2C = check(r)
        if r < r2C:
            r0, r1 = r, r1
        else:
            r0, r1 = r0, r
    #r = solve(r, r2C)
    r = r0
    AT = rA+r
    BT = rB+r
    pT = ABfindCby2e(pA, pB, AT, BT)
    print(r)
    polygon([pA,pB,pC], color="b")
    ax.add_patch(matplotlib.patches.Circle(pA, rA, color="r", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pB, rB, color="r", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pC, rC, color="r", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pT, r, color="g", fill=False))

def test():
    pA = np.array([0.,0.])
    pB = np.array([1.,0.])
    aA = rad(60)
    aB = rad(35)
    pC = ABfindCby2ang(pA, pB, aA, aB)
    pO, rO = p3circle(pA, pB, pC)
    pI, rI = p3_in_circle(pA, pB, pC)
    phA,phB,phC = half_triangle(pA, pB, pC)
    phiA,phiB,phiC = half_triangle_inv(pA, pB, pC)
    pxA,pxB,pxC = vfoot3_inv(pA, pB, pC)
    pOI,rOI = p3_in_circle(pxA, pxB, pxC)

    polygon([pA,pB,pC], color="b")
    #polygon([phA,phB,phC], color="r")
    #polygon([phiA,phiB,phiC], color="g")
    polygon([pxA,pxB,pxC], color="g")
    ax.add_patch(matplotlib.patches.Circle(pO, rO, color="g", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pI, rI, color="r", fill=False))
    ax.add_patch(matplotlib.patches.Circle(pOI, rOI, color="r", fill=False))

    plt.show()

def test2():
    pA = np.array([0.,0.])
    pB = np.array([1.,0.])
    aA = rad(60)
    aB = rad(35)
    pC = ABfindCby2ang(pA, pB, aA, aB)
    yP = -.2
    xP,_ = (pA+pB)/2
    pP = np.array([xP,yP])
    xQ = xP
    theta = vangle(pP-pA, pC-pA)
    vuBQ = rotate(pC-pB,math.pi-theta)
    pQ = xpv2(pP, vrotate(pB-pA), pB, vuBQ)
    aP = vangle(pA-pC, pP-pC)
    aQ = vangle(pB-pC, pQ-pC)
    b,c,t = math.tan(aA),math.tan(aB),math.tan(theta)
    print(math.tan(aP))
    print(t*(b+c)/(c+2*t*b*c-b))
    print(aP, aQ)
    polyline([pA,pP,pC], color="g")
    polyline([pC,pQ,pB], color="g")
    polygon([pA,pB,pC], color="r")
    plt.show()

#export
def line_circle(pA, pB, pO, rO):
    xA,yA = pA
    xB,yB = pB
    xO,yO = pO
    if xA == xB:
        x1,x2 = xA,xA
        d = math.sqrt(rO**2-(xO-xA)**2)
        y1,y2 = yO+d,yO-d
    else:
        m = (yB-yA)/(xB-xA)
        t = yA-m*xA-yO
        a = m**2+1
        b = 2*m*t-2*xO
        c = t**2+xO**2-rO**2
        d = b**2-4*a*c
        x1 = (-b+math.sqrt(d))/(2*a)
        x2 = (-b-math.sqrt(d))/(2*a)
        y1 = yA+m*(x1-xA)
        y2 = yA+m*(x2-xA)
    return np.array([x1,y1]), np.array([x2,y2])

def test3():
    pB = np.array([0.,0.])
    pC = np.array([1.,0.])
    #aB = rad(55)
    #aC = rad(75)
    aB = rad(45)
    aC = rad(75)
    pA = ABfindCby2ang(pB, pC, aB, aC)
    pD = (pA+pB)/2
    pE = (pA+pC)/2
    pM = (pB+pC)/2
    pO,rO = p3circle(pA, pB, pC)
    pH,rI = p3_in_circle(*vfoot3(pA, pB, pC))
    pF,pF1 = line_circle(pD, pH, pO, rO)
    pG1,pG = line_circle(pE, pH, pO, rO)
    pN,pN1 = line_circle(pM, pH, pO, rO)
    pI = xpv2(pD, pE-pD, pG, pF-pG)
    pX,rX = p3circle(pD, pF, pE)
    #print(vlen(pF1-pG1)-vlen(pC-pB))
    polygon([pA,pB,pC], color="k")
    lseg(pD, pI)
    lseg(pG, pI)
    #lseg(pD, pF, color="b")
    #lseg(pE, pG, color="b")
    lseg(pF1, pF, color="b")
    lseg(pG1, pG, color="b")
    lseg(pG1, pF1, color="b")
    #lseg(pN, pN1, color="b")
    #lseg(pO, pA, color="b")
    lseg(pI, pA, color="b")
    #lseg(pO, pD, color="b")
    #lseg(pO, pE, color="b")
    ax.add_patch(matplotlib.patches.Circle(pO, rO, color="g", fill=False))
    #ax.add_patch(matplotlib.patches.Circle(pX, rX, color="r", fill=False))
    plt.show()

def test4():
    pA = np.array([0,0])
    pB = np.array([1,0])
    aA = rad(50)
    aB = rad(50)
    pC = ABfindCby2ang(pA, pB, aA, aB)
    pD = ABfindCby2ang(pA, pB, rad(40), rad(30))
    a = vangle(pC-pD, pA-pD)
    print(a/math.pi*180)

if __name__ == "__main__":
    test4()
