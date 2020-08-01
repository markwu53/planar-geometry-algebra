import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import sympy
import numpy as np


fig, ax = plt.subplots()
plt.xlim(-1.5,2)
plt.ylim(-1.5,1.5)
ax.set_aspect(1)
#plt.grid(linestyle='--')

C1 = matplotlib.patches.Circle((0,0), .5, color="r", fill=False)
#ax.add_patch(C1)
#circle1 = plt.Circle((0, 0), 0.5, color='r')
#ax.add_artist(circle1)

aAD = 30
aBD = 70
aCD = 180-aAD-aBD
aA = math.pi*aAD/180
aB = math.pi*aBD/180
aC = math.pi*aCD/180
xA = 0
yA = 0
pA = (xA,yA)
xB = 1
yB = 0
pB = (xB,yB)
xC = math.tan(aB)/(math.tan(aA)+math.tan(aB))
yC = xC*math.tan(aA)
pC = (xC,yC)
x,y = list(zip(pA,pB))
plt.plot(x,y, color="b")
x,y = list(zip(pB,pC))
plt.plot(x,y, color="b")
x,y = list(zip(pC,pA))
plt.plot(x,y, color="b")

xI = math.tan(aB/2)/(math.tan(aA/2)+math.tan(aB/2))
yI = xI*math.tan(aA/2)
pI = (xI,yI)
x,y = list(zip(pA,pI))
#plt.plot(x,y, color="r")
x,y = list(zip(pB,pI))
#plt.plot(x,y, color="r")
x,y = list(zip(pC,pI))
#plt.plot(x,y, color="r")
rI = yI
cI = matplotlib.patches.Circle(pI, rI, color="r", fill=False)
ax.add_patch(cI)

xO = 1./2
yO = xO/math.tan(aC)
pO = (xO,yO)
x,y = list(zip(pA,pO))
#plt.plot(x,y, color="g")
x,y = list(zip(pB,pO))
#plt.plot(x,y, color="g")
x,y = list(zip(pC,pO))
#plt.plot(x,y, color="g")
rO = xO/math.sin(aC)
cO = matplotlib.patches.Circle(pO, rO, color="g", fill=False)
ax.add_patch(cO)

xHAB = .5*xA+.5*xB
yHAB = .5*yA+.5*yB
pHAB = (xHAB,yHAB)
xHBC = .5*xB+.5*xC
yHBC = .5*yB+.5*yC
pHBC = (xHBC,yHBC)
xHCA = .5*xC+.5*xA
yHCA = .5*yC+.5*yA
pHCA = (xHCA,yHCA)
x,y = list(zip(pHAB,pHBC))
#plt.plot(x,y, color="y")
x,y = list(zip(pHBC,pHCA))
#plt.plot(x,y, color="y")
x,y = list(zip(pHCA,pHAB))
#plt.plot(x,y, color="y")

xN = .25+.5*xC
yN = .5*yC-.25/math.tan(aC)
pN = (xN,yN)
rN = .25/math.sin(aC)
cN = matplotlib.patches.Circle(pN, rN, color="y", fill=False)
ax.add_patch(cN)

xHC = xC
yHC = 0
pHC = (xHC,yHC)
xHA = math.sin(aB)*math.sin(aB)
yHA = math.sin(aB)*math.cos(aB)
pHA = (xHA,yHA)
xHB = math.cos(aA)*math.cos(aA)
yHB = math.cos(aA)*math.sin(aA)
pHB = (xHB,yHB)
x,y = list(zip(pA,pHA))
#plt.plot(x,y, color="b")
x,y = list(zip(pB,pHB))
#plt.plot(x,y, color="b")
x,y = list(zip(pC,pHC))
#plt.plot(x,y, color="b")
x,y = list(zip(pHA,pHB))
#plt.plot(x,y, color="b")
x,y = list(zip(pHB,pHC))
#plt.plot(x,y, color="b")
x,y = list(zip(pHC,pHA))
#plt.plot(x,y, color="b")

def ccrot(p, theta):
    x,y = p
    nx = x*math.cos(theta)-y*math.sin(theta)
    ny = x*math.sin(theta)+y*math.cos(theta)
    return nx,ny

def trans(p, v):
    x,y = p
    v1,v2 = v
    return x+v1,y+v2

AB = 1
xOA = AB/(1.-math.tan(aA/2)*math.tan(aB/2))
yOA = xOA*math.tan(aA/2)
rOA = yOA
pOA = (xOA,yOA)
x,y = list(zip(pA,pOA))
#plt.plot(x,y, color="b")
x,y = list(zip(pB,pOA))
#plt.plot(x,y, color="b")
cOA = matplotlib.patches.Circle(pOA, rOA, color="b", fill=False)
ax.add_patch(cOA)

BC = math.sqrt((xB-xC)**2+(yB-yC)**2)
xOB = BC/(1.-math.tan(aB/2)*math.tan(aC/2))
yOB = xOB*math.tan(aB/2)
rOB = yOB
pOB = (xOB,yOB)
pOB = ccrot(pOB, math.pi-aB)
pOB = trans(pOB, (xB-xA,yB-yA))
x,y = list(zip(pB,pOB))
#plt.plot(x,y, color="b")
x,y = list(zip(pC,pOB))
#plt.plot(x,y, color="b")
cOB = matplotlib.patches.Circle(pOB, rOB, color="b", fill=False)
ax.add_patch(cOB)

CA = math.sqrt((xC-xA)**2+(yC-yA)**2)
xOC = CA/(1.-math.tan(aC/2)*math.tan(aA/2))
yOC = xOC*math.tan(aC/2)
rOC = yOC
pOC = (xOC,yOC)
pOC = ccrot(pOC, aA-math.pi)
pOC = trans(pOC, (xC-xA,yC-yA))
x,y = list(zip(pC,pOC))
#plt.plot(x,y, color="b")
x,y = list(zip(pA,pOC))
#plt.plot(x,y, color="b")
cOC = matplotlib.patches.Circle(pOC, rOC, color="b", fill=False)
ax.add_patch(cOC)

p2AC = pA
p2AC = trans(p2AC, (xA-xC,yA-yC))
p2AC = trans(p2AC, (xA-xC,yA-yC))
x,y = list(zip(pA,p2AC))
plt.plot(x,y, color="b")
p2CA = pC
p2CA = trans(p2CA, (xC-xA,yC-yA))
p2CA = trans(p2CA, (xC-xA,yC-yA))
x,y = list(zip(pC,p2CA))
plt.plot(x,y, color="b")

p2CB= pC
p2CB= trans(p2CB, (xC-xB,yC-yB))
p2CB= trans(p2CB, (xC-xB,yC-yB))
x,y = list(zip(pC,p2CB))
plt.plot(x,y, color="b")
p2BC = pB
p2BC = trans(p2BC, (xB-xC,yB-yC))
p2BC = trans(p2BC, (xB-xC,yB-yC))
x,y = list(zip(pB,p2BC))
plt.plot(x,y, color="b")

p2BA= pB
p2BA= trans(p2BA, (xB-xA,yB-yA))
p2BA= trans(p2BA, (xB-xA,yB-yA))
x,y = list(zip(pB,p2BA))
plt.plot(x,y, color="b")
p2AB = pA
p2AB = trans(p2AB, (xA-xB,yA-yB))
p2AB = trans(p2AB, (xA-xB,yA-yB))
x,y = list(zip(pA,p2AB))
plt.plot(x,y, color="b")

x,y = list(zip(pOA,pOB))
plt.plot(x,y, color="r")
x,y = list(zip(pOB,pOC))
plt.plot(x,y, color="r")
x,y = list(zip(pOC,pOA))
plt.plot(x,y, color="r")

plt.show()