import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import sympy
import numpy as np

x,y,r,x1,y1,r1,x2,y2,r2,x3,y3,r3 = sympy.symbols("x y r x1 y1 r1 x2 y2 r2 x3 y3 r3")
A13 = x1**2-x3**2+y1**2-y3**2-(r1-r3)**2
A = 4*(x1-x3)**2-4*(r1-r3)**2
B = 8*(x1-x3)*(y1-y3)
C = 4*(y1-y3)**2-4*(r1-r3)**2
D = -4*(x1-x3)*A13+8*(r1-r3)**2*x3
E = -4*(y1-y3)*A13+8*(r1-r3)**2*y3
F = A13**2-4*(r1-r3)**2*x3**2-4*(r1-r3)**2*y3**2

print(sympy.expand(A))
print(sympy.expand(B))
print(sympy.expand(C))
print(sympy.expand(D))
print(sympy.expand(E))
print(sympy.expand(F))