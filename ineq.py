import math
import matplotlib
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as np
import sympy

def run():
    x = np.linspace(0.1,1,100)
    x = x*np.pi/2
    y = np.sin(x)
    y = (y**3+y**2-1)/(2*y**4)
    a = -np.sin(x)**2*np.tan(x)**2
    b = np.sin(x)*np.tan(x)**2-np.cos(x)
    c = 1
    d = b**2-4*a*c
    u = -b/(2*a)
    v = np.sqrt(d)/(2*a)
    #plt.plot(x,u+v, color="b")
    plt.plot(x,(u-v)**2, color="r")
    plt.show()

if __name__ == "__main__":
    run()