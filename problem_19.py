import matplotlib.pyplot as plt
import numpy as np
import math
import sympy

def run():
    a, b = 1., 3.
    x = np.linspace(.6062, 0.6063, 100)
    y = a/np.sin(x) + b/np.cos(x)
    x0 = math.atan((a/b)**(1/3))
    print(x0)
    plt.plot(x, y, color="r")
    plt.show()

def run2():
    a,b,x = sympy.symbols("a b x")
    y = a/sympy.sin(x) + b/sympy.cos(x)
    d1 = sympy.diff(y, x)
    print(d1)
    d1f = sympy.fraction(d1)
    print(d1f)

if __name__ == "__main__":
    run()
