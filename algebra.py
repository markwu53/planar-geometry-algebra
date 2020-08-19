import sympy

def run():
    x,y,z = sympy.symbols("x y z")
    ex = x**4+y**4+z**4-x*y*z
    ex2 = sympy.expand(ex.subs({z:1-x-y}))
    #print(ex2)
    homo_ex2 = 2*x**4 + 4*x**3*y - 4*x**3*z + 6*x**2*y**2 - 11*x**2*y*z + 6*x**2*z**2 + 4*x*y**3 - 11*x*y**2*z + 11*x*y*z**2 - 4*x*z**3 + 2*y**4 - 4*y**3*z + 6*y**2*z**2 - 4*y*z**3 + z**4
    homo_ex2_a = 2*x**4 + 4*x**3*y + 6*x**2*y**2 + 6*x**2*z**2 + 4*x*y**3 + 11*x*y*z**2 + 2*y**4 + 6*y**2*z**2 + z**4
    homo_ex2_b = 4*x**3*z + 11*x**2*y*z + 11*x*y**2*z + 4*x*z**3 + 4*y**3*z + 4*y*z**3
    print(homo_ex2_a)
    print(homo_ex2_b)
    homo_ex2_1 = 2*x**4 + 4*x**3*y + 6*x**2*y**2 + 6*x**2 + 4*x*y**3 + 11*x*y + 2*y**4 + 6*y**2 + 1
    homo_ex2_2 = + 4*x**3 + 11*x**2*y + 11*x*y**2 + 4*x + 4*y**3 + 4*y
    print(homo_ex2_1)
    print(homo_ex2_2)

if __name__ == "__main__":
    run()
