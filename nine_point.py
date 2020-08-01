from sympy import *
import math

def test1():
    x, y = symbols('x y')
    xc = tan(y)/(tan(x)+tan(y))
    yc = xc*tan(y)
    xn = xc/2+1/4
    tz = yc/(sin(x)*sin(x))-cot(x)
    yn = yc/2-tz/4
    xo = tan(y/2)/(tan(y/2)-tan(x/2))
    yo = xo*tan(x/2)
    ro = yo
    rn = sqrt(1+tz*tz)/4
    d2 = (xo-xn)**2+(yo-yn)**2
    r2 = (ro+rn)**2

    result = diff(d2-r2, x)
    print(result)

def test2():
    x, y = symbols('x y')

    sinHa=x
    cosHa=sqrt(1-sinHa**2)
    tanHa=sinHa/cosHa
    sina=2*sinHa*cosHa
    cosa=sqrt(1-sina**2)
    tana=sina/cosa

    sinHb=y
    cosHb=sqrt(1-sinHb**2)
    tanHb=sinHb/cosHb
    sinb=2*sinHb*cosHb
    cosb=sqrt(1-sinb**2)
    tanb=sinb/cosb

    xc=tanb/(tana+tanb)
    yc=tana*tanb/(tana+tanb)
    xo=1/(1-tanHa*tanHb)
    yo=xo*tanHa
    ro=yo
    xn=1/4+xc/2
    tann=(xc-cosa**2)/(sina*cosa)
    yn=yc/2-tann/4
    rn=sqrt(1+tann**2)/4
    don2=(xo-xn)**2+(yo-yn)**2
    rorn2=(ro+rn)**2
    print(don2)
    print(rorn2)

    subs={x:math.sin(math.pi*35/180), y:math.sin(math.pi*30/180)}
    print(don2.evalf(subs=subs))
    print(rorn2.evalf(subs=subs))

def test3():
    x, y= symbols("x y")
    ex=(x-1)**2
    fx=x**2-2*x+1
    print(ex)
    print(expand(ex)==fx)

def test4():
    x, y= symbols("x y")
    ynu=4*x*y+(1-x**2)*(1-y**2)
    xnu=2+4*y*(1-x**2)
    nl=8*x*(1-y**2)+8*y*(1-x**2)
    rnl=4*(1+x**2)*(1+y**2)
    ol=1-x*y
    lhs=rnl**2*(xnu*ol-nl)**2+rnl**2*(ynu*ol-x*nl)**2
    rhs=nl**2*(ol+x*rnl)**2
    #print(lhs)
    #print(expand(lhs))
    #print(rhs)
    #print(expand(rhs))
    #print(expand(lhs)==expand(rhs))

    subs={x:math.tan(math.pi*35/180), y:math.tan(math.pi*30/180)}
    print(lhs.evalf(subs=subs))
    print(rhs.evalf(subs=subs))

def test5():
    x, y= symbols("x y")
    xc1=y*(1-x**2)
    xc2=x*(1-y**2)+y*(1-x**2)
    xc=xc1/xc2
    yc1=2*x*y
    yc2=xc2
    yc=yc1/yc2

    tn1=y*(1+x**2)**2-(1-x**2)*(x*(1-y**2)+y*(1-x**2))
    tn2=2*x**2*(1-y**2)+2*x*y*(1-x**2)
    #tann=xc*(1+x**2)**2/2/x/(1-x**2)-(1-x**2)/2/x
    tn=tn1/tn2
    #xn=1/4+xc/2
    xn1=xc2+2*xc1
    xn2=4*xc2
    xn=xn1/xn2
    #yn=yc/2-tn/4
    yn1=2*yc1*tn2-yc2*tn1
    yn2=4*yc2*tn2
    yn=yn1/yn2
    #rn=sqrt(1+tn**2)/4
    rn1=x*(1+x**2)*(1+y**2)
    rn2=4*tn2
    rn=rn1/rn2
    #lhs=(xo-xn)**2+(yo-yn)**2
    #rhs=(ro+rn)**2

    def in_circle():
        xo1=y
        xo2=x+y
        xo=xo1/xo2
        yo1=x*y
        yo2=xo2
        yo=yo1/yo2
        ro=yo

        lhs1=(xo1*xn2-xo2*xn1)**2*(yo2*yn2)**2
        lhs2=(yo1*yn2-yo2*yn1)**2*(xo2*xn2)**2
        lhs3=(yo2*rn2)**2
        lhs=(lhs1+lhs2)*lhs3
        rhs1=(yo2*rn1-yo1*rn2)**2
        rhs2=(xo2*xn2*yo2*yn2)**2
        rhs=rhs1*rhs2

        subs={x:math.tan(math.pi*35/180), y:math.tan(math.pi*30/180)}
        print(lhs.evalf(subs=subs))
        print(rhs.evalf(subs=subs))
        print(expand(lhs)==expand(rhs))
        #print(expand(lhs))

    def out_circle():
        xo1=1
        xo2=1-x*y
        xo=xo1/xo2
        yo1=x
        yo2=xo2
        yo=yo1/yo2
        ro=yo

        lhs1=(xo1*xn2-xo2*xn1)**2*(yo2*yn2)**2
        lhs2=(yo1*yn2-yo2*yn1)**2*(xo2*xn2)**2
        lhs3=(yo2*rn2)**2
        lhs=(lhs1+lhs2)*lhs3
        rhs1=(yo1*rn2+yo2*rn1)**2
        rhs2=(xo2*xn2*yo2*yn2)**2
        rhs=rhs1*rhs2

        subs={x:math.tan(math.pi*35/180), y:math.tan(math.pi*30/180)}
        print(lhs.evalf(subs=subs))
        print(rhs.evalf(subs=subs))
        print(expand(lhs)==expand(rhs))
        #print(expand(lhs))

    in_circle()
    out_circle()

test5()