import math
pi = math.pi
a = 6378
import matplotlib.pyplot as plt

def Ke(dn_dh):
    return 1/(1+a*dn_dh)
def eq2_28a(dn_dh, theta, s):
    ke = Ke(dn_dh)
    h = ke*a*(math.cos(theta)/math.cos(theta+s/(ke*a))-1)
    return h
def eq2_30(dn_dh, theta, s): # beta0 == dn_dh
    h = ((((s+s*dn_dh*a)/math.cos(theta))+a*math.sin(theta))**2 - a**2 * (math.sin(theta)**2)) / (2*a*(1+dn_dh*a))
    return h

def q2():
    dn_dh = -(1/(4*a))
    x = [s for s in range(200+1)]
    delta_list = []
    for theta in [0.0*2*pi/360, 0.5*2*pi/360, 1.0*2*pi/360, 3.0*2*pi/360]:
        y1, y2, delta = [], [], []
        for s in range(200+1):
            y1.append(eq2_28a(dn_dh, theta, s))
            y2.append(eq2_30(dn_dh, theta, s))
            delta.append(y1[-1]-y2[-1])
        delta_list.append(delta)
        plt.plot(x, y1, label='2_28a {}'.format(str(theta*360/2/pi)))
        plt.plot(x, y2, label='2_30 {}'.format(str(theta*360/2/pi)))
    plt.legend()
    plt.title('Propagation trajectory of two formulas at different elevation angles')
    plt.xlabel('s(km)')
    plt.ylabel('h(km)')
    plt.show()
    plt.close()

    theta = [0.0, 0.5, 1.0, 3.0]
    for i in range(4):
        plt.plot(x, delta_list[i], label=str(theta[i]))
    plt.legend()
    plt.title('The difference in the propagation trajectory of the two formulas at different elevation angles')
    plt.xlabel('s(km)')
    plt.ylabel('delta h(km)')
    plt.show()
    plt.close()

def q3():
    dn_dh = -(3*10**(-4))
    x = [s for s in range(50+1)]
    ground = [0 for _ in range(50+1)]
    plt.plot(x, ground, label='ground')
    delta_list = []
    for theta in [0.1*2*pi/360, 0.2*2*pi/360, 0.3*2*pi/360]:
        y1, y2, delta = [], [], []
        for s in range(50+1):
            y1.append(eq2_28a(dn_dh, theta, s))
            y2.append(eq2_30(dn_dh, theta, s))
            delta.append(y1[-1]-y2[-1])
        delta_list.append(delta)
        plt.plot(x, y1, label='2_28a {}'.format(str(theta*360/2/pi)))
        plt.plot(x, y2, label='2_30 {}'.format(str(theta*360/2/pi)))
    plt.legend()
    plt.title('Propagation trajectory of two formulas at different elevation angles')
    plt.xlabel('s(km)')
    plt.ylabel('h(km)')
    plt.show()
    plt.close()

    theta = [0.1, 0.2, 0.3]
    ground = [0 for _ in range(50+1)]
    plt.plot(x, ground, label='ground')
    for i in range(len(theta)):
        plt.plot(x, delta_list[i], label=str(theta[i]))
    plt.legend()
    plt.title('The difference in the propagation trajectory of the two formulas at different elevation angles')
    plt.xlabel('s(km)')
    plt.ylabel('delta h(km)')
    plt.show()
    plt.close()

def q5():
    #import scipy
    from scipy.special import jv
    import numpy as np
    l = 0.1 # lambda, m
    d = 11.89 # antenna diameter, m
    
    theta = np.linspace(-3*2*pi/360, 3*2*pi/360, 100)
    f = []
    for t in theta:
        item = pi*d*math.sin(t)/l
        f.append((8*jv(2, item)/(item)**2)**2)
    plt.plot(theta, f, label='f^2 -- theta map')
    plt.title('f^2 -- theta map')
    plt.xlabel('theta')
    plt.ylabel('f^2(theta)')
    plt.show()
    plt.close()

    plt.plot(theta, f, label='f^2 -- theta map')
    plt.title('f^2 -- theta map')
    plt.xlabel('theta')
    plt.ylabel('f^2(theta)')
    plt.xlim(0.01, 0.02)
    plt.ylim(0.0, 0.1)
    plt.show()
    plt.close()


if __name__ == "__main__":
    #q2()
    q3()
    #q5()