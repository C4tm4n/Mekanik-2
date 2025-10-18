import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
import math


def simpleUnhinged(t, y ):
    u1, u2 = y
    return [u2, k*math.cos(u1)]



    
    



if __name__ == "__main__":
    L = 15
    M = 500
    m = 50
    a = 0.25
    g = 9.82
    rk = 2
    k = g*L*(a*M + (0.5-a)*m)/(L**2*(M*a+m/3*(a**3 + (1-a)**3))+M*rk**2)
    K   =   g*(M + (m**2)/2 - (m*(1-a)**2)/2)   /   (L*(M**2 + (m**3)/3 + (m*(1-a)**3)/3))
    print(K)
    h = 4
    y0 = [math.cos(h/(1-a)*L),0.0]
    t = np.linspace(0,3, 50)
    sol = solve_ivp(simpleUnhinged,[0,3], y0, t_eval=t)
    print(k)

    plt.subplot(1,2,1)
    plt.plot(sol.t, sol.y[0], label='Angle $\phi$')
    plt.xlabel('Time [s]')

    plt.subplot(1,2,2)
    plt.plot(sol.t, sol.y[1], label='Angular Velocity $\omega$', color='orange')

    plt.legend()

    plt.show()

    


