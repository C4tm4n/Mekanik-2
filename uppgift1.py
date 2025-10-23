import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math


def simpleUnhinged(t, y ):
    u1, u2 = y
    return [u2, k*math.cos(u1)]


def energyLaw():
    startEnergy =L*g*math.sin(math.pi/2+startAngle) *( M+ m/2)
    positionalEnergy = g*(M*(h-a*L) + m*(L*(0.5-a)+h) )
    rotationEnergy = startEnergy-positionalEnergy
    angleSpeed = math.sqrt(rotationEnergy*2/inertia)
    print(f"radians /s: {angleSpeed}")
    print(f"deg/s: {angleSpeed/0.017453}")
    
    



if __name__ == "__main__":
    L = 15
    M = 500
    m = 50
    a = 0.25
    g = 9.82
    rk = 2
    inertia = L**2*(M*a**2+m/3*(a**3 + (1-a)**3)) # hinged I believe
  #  inertia = (L**2*(M*a**2+m/3*(a**3 + (1-a)**3))+M*rk**2) # unhinged 
    k = g*L*(a*M - (0.5-a)*m)/(inertia)
    print(f"k:{k}")
    h = 4
    startAngle = -math.acos(h/((1-a)*L))
    y0 = [startAngle,0.0]
    t = np.linspace(0,3, 50)
    energyLaw()
    sol = solve_ivp(simpleUnhinged,[0,3], y0, t_eval=t)
    
    for i in range(len(sol.t)):
        if(sol.t[i]>2.3):
            print(f"t: {sol.t[i]} deg: {sol.y[0][i]/(2*math.pi)*360} speed: {sol.y[1][i]/0.017453} ")
    plt.subplot(1,2,1)
    plt.plot(sol.t, sol.y[0], label='Angle $\phi$')
    plt.xlabel('Time [s]')

    plt.subplot(1,2,2)
    plt.plot(sol.t, sol.y[1], label='Angular Velocity $\omega$', color='orange')

    plt.legend()

    plt.show()


    


