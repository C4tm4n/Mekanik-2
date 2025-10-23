import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

L = 15
l = 5
M = 500
m = 50
mp = 5
a = 0.25
g = 9.82
rk = 2

def advanced(t, y ):
    phi, phi_prim, theta, theta_prim = y
    phi_biz =   (
                    (
                        (M*a*L + m*a**2*L/2)*g*math.cos(phi) - 
                        mp*L*math.cos(phi)* 
                        (
                    - (
                        L**2*phi_prim**2*(1-a) * L*math.sin(phi) + l*math.sin(phi+ theta)
                        )
                    /
                        (
                            ((1-a)*L*math.cos(phi) + l*math.cos(phi + theta))**2+
                            ((1-a)*L*math.sin(phi) +l*math.sin(phi+theta))**2
                        )
                        -l*theta_prim**2*math.sin(phi+theta)+g*math.sin(phi+theta)**2
                    )
                    +mp*L*math.sin(phi)*
                    (
                        (-L**2*phi**2 *((1-a) *L* math.cos(phi) + l*math.cos(phi+theta)))
                        /
                        (
                            ((1-a)*L*math.cos(phi) + l*math.cos(phi + theta))**2+
                            ((1-a)*L*math.sin(phi) +l*math.sin(phi+theta))**2
                        )
                        + g*math.cos(phi+theta)*math.sin(phi+theta) - l* theta_prim**2 * math.cos(phi+ theta)
                    )
                    )
                    /
                    (
                        a**3*m*L**2/3 + 1/3*(1-a)**3*m*L**2+M*a**2*L**2+
                        mp*(
                            (1-a)*L*math.cos(phi)+l*math.cos(phi+theta)**2 + ((1-a)*L*math.sin(phi) + l*math.sin(phi+theta))**2
                        )
                        + mp*L**2*(1-math.sin(phi)*math.cos(theta)* math.sin(phi+theta) - math.cos(phi)*math.cos(theta)*math.cos(phi+theta))
                    ) 
                )

    theta_biz = (-g*math.cos(phi+theta) - L*phi_biz*math.cos(theta))/l
    return [phi_prim, phi_biz, theta_prim, theta_biz]


if __name__ == "__main__":
    h = 4
    phi = -math.acos(h/((1-a)*L))
    theta = math.acos(h/((1-a)*L)) -math.pi
    y0 = [phi,  0.0, theta, 0.0]
    t = np.linspace(0,3, 50)
    sol = solve_ivp(advanced,[0,3], y0, t_eval=t)
    print(sol)



    positions = np.zeros((len(sol.t),2))
    for i in range(len(sol.t)):
        positions[i][0] = (1-a) * L* math.cos(sol.y[0][i])+ math.cos(sol.y[2][i]+ sol.y[0][i])* l
        positions[i][1] = (1-a) * L* math.sin(sol.y[0][i])+ math.sin(sol.y[2][i]+ sol.y[0][i])* l + h
    
    plt.subplot(2,3,1)
    plt.plot(positions[:,0], positions[:,1], label='Pendulum Path')
    plt.subplot(2,3,2)
    plt.plot(sol.t, sol.y[0], label='Angle $\phi$')
    plt.xlabel('Time [s]')
    plt.ylabel('angle $\phi$')

    plt.subplot(2,3,3)
    plt.plot(sol.t, sol.y[1], label='Angular Velocity $\phi$', color='orange')
    plt.xlabel('Time [s]')
    plt.ylabel('velocity $\phi$')
    plt.subplot(2,3,4)
    plt.plot(sol.t, sol.y[2], label='$\\theta$', color='blue')
    plt.xlabel('Time [s]')
    plt.ylabel('angle $\\theta$')
    plt.subplot(2,3,5)
    plt.xlabel('Time [s]')
    plt.ylabel('velocity $\\theta$')
    plt.plot(sol.t, sol.y[3], label='Angular Velocity $\\theta $', color='green')

#    plt.legend()

    plt.show()