#######################
#Annalisa Cognigni
#Simple Harmonic Oscillator
#UTAT Space Systems
#Orbit Propagator Project (ADCS)
#Dec 3, 2020
#######################

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


def harm_osc(t_int, w, x0, v0):
    initial_conditions=[x0,v0]
    print('######## INITIAL CONDITIONS ########')
    print(initial_conditions)
    print('####################################')

    t_list=[]
    for i in range(0,t_int,1):
        t_list=t_list+[i*0.01]

    print('###### TIME POINTS ######')
    print(t_list)
    print('#########################')

    #Equation:
    # x''+(w^2)x=0
    # x'' = -(w^2)x

    A=np.array([[0, 1],[-1*(w**2), 0]])
    print(A)

    def ode(t,x):  # x=[x, v]
        x_double_dot = np.matmul(A,x)
        return x_double_dot  # x'=[v, a]
    #x'=Ax

    data=solve_ivp(ode,[0,t_int],initial_conditions,'RK45',t_list) #integrate the equation defined above

    print(data.y[0])
    print(data.t)

    #plot position
    plt.plot(data.t,data.y[0])
    plt.xlabel('time (s)')
    plt.ylabel('position (m)')
    plt.show()

    #plot velocity
    plt.plot(data.t,data.y[1])
    plt.xlabel('time (s)')
    plt.ylabel('velocity (m/s)')
    plt.show()

    return data

# harm_osc(t_int, w, x0, v0)
#harm_osc(60, 20, 0, 2)
#harm_osc(60, 40, 0, 0)
harm_osc(60, 40, 10, 0)

#### References ######
# [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
# [2] https://scicomp.stackexchange.com/questions/34304/solving-coupled-differential-equations-in-python-2nd-order



