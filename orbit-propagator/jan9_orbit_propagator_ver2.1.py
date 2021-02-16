#######################
#Annalisa Cognigni
#Orbit Propagator Ver2
#UTAT Space Systems
#Orbit Propagator Project (ADCS)
#Dec 28th, 2020
#Version 2.1: Jan 9, 2020
#######################

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from convert_to_vectors2 import vectors

#Define constants##########################
radius=6378.0
mu=6.67408*(10**-20)*5.972*(10**24) #G*m
###########################################

#Define inputs#############################
eccentricity=0.7
semi_maj_axis=9000
inclination=0
raan=0
periapsis=3.49
t0=0
time_step=5
time_interval=20000
###########################################

def magnitude(a,b,c):
    r=(a)**2+(b)**2+(c)**2
    r=np.sqrt(r)
    return r



def orbit_propagator(t_int, mu, r0, v0):
    print(r0)
    print(v0)
    initial_conditions=np.concatenate((r0, v0),axis=None)

    print('######## INITIAL CONDITIONS ########')
    print(initial_conditions)


    t_list=[]
    for i in range(t0,t_int,time_step): #time step
        t_list=t_list+[i]

    print('###### TIME POINTS ######')
    print(t_list)

    #Equation:
    # x' = (-mu/|x|^3)*x

    def ode(t,x):  # x=[x v]
        #print('working...')
        q = -mu/(magnitude(x[0], x[1], x[2])**3)

        A = np.array([[0,0,0,1,0,0], [0,0,0,0,1,0], [0,0,0,0,0,1],[q, 0,0,0,0,0], [0,q,0,0,0,0],[0,0,q,0,0,0]])

        x_double_dot = np.matmul(A,x)

        return x_double_dot  # x'=[v, a]
    #x'=Ax

    data=solve_ivp(ode,[t0,t_int],initial_conditions,'DOP853',t_list) #integrate the equation defined above

    print('##### DATA #############')
    print(data)


    position=np.concatenate(([data.y[0]],[data.y[1]],[data.y[2]]),axis=0)

    print('| POSITION |')
    print(position)

    velocity=np.concatenate(([data.y[3]],[data.y[4]],[data.y[5]]),axis=0)

    print('| VELOCITY |')
    print(velocity)

    print('| TIME |')
    print(t_list)


    return position, velocity, t_list

#### References ######
# [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
# [2] https://scicomp.stackexchange.com/questions/34304/solving-coupled-differential-equations-in-python-2nd-order
# [3] H. D. Curtis, Orbital Mechanics for Engineering Students. Pg. e20.
# [4] https://numpy.org/doc/stable/reference/generated/numpy.concatenate.html

##########
#PLOT

def plots(position, velocity, time_step):
    x=position[0]
    y=position[1]
    z=position[2]

    u=velocity[0]
    v=velocity[1]
    w=velocity[2]

    testing_plot = plt.figure()
    l = len(x)

    ########EARTH
    rho = 6378  # radius of earth
    theta = np.linspace(0, 2 * np.pi, l)
    phi = np.linspace(0, np.pi, l)

    a = rho * np.outer(np.cos(theta), np.sin(phi))
    b = rho * np.outer(np.sin(theta), np.sin(phi))
    c = rho * np.outer(np.ones(l), np.cos(phi))

    earth = testing_plot.add_subplot(1, 2, 1, projection='3d')
    earth.plot_surface(a, b, c)

    #########PLOT DATA
    earth.scatter(x, y, z, c='red')
    earth.set_title('Satellite Position (Ver2)')
    zoom = 15000
    earth.set(xlim=(-zoom, zoom), ylim=(-zoom, zoom), zlim=(-zoom, zoom))

    velocity = testing_plot.add_subplot(1, 2, 2)
    speed = magnitude(u,v,w)
    velocity.plot(time_step, speed)
    velocity.set_title('Satellite Speed (Ver2)')
    velocity.set_xlabel('Time (s)')
    velocity.set_ylabel('Speed (km/s)')

    plt.show()
    return True

r0, v0=vectors(eccentricity, semi_maj_axis, inclination, raan, periapsis, t0, time_step)
print(r0,v0)

position, velocity, t_list = orbit_propagator(time_interval, mu, r0, v0)
plots(position, velocity, t_list)
