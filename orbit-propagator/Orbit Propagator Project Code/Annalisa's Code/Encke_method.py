#######################
#Annalisa Cognigni
#Encke's Method Code
#UTAT Space Systems
#Orbit Propagator Project (ADCS)
#######################


import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from convert_to_vectors2_copy import vectors

########################################
#define orbital parameters
eccentricity=0.007      #0.007
semi_maj_axis=6920    #6920
inclination=1.7
raan=4
periapsis=0
t0=0

print('############### initial elements #################')
print([eccentricity, semi_maj_axis, inclination, raan, periapsis, t0])
#########################################

#########################################
#Initial conditions
dt=1
error=0.0000000001
iterations= 60*90       #2592000 #1 month
r0, v0 =vectors(eccentricity, semi_maj_axis, inclination, raan, periapsis, t0, dt)
print('################### INITIAL CONDITIONS #################')
print(r0)
print(v0)
#orbit=Orbit(eccentricity, semi_maj_axis, inclination, raan, periapsis, t0, time_step, duration)


#Define constants
J2=0.00108263
radius=6378.0
mu=6.67408*(10**-20)*5.972*(10**24) #G*m
#########################################


#Calculate magnitude of vector
def magnitude(vector):
    r=(vector[0])**2+(vector[1])**2+(vector[2])**2
    r=np.sqrt(r)
    return r

#Stumpff Coefficients (pg. 174 Curtis)
def S(z):
    if z>0:
        r=(np.sqrt(z)-np.sin(np.sqrt(z)))/(np.sqrt(z)**3)
    elif z<0:
        r=(np.sinh(np.sqrt(-1*z))-np.sqrt(-1*z))/(np.sqrt(-1*z)**3)
    elif z==0:
        r=1/6
    else:
        print('ERROR')
        print(z)
        return 'error'
    return r

def C(z):
    if z>0:
        r=(1-np.cos(np.sqrt(z)))/z
    elif z<0:
        r=(np.cosh(np.sqrt(-1*z))-1)/(-1*z)
    elif z==0:
        r=1/2
    else:
        print('ERROR')
        print(z)
        return 'error'
    return r
#############

def get_anomaly(magnitude_r0,v_r0,mu,alpha,error,time_passed):
    X_i = np.sqrt(mu) * np.abs(alpha) * time_passed
    #print('anomaly:')
    #print(X_i)
    while True:

        f_Xi = ((magnitude_r0 * v_r0) / np.sqrt(mu)) * (X_i ** 2) * C(alpha * (X_i ** 2)) + (
                    1 - alpha * magnitude_r0) * (X_i ** 3) * S(alpha * (X_i ** 2)) + magnitude_r0 * X_i - np.sqrt(
            mu) * time_passed

        f_prime_Xi = ((magnitude_r0 * v_r0) / np.sqrt(mu)) * (X_i) * (
                    1 - alpha * (X_i ** 2) * S(alpha * (X_i ** 2))) + (1 - alpha * magnitude_r0) * X_i ** 2 * C(
            alpha * (X_i ** 2)) + magnitude_r0

        ratio = np.abs(f_Xi / f_prime_Xi)
        if (ratio > error):
            X_i = X_i - f_Xi / f_prime_Xi
        else:
            return X_i


def step_forward_orbit(r0,v0,time_passed,mu,error):

    v_r0=np.dot(r0,v0)/magnitude(r0)
    alpha=(2/magnitude(r0))-(magnitude(v0)**2)/mu

    #Algorithm 3.3 (pg.178)
    X=get_anomaly(magnitude(r0), v_r0, mu, alpha, error, time_passed)

    #Lagrange Coefficients (Algorithm 3.4)
    f=1-(X**2)/magnitude(r0)*C(alpha*(X**2))
    g=time_passed-(1/np.sqrt(mu))*(X**3)*S(alpha*(X**2))

    #r_osc
    r1=f*r0+g*v0
    magnitude_r1=np.sqrt(np.dot(r1,r1))
    #####################################

    f_dot=((np.sqrt(mu))/(magnitude_r1*magnitude(r0)))*(alpha*(X**3)*S(alpha*(X**2))-X)
    g_dot=1-((X**2)/magnitude_r1)*C(alpha*(X**2))

    #v_osc
    v1=f_dot*r0+g_dot*v0

    return r1, v1

#Reduce numerical computation error:
def f_expansion(q):
    r=q*(3+3*q+q**2)/(((1+q)**(3/2))+((1+q)**3))
    return r



def encke(r0, v0, dt, mu, error, iterations):
    data=np.zeros((iterations,6))
    ic_r0 = r0
    ic_v0 = v0

    #initial conditions
    dr=np.array([0, 0, 0]) #delta

    #r0 and v0 are initial position and velocity vectors
    time=np.zeros((iterations,1))

    data[0, 0] = r0[0]
    data[0, 1] = r0[1]
    data[0, 2] = r0[2]
    data[0, 3] = v0[0]
    data[0, 4] = v0[1]
    data[0, 5] = v0[2]
    time[0, 0] = 0

    r_i=r0
    r_osc_i=r0


    for i in range(1,iterations):
        time_passed=i*dt

        ######################################################
        # define J2 perturbation/acceleration vector: (pg. 664, Curtis)
        accel = np.zeros((1, 3))
        constant = (3 / 2) * (J2 * mu * (radius ** 2)) / (magnitude(r_i) ** 4)

        x = r_i[0]
        y = r_i[1]
        z = r_i[2]

        accel[0, 0] = constant * (x / magnitude(r_i)) * (5 * (z ** 2) / (magnitude(r_i) ** 2) - 1)
        accel[0, 1] = constant * (y / magnitude(r_i)) * (5 * (z ** 2) / (magnitude(r_i) ** 2) - 1)
        accel[0, 2] = constant * (z / magnitude(r_i)) * (5 * (z ** 2) / (magnitude(r_i) ** 2) - 3)
        # print('acceleration:')
        # print(accel)

        # calculate new dr, dv:
        # define differential equation
        q = (np.dot((dr + 2 * r_osc_i), dr)) / (magnitude(r_osc_i) ** 2)
        a = mu / (magnitude(r_osc_i) ** 3)
        b = a * f_expansion(q) * r_i + accel

        # differential eq: dr_double_dot + a*dr + b = 0
        def diff_eq(t, P):  # P=[dx, dy, dz, du, dv, dw] See source[4]

            rx = P[0]
            ry = P[1]
            rz = P[2]

            ru = P[3]
            rv = P[4]
            rw = P[5]

            P_prime = [0, 0, 0, 0, 0, 0]

            P_prime[0] = ru
            P_prime[1] = rv
            P_prime[2] = rw

            P_prime[3] = -a * rx + b[0, 0]
            P_prime[4] = -a * ry + b[0, 1]
            P_prime[5] = -a * rz + b[0, 2]
            return P_prime

        initial_conditions = [0, 0, 0, 0, 0, 0]
        time_steps = np.linspace(time_passed, time_passed + 2 * dt, 2)

        result = solve_ivp(diff_eq, [time_passed, time_passed + 2 * dt], initial_conditions, 'DOP853', time_steps)
        # solve_ivp(ode,[0,t_int],initial_conditions,'LSODA',t_list)
        # print(result.y)
        dr = [result.y[0][1], result.y[1][1], result.y[2][1]]
        dv = [result.y[3][1], result.y[4][1], result.y[5][1]]

        # print('dr:')
        # print(dr)

        r_osc_i, v_osc_i=step_forward_orbit(ic_r0, ic_v0, dt, mu, error)
        #print('r_osc_i')
        #print(r_osc_i)

        r_i=r_osc_i+dr #new position vector to be graphed
        v_i=v_osc_i+dv

        #print('i:')
        #print(i)
        #print('r_i')
        #print(r_i)

        ic_r0 = r_i # becomes initial condition
        ic_v0 = v_i

        # Save data
        data[i, 0] = r_i[0]
        data[i, 1] = r_i[1]
        data[i, 2] = r_i[2]
        data[i, 3] = v_i[0]
        data[i, 4] = v_i[1]
        data[i, 5] = v_i[2]
        time[i, 0] = time_passed
        #print('time')
        #print(time)
        #print('pos')
        #print(data)


    print('############# DATA #################')
    print(data)
    print('#######################')

    return data, time

#Plotting function
def plot_position(data, time):
    testing_plot = plt.figure()
    l=iterations

    velocity_input=np.concatenate(([data[:,3]],[data[:,4]],[data[:,5]]),axis=0)
    position_input=np.concatenate(([data[:,0]],[data[:,1]],[data[:,2]]),axis=0)
    velocity_input=np.transpose(velocity_input)
    position_input=np.transpose(position_input)

    ########EARTH
    rho = radius  # radius of earth
    theta = np.linspace(0, 2 * np.pi, l)
    phi = np.linspace(0, np.pi, l)

    a = rho * np.outer(np.cos(theta), np.sin(phi))
    b = rho * np.outer(np.sin(theta), np.sin(phi))
    c = rho * np.outer(np.ones(l), np.cos(phi))

    earth = testing_plot.add_subplot(1,2,1, projection='3d')
    earth.plot_surface(a, b, c)

    #########PLOT DATA
    earth.scatter(position_input[:,0], position_input[:,1], position_input[:,2], s=1, c='red', linestyle='-')
    earth.set_title('Satellite Position')
    zoom=10000
    earth.set(xlim=(-zoom, zoom), ylim=(-zoom, zoom), zlim=(-zoom, zoom))

    velocity = testing_plot.add_subplot(1,2,2)
    speed=np.sqrt((np.abs(velocity_input[:,0]))**2+(np.abs(velocity_input[:,1]))**2+(np.abs(velocity_input[:,2]))**2)
    velocity.plot(time,speed)
    velocity.set_title('Satellite Speed')
    velocity.set_xlabel('Time (s)')
    velocity.set_ylabel('Speed (km/s)')

    plt.show()
    return True


def orbital_elements(r_vector, v_vector, time):

    r=magnitude(r_vector)
    v=magnitude(v_vector)

    v_r=(1/r)*np.dot(r_vector,v_vector)

    h_vector=np.cross(r_vector,v_vector)
    h=magnitude(h_vector)

    energy=(np.dot(v_vector,v_vector)/2)-mu/r
    semi_maj_axis=-mu/(2*energy)

    inclination=np.arccos(h_vector[2]/h)*(180/np.pi)

    K_hat=np.array([0, 0, 1])
    print(h_vector)
    node_line=np.cross(np.transpose(K_hat),np.transpose(h_vector))
    node=magnitude(node_line)

    if node_line[1]>=0:
        raan=np.arccos(node_line[0]/node)*(180/np.pi)
    else:
        raan=360-np.arccos(node_line[0]/node)*(180/np.pi)

    ecc_vector=(1/mu)*((v*v-mu/r)*r_vector-r*v_r*v_vector)
    eccentricity=magnitude(ecc_vector)

    if ecc_vector[2]>=0:
        periapsis=np.arccos(np.dot(node_line,ecc_vector)/(node*eccentricity))*(180/np.pi)
    else:
        periapsis = 360-np.arccos(np.dot(node_line, ecc_vector) / (node * eccentricity))*(180/np.pi)

    if v_r>=0:
        theta=np.arccos(np.dot(ecc_vector/eccentricity,r_vector/r))*(180/np.pi)
    else:
        theta = 360-np.arccos(np.dot(ecc_vector / eccentricity, r_vector / r))*(180/np.pi)

    theta=theta/(180/np.pi)
    E=2*np.arctan((np.sqrt((1-eccentricity)/(1+eccentricity)))*np.tan(theta/2))
    t0=time-(E-eccentricity*np.sin(E))*np.sqrt((semi_maj_axis**3)/mu)


    return [eccentricity, semi_maj_axis, inclination, raan, periapsis, t0, h]

def plot_orb(orb, time):
    elements, x=plt.subplots(nrows=2,ncols=3)
    print(orb[1, :])
    x[0,0].plot(time, orb[:,0])  #eccentricity
    x[0, 0].set_ylabel('eccentricity')
    x[0, 0].set_xlabel('time (s)')

    x[0, 1].plot(time, orb[:, 1]) #semi_maj_axis
    x[0, 1].set_ylabel('semi-major axis [km]')
    x[0, 1].set_xlabel('time (s)')
    print('sma:')
    print(orb[:, 1])


    x[0, 2].plot( time,orb[:, 2] ) #inclination
    x[0, 2].set_ylabel('inclination [deg]')
    x[0, 2].set_xlabel('time (s)')
    print('inclination:')
    print(orb[:, 2])

    x[1, 0].plot( time, orb[:, 3]) #raan
    x[1, 0].set_ylabel('raan [deg]')
    x[1, 0].set_xlabel('time (s)')
    print('raan')
    print(orb[:, 3])


    x[1, 1].plot( time, np.unwrap(orb[:, 4],discont=180)) #periapsis
    x[1, 1].set_ylabel('argument of periapsis [deg]')
    x[1, 1].set_xlabel('time (s)')
    print(orb[:, 4])

    x[1, 2].plot( time, np.unwrap(orb[:, 5], discont=180)) #t0
    x[1, 2].set_ylabel('t0 [deg]')
    x[1, 2].set_xlabel('time (s)')


    elements.tight_layout()
    plt.show()

    return True


data,time=encke(r0, v0, dt, mu, error, iterations)
plot_position(data, time)
#convert position and velocity to elements

velocity_input = np.concatenate(([data[:, 3]], [data[:, 4]], [data[:, 5]]), axis=0)
position_input = np.concatenate(([data[:, 0]], [data[:, 1]], [data[:, 2]]), axis=0)
velocity_input = np.transpose(velocity_input)
position_input = np.transpose(position_input)

orb=np.zeros((iterations,7))

for i in range(0,iterations):
    orb[i,:]=orbital_elements(position_input[i,:], velocity_input[i,:], time[i,0])

plot_orb(orb, time)

print('end')



"""
RESOURCES
[1] Curtis, H.D. Orbital Mechanics for Engineering Students, 2014.
    (Algorithm 3.4, 3.3, 12.1, 4.2)

[2] https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-346-astrodynamics-fall-2008/lecture-notes/lec_23.pdf

[3] Vallado, D. A. Fundamentals of Astrodynamics and Applications
    (pg. 473)
    
[4] https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html

[5] https://matplotlib.org/3.3.3/gallery/lines_bars_and_markers/spectrum_demo.html#sphx-glr-gallery-lines-bars-and-markers-spectrum-demo-py
"""

### ALGORITHM 4.2 pg. 197 to graph orbital elements


