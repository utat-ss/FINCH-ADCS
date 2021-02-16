import matplotlib.pyplot as plt
import numpy as np


class Orbit():
    def __init__(self,eccentricity, semi_maj_axis, inclination, raan, periapsis, t0, time_step, duration):
        self.velocity=np.zeros(3)
        self.position=np.zeros(3)

        self.time_step=time_step
        self.duration=duration
        self.current_t=t0

        self.e=eccentricity
        self.a=semi_maj_axis
        self.i=inclination
        self.raan=raan
        self.w=periapsis
        self.t0=t0

        self.G=6.67408*(10**-20)
        self.mass_earth=5.972*(10**24)
        self.mu=self.G*self.mass_earth


        ###Euler Angles
        self.C3_w=np.array([[np.cos(self.w),  np.sin(self.w), 0],
                          [-np.sin(self.w),  np.cos(self.w), 0],
                          [0,            0,          1]])

        self.C1_i=np.array([[1,     0,           0],
                          [0,     np.cos(self.i),  np.sin(self.i)],
                          [0,     -np.sin(self.i), np.cos(self.i)]])

        self.C3_omega=np.array([[np.cos(self.raan),  np.sin(self.raan), 0],
                          [-np.sin(self.raan),  np.cos(self.raan), 0],
                          [0,            0,          1]])
        t=np.matmul(self.C3_w, self.C1_i)
        self.geo_to_perifocal = np.matmul(t,self.C3_omega)
        self.perifocal_to_geo = np.linalg.inv(self.geo_to_perifocal)

    def euler_method(self, est, e, M, acc):
        while True:
            g_est= est-e*np.sin(est)-M
            g_prime_est=1-e*np.cos(est)
            new_est=est-(g_est)/(g_prime_est)
            if (np.abs(new_est-est)<acc):
                return new_est
            else:
                est=new_est

    def get_eccentric_anomaly(self):
        # M = E - e*sin(E)
        M = np.sqrt(self.mu/self.a**3)*(self.current_t-self.t0)
        E= self.euler_method(M, self.e, M, 0.0001)
        return E


    def get_true_anomaly(self, E):
        theta = 2*np.arctan(np.sqrt((1+self.e)/(1-self.e))*np.tan(E/2))
        return theta

    def get_distance(self, theta):
        l = self.a*(1-self.e**2)
        r = l/(1+self.e*np.cos(theta))
        return r

    def get_next_pos_vector(self, r, theta):
        r_vector = np.matmul(self.perifocal_to_geo,([r*np.cos(theta), r*np.sin(theta), 0]))
        return r_vector

    def get_next_velocity_vector(self, theta):

        b=np.sqrt(self.mu/(self.a*(1-self.e**2)))
        a=[-np.sin(theta)*b, (np.cos(theta)+self.e)*b, 0]
        v_vector = np.matmul(self.perifocal_to_geo,a)
        return v_vector

    def concat(self, pos_vector):
        self.position = np.append(self.position, pos_vector)
        return True

    def plot_position(self, x, y, z):
        testing_plot = plt.figure()
        l=len(x)

        ########EARTH
        rho = 6378  # radius of earth
        theta = np.linspace(0, 2 * np.pi, l)
        phi = np.linspace(0, np.pi, l)

        a = rho * np.outer(np.cos(theta), np.sin(phi))
        b = rho * np.outer(np.sin(theta), np.sin(phi))
        c = rho * np.outer(np.ones(l), np.cos(phi))

        earth = testing_plot.add_subplot(1,2,1, projection='3d')
        earth.plot_surface(a, b, c)

        #########PLOT DATA
        earth.scatter(x, y, z, c='red')
        earth.set_title('Satellite Position')
        zoom=9000
        earth.set(xlim=(-zoom, zoom), ylim=(-zoom, zoom), zlim=(-zoom, zoom))

        velocity = testing_plot.add_subplot(1,2,2)
        speed=np.sqrt((np.abs(self.velocity[:,0]))**2+(np.abs(self.velocity[:,1]))**2+(np.abs(self.velocity[:,2]))**2)
        velocity.plot(speed)
        velocity.set_title('Satellite Speed')
        velocity.set_xlabel('Time (s)')
        velocity.set_ylabel('Speed (km/s)')

        plt.show()
        return True

def vectors(eccentricity, semi_maj_axis, inclination, raan, periapsis, t0, time_step):
    duration=1
    orbit=Orbit(eccentricity, semi_maj_axis, inclination, raan, periapsis, t0, time_step, duration)


    orbit.current_t=orbit.current_t+orbit.time_step
    E=orbit.get_eccentric_anomaly()
    theta=orbit.get_true_anomaly(E)
    r=orbit.get_distance(theta)
    r_vector=orbit.get_next_pos_vector(r, theta)

    v_vector=orbit.get_next_velocity_vector(theta)

    return r_vector, v_vector



"""
Resources:
[1] https://matplotlib.org/gallery/mplot3d/scatter3d.html#sphx-glr-gallery-mplot3d-scatter3d-py
[2] https://matplotlib.org/gallery/mplot3d/surface3d_2.html#sphx-glr-gallery-mplot3d-surface3d-2-py
[3] Spacecraft Dynamics and Control
[4] Course Notes
[5] http://www.bogan.ca/orbits/kepler/keplerex.html
[6] https://www.school-for-champions.com/science/gravitation_universal_equation.htm#.XzyR8S0ZPUp

"""