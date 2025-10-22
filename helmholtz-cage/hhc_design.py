import magpylib as magpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import scipy.optimize as opt

### Constants ###
################################################
MU_0 = 1.256637061E-6 # H/m
################################################

### Design Requirements ###
################################################
MAX_OUTPUT = 150E-6 # T
FIELD_WIDTH = 0.5 # m
FIELD_LIMIT = 0.02 # decimal %
################################################

def B_sq_coil(n, I, a, z, p=0):
    """
    Calculate the magnetic field output of a square electromagnetic coil.
    Args:
        n (int): Number of loops
        I (float): Supplied current [A]
        a (float): Side length [m]
        z (float): (nx1) Sampling positions along centerline [m]
        p (float): Position of coil along centerline [m]
    Returns:
        float: (nx1) Output field magnitude along coil centerline [T]
    """
    z = z - p
    B = (MU_0 * n * I * a**2) / (np.pi * (z**2 + ((a**2)/4)) * np.sqrt(z**2 + ((a**2)/2)))
    return B

def s_sq_pair(a):
    """
    Calculate the separation distance for a square Helmholtz pair.
    Args:
        a (float): Side length [m]
    Returns:
        float: Separation between coils along the centerlines [m]
    """
    s = 0.5445 * a
    return s

def B_sq_pair(n, I, a, z):
    """
    Calculate the magnetic field output of a pair of square electromagnetic coils.
    Args:
        n (int): Number of loops per coil
        I (float): Supplied current [A]
        a (float): Side length [m]
        z (float): (nx1) Sampling positions along centerline [m]
    Returns:
        float: (nx1) Output field magnitude along centerline [T]
    """
    s = s_sq_pair(a)
    Bl = B_sq_coil(n=n, I=I, a=a, z=z, p=-s/2)
    Br = B_sq_coil(n=n, I=I, a=a, z=z, p=s/2)
    B = Bl + Br
    return B
 
def sq_pair_geom(w, db):
    """
    Calculate the geometry of a square Helmholtz pair.
    Args:
        w (float): Width of the homogeneous field [m]
        db (float): Maximum change of magnetic field magnitude at edge of homogeneous field [dec %]
    Returns:
        float: Side length of coils [m]
        float: Separation between coils along the centerlines [m]
    """
    def f(a, z, db):
        B_c = B_sq_pair(n=10, I=1, a=a, z=0)
        B_e = B_sq_pair(n=10, I=1, a=a, z=z)
        delta = B_c * (1 - db)
        return B_e - delta

    z = w/2 # forward edge of the homogeneous limit [m]
    a = opt.brentq(f, a=0.5, b=3, xtol=1e-2, args=(z, db))
    s = s_sq_pair(a)
    return a, s
 
def sq_pair_scale(B, a):
    """
    Calculate the n*I scale factor required to meat the center field requirement
    Args:
        B (float): Central magnetic field magnitude target [T]
        a (float): Side length of the coils [m]
    Returns:
        float: Scale factor n*I required to achieve the desired output
    """
    B0 = B_sq_pair(n=1, I=1, a=a, z=0)
    k = B / B0
    return k

def sq_pair_polyline(n, I, a):
    """
    Construct the magpylib polyline collection corresponding to the square coil pair.
    Args:
        n (int): Number of loops.
        I (float): Supply current [A].
        a (float): Side length.
    Return:
        magpy.Collection: magpylib Collection object for calculating field output.
    """
    s = s_sq_pair(a)
    coil1 = magpy.Collection()
    for loop in range(n):
        vertices = np.array([[a/2, a/2, 0],
                             [a/2, -a/2, 0],
                             [-a/2, -a/2, 0],
                             [-a/2, a/2, 0],
                             [a/2, a/2, 0]])
        winding = magpy.current.Polyline(current=I, vertices=vertices)
        coil1.add(winding)
    coil2 = coil1.copy(position=(0, 0, s/2))
    coil1.position = (0, 0, -s/2)
    pair = magpy.Collection(coil1, coil2)
    return pair

### Demo: Field Taper Invariance to Num Loops ###
################################################
s = s_sq_pair(a=2)
fig, ax = plt.subplots()
for n in range(1,100,10):
    z = np.linspace(-1, 1, 100)
    B = B_sq_pair(n=n, I=1, a=2, z=z)
    homog_i = np.where(B>np.max(B) * (1-FIELD_LIMIT))
    ax.plot(z, B)
    ax.scatter(z[homog_i], B[homog_i], s=1.5)
ax.set(title='Field Width Invariance to Num Loops',
       xlabel='Centerline Position [m]',
       ylabel='|B| [T]')
#plt.show()
################################################

### Demo: Field Taper Invariance to Supply Current ###
################################################
s = s_sq_pair(a=2)
fig, ax = plt.subplots()
for I in range(1,20,1):
    z = np.linspace(-1,1,100)
    B = B_sq_pair(n=10, I=I, a=2, z=z)
    homog_i = np.where(B>np.max(B) * (1-FIELD_LIMIT))
    ax.plot(z, B)
    ax.scatter(z[homog_i], B[homog_i], s=1.5)
ax.set(title='Field Width Invariance to Supply Current',
       xlabel='Centerline Position [m]',
       ylabel='|B| [T]')
#plt.show()
################################################

### Find Required Side Length ###
################################################
a, s = sq_pair_geom(FIELD_WIDTH, FIELD_LIMIT)
print(f"Side length: {a} m.")
print(f"Separation: {s} m.")
z = np.linspace(-1,1,100)
B = B_sq_pair(n=10, I=I, a=a, z=z)
homog_i = np.where(B>np.max(B) * (1-FIELD_LIMIT))
fig, ax = plt.subplots()
ax.plot(z, B)
ax.hlines(np.max(B) * (1-FIELD_LIMIT), -1, 1, color='k', linestyle='--')
ax.vlines(z[homog_i[0][0]], 0, 0.001, color='k', linestyle='--')
ax.vlines(z[homog_i[0][-1]], 0, 0.001, color='k', linestyle='--')
ax.vlines(-s/2, 0, 0.001, color='k')
ax.vlines(s/2, 0, 0.001, color='k')
ax.set(title='Centerline Field Shape for Selected Geometry',
       xlabel='Centerline Position [m]',
       ylabel='|B| [T]')
#plt.show()
################################################

### Find Required Scale Factor ###
################################################
nI = sq_pair_scale(MAX_OUTPUT, a)
print(f"Scale Factor nI: {nI} A")
n = np.arange(1,100)
I = nI / n
l = 8 * n
# Create a cost function to optimize between supply current and wire length
W = (1, 1) # weights: can be tuned to reflect different considerations (finances, FS precision, etc)
# use normalized values for better comparability - assume 100 loop limit
J = lambda n, nI, W : W[0] * (nI/np.round(n)) / (nI/100) + W[1] * (8 * np.round(n)) / 100
n_select = int(opt.minimize(J, x0=20, args=(nI, W), method='Nelder-Mead', bounds=[(1, 100)]).x[0])
I_select = nI/n_select
print(f"Optimized num loops: {n_select}")
print(f"Optimized current supply: {I_select} A")
print(f"Optimized wire length: {8*n_select} m.")
fig, ax = plt.subplots()
twin = ax.twinx()
ax.scatter(n, I, s=2)
twin.scatter(n, l, s=2)
ax.vlines(n_select, 0, np.max(I), color='k')
ax.set(title='Num Loops vs Current & Winding Length',
       xlabel='Number of Windings',
       ylabel='Supply Current Required [A]')
twin.set_ylabel('Length of Wire Required [m]')
#plt.show()
################################################

### Extra Analysis & Validation with magpylib ###
################################################
pair = sq_pair_polyline(n_select, I_select, a)
grid = np.mgrid[0:0:1j, -1:1:100j, -1:1:100j].T[:,:,0]
_, Y, Z = np.moveaxis(grid, 2, 0)
B = pair.getB(grid)
_, By, Bz = np.moveaxis(B, 2, 0)
Bmag = np.linalg.norm(B, axis=2)
B0 = pair.getB((0,0,0))
B0amp = np.linalg.norm(B0)
err = np.linalg.norm((B-B0)/B0amp, axis=2)
err[np.where(err>1)] = 1 # place ceiling to improve data display 
# Show streamplot
fig, ax = plt.subplots(1, 1)
sp = ax.streamplot(Y, Z, By, Bz,
                   density=2,
                   color=Bmag,
                   norm=Normalize(vmin=0, vmax=MAX_OUTPUT*1.2, clip=True),
                   linewidth=np.sqrt(Bmag)*100,
                   cmap='inferno')
plt.colorbar(sp.lines, ax=ax, label='[T]')
plt.tight_layout()
ax.set(title='Magnetic Streamplot',
       xlabel='x [m]',
       ylabel='y [m]')
# Show magnitude homogeneity
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(Y, Z, err*100, levels=int(100/(FIELD_LIMIT*100)), cmap='inferno_r')
plt.colorbar(cp, ax=ax, label='% of B0')
plt.tight_layout()
ax.set(title='Magnetic Field Magnitude Error',
       xlabel='x [m]',
       ylabel='y [m]')
# Show 3D representation of coils
pair.show()
################################################

