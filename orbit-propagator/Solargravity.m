function[r_sun,u] = Solargravity(jd)
%...Astronomical unit (m):
AU = 149597870691;
%...Julian days since J2000:
n = jd - 2451545;
%...Julian centuries since J2000:
cy = n/36525;
%...Mean anomaly (deg{:
M = 357.528 + 0.9856003*n;
M = mod(M,360);
M = deg2rad(M);
%...Mean longitude (deg):
L = 280.460 + 0.98564736*n;
L = mod(L,360);
L = deg2rad(L);
%...Apparent ecliptic longitude (deg):
lamda = L + 1.915*sin(M) + 0.020*sin(2*M);
lamda = mod(lamda,360);
lamda = deg2rad(lamda);
%...Obliquity of the ecliptic (deg):
eps = deg2rad(23.439 - 0.0000004*n);
%...Unit vector from earth to sun:
u = [cos(lamda); sin(lamda)*cos(eps); sin(lamda)*sin(eps)];
%...Distance from earth to sun (m):
rS = (1.00014 - 0.01671*cos(M) - 0.000140*cos(2*M))*AU;
%...Geocentric position vector (m):
r_sun = rS*u;
end