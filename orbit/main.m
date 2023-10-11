% this is just a simple calculation for HERON Mk II

Apogee = 527.881;
Perigee = 514.644;
M = -84.795;
M = deg2rad(M);
eps = 0.000000008;

e = eccentricity(Apogee,Perigee); %
fprintf("eccentricity: %d :  \n",e);

p = semimajor_axis(Apogee,Perigee);
fprintf("semimajor axis: %d :  \n",p);

E = eccentric_anomaly(M,e,eps);
T = eccentric_to_true(E,e);
T = rad2deg(T);
fprintf("Mean Anomaly: %d :  \n", T);
