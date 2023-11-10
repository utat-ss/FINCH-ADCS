% this is just a simple calculation for HERON Mk II

Apogee = 522.784; % mean_apogee_altitude_km:
Perigee = 518.649; % mean_perigee_altitude_km
M = -53.039;
M = deg2rad(M);
eps = 0.000000001;

e = eccentricity(Apogee,Perigee); %
fprintf("eccentricity: %d :  \n",e);

p = semimajor_axis(Apogee,Perigee);
fprintf("semimajor axis: %d :  \n",p);

E = eccentric_anomaly(M,e,eps);
T = eccentric_to_true(E,e);
T = rad2deg(T);
fprintf("True Anomaly: %d :  \n", T);
