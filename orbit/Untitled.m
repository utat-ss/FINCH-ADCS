% this is just a simple calculation for HERON Mk II

Apogee = 508.987;
Perigee = 504.112;

e = eccentricity(Apogee,Perigee); %
fprintf("eccentricity: %d :  \n",e);

p = semimajor_axis(Apogee,Perigee);
fprintf("semimajor axis: %d :  \n",p);