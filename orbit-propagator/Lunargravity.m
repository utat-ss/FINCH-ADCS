function[r_moon] = Lunargravity(jd)
% Calculates the geocentric position of the moon 
% Reference : Howard Curtis, Chp 12. pages 707 - 708
% Moon Vector
RE = 6378100; % in m

T = (jd - 2451545.0)/36525; % Number of julian centuries since J2000

e = 23.439 - 0.0130042*T; % Obliquity of the ecliptic plane

% Ecliptic longitude of the moon
elong = 218.32 + 481267.881*T ...
+ 6.29*sin(135.0 + 477198.87*T) - 1.27*sin(259.3 - 413335.36*T)...
+ 0.66*sin(235.7 + 890534.22*T) + 0.21*sin(269.9 + 954397.74*T)...
- 0.19*sin(357.5 + 35999.05*T) - 0.11*sin(186.5 + 966404.03*T);
e_long = mod(elong,360);

% Ecliptic latitude of the moon
elat = 5.13*sin( 93.3 + 483202.02*T) + 0.28*sin(228.2 + 960400.89*T)...
- 0.28*sin(318.3 + 6003.15*T) - 0.17*sin(217.6 - 407332.21*T);
e_lat = mod(elat,360);

% Horizontal parallax of the moon (deg):
hpar =0.9508 ...
+ 0.0518*cos(135.0 + 477198.87*T) + 0.0095*cos(259.3 - 413335.36*T)...
+ 0.0078*cos(235.7 + 890534.22*T) + 0.0028*cos(269.9 + 954397.74*T);
h_par = mod(hpar,360);

%...Direction cosines of the moon’s geocentric equatorial position vector:
l = cos(e_lat) * cos(e_long);
m = cos(e)*cos(e_lat)*sin(e_long) - sin(e)*sin(e_lat);
n = sin(e)*cos(e_lat)*sin(e_long) + cos(e)*sin(e_lat);

%...Earth-moon distance (km):
dist = RE/sin(h_par);

%...Moon’s geocentric equatorial position vector (km):
r_moon = dist*[l;m;n];

end
