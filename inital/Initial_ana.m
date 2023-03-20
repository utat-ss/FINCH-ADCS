% inital ADCS analysis.
GSD = 82.5; % Ground sampling distance
h = 500000; % height above ground m
mu_E = 3.986e14; % Earth’s gravitational constant unit m^3/s^2
arcsec = 3600;
R_E = 6378.1363*1000; % in m
nf = 4; %assuming FMC factor is 4

 % FINCH satellite velocity
disp("angular distance (arcsec)");
theta_pixel = rad2deg(atan(GSD/h));
pixel_arcsec = theta_pixel*arcsec;
disp( pixel_arcsec);

v_sat = sqrt(mu_E/(R_E+h)); % Satellite speed;
disp("FINCH speed (m/s)");
disp(v_sat)

vg = v_sat*(R_E/(R_E+h));
disp("ground speed (m/s)");
disp(vg)

vg_eff = vg/nf;

% exposure time
t_exp = GSD/vg_eff;

fprintf("exposure time for %4.2f meters in seconds \n", GSD);
disp(t_exp);
 
