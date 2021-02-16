%% Orbit Inputs
a0 = 6721800;
e0 = 0.007;
TA0 = deg2rad(0);
RAAN_0 = deg2rad(45);
Inc = deg2rad(96.85);
AoP = deg2rad(90);
t0 = 0;
coe = [a0;e0;Inc;RAAN_0;AoP;t0];
E_prev = 0;
J2 = 1.083e-03;
RE = 6371800; 
mu_earth = 3.986e+14;
mu_moon = 4.913e+12;
mu_sun = 1.32712e+20; 
mass = 0.6;
w_earth = [0;0;7.2e-05];

%% Disturbance Parameters
R_pm = [0.0045 0.002 -0.0082]'; % in m
Cd = 2.2;
S = [0.00921 0.01229 0.00252]'; % in m^2
S_norm = norm(S);
c_rk = 1.5;
F_solar = 1366; % in W/m^2
altitude = a0 - RE; 
%[T, a, P, rho] = atmosisa(altitude);
year = 2020;
month = 06;
day = 01;
hours = 00;
minutes = 00;
seconds = 00;
jd0 = JDnumber(year,month,day,hours,minutes,seconds);
%% Dynamics Model Inputs
I = diag([0.001731,0.001726,0.000264]); %Moment of inertia [kg-m2]
I_inv = inv(I); % Inverse of Interia matrix
W0 = [0 0 -1]'; % Initial angular velocity [rad/s]
q0 = [1 0 0 0]';
%torques = [0 0 0]';
K_noise = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added for Aero/Solar Torque
x_dim=0.1
y_dim=0.1
z_dim=0.3

