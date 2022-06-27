%% Orbit Inputs
a0 = 6903000;
e0 = 0.003621614;
TA0 = deg2rad(0);
RAAN_0 = deg2rad(227.05566156);
Inc = deg2rad(97.49588373);
AoP = deg2rad(0.001);
t0 = 0;
coe = [a0;e0;Inc;RAAN_0;AoP;t0];
E_prev = 0;
J2 = 1.083e-03;
RE = 6371800; 
mu_earth = 3.986e+14;
mu_moon = 4.913e+12;
mu_sun = 1.32712e+20; 
mass = 4;
w_earth = [0;0;7.2e-05];

%[0.0378      0               0         ]
%[0               0.0375    -0.0001]
%[0             -0.0001    0.0162  ]

%% Disturbance Parameters
R_pm = [0.01 0.01 0.1]'; % in m
Cd = 2.1; %drag coeff
S = [0.0300  0.0300, 0.0100]'; % in m^2, cross sectional area
S_norm = norm(S);
c_rk = 1.5;
F_solar = 1366; % in W/m^2 %solar irradiance
altitude = a0 - RE; 
%[T, a, P, rho] = atmosisa(altitude);
year = 2024;
month = 06;
day = 01;
hours = 00;
minutes = 00;
seconds = 00;
jd0 = JDnumber(year,month,day,hours,minutes,seconds);
%% Dynamics Model Inputs
%I = diag([0.001731,0.001726,0.000264]); %Moment of inertia [kg-m2]
%I_inv = inv(I); % Inverse of Interia matrix
%W0 = [0 0 -1]'; % Initial angular velocity [rad/s]
%q0 = [1 0 0 0]';
%torques = [0 0 0]';
K_noise = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added for Aero/Solar Torque
%x_dim=0.1
%y_dim=0.1
%z_dim=0.3

wc_solar_area=0.12
wc_drag_area=0.042

sun_to_earth=[151*10^9, 151*10^9, 151*10^9]
speed_light_inverse=1/299792458
% Density of air at 500km
rho=10^(-12)

% DYNAMICS INPUTS
W0=[0 0 0]
I=[[37824169.23*10^-9 -25788.01*10^-9 32047.93*10^-9],
    [-25788.01*10^-9 37537377.89*10^-9 -51485.52*10^-9],
    [-32047.93*10^-9 -51485.52*10^-9, 16200348.41*10^-9]]
q0=[0 0 0 1]

function [qdot] = QDotSolver(w,q)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q_product = [
        q4, -q3, q2;
        q3, q4, -q1;
        -q2, q1, q4;
        -q1, -q2, -q3    
    ];
    qdot = 0.5*q_product*w;
end

function[jd] = JDnumber(year,month,day,hour,minute,second)
% Calculates the JD number for a given year,month and day

j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) + fix(275*month/9) + day + 1721013.5; 

ut = hour + minute/60 + second/3600; % UTC time in hour format

jd = j0 + ut/24; % JD number calculation
end
