%%% Clear workspace
clear
clc
close all

%%%Simulated noise parameters
% Process noise covariance
Q = 1e-7;
p_noise = [Q;Q;Q;Q;Q;Q;Q];
% Measurement noise covariance
Re = 5e-6;
m_noise = [Re;Re;Re;Re;0.0001*Re;0.0001*Re;0.0001*Re];
% Sampling time
Ts = 0.1; % [s]

%%%Get Earth Parameters for orbit
planet

%%%Cubesat parameters
global m invI I
m = 2.6; %%mass in kilograms
inertia

%%%Initial Conditions Position and Velocity
altitude = 408*1000; %%meters
x0 = R + altitude;
y0 = 0;
z0 = 0; 
xdot0 = 0;
inclination = deg2rad(56);
semi_major = norm([x0;y0;z0]); %%Orbit Radius for circular orbit
vcircular = sqrt(mu/semi_major); %%Orbital Speed from Vis-Viva Equation with a=r
ydot0 = vcircular*cos(inclination);
zdot0 = vcircular*sin(inclination);

%%%Initial Conditions Attitude and Angular velocity
%Initialize Euler Angles
phi0 = 0;
theta0 = 0;
psi0 = 0;
ptp0 = [phi0;theta0;psi0];
q0123_0 = Eu2Quat(ptp0); %Transform to quat initial conditions

%Initial Angular Velocity (Body Frame)
p0=0.01;
q0=0.01;
r0=0;

%Initial Conditions state vector
stateinitial = [q0123_0;p0;q0;r0];

%%%Orbit time parameters Circular Orbit
period = 2*pi/sqrt(mu)*semi_major^(3/2); %%Tcircular
number_of_orbits = 1;
tfinal = period*number_of_orbits;
tspan = [0,tfinal];

%%%This is where we integrate the equations of motion
[tout,stateout] = ode45(@Satellite,tspan,stateinitial);
%{
%%%Convert state to kilometers
stateout(:,1:6) = stateout(:,1:6)/1000;

%%%Extract the state vector
xout = stateout(:,1);
yout = stateout(:,2);
zout = stateout(:,3);
%}
q0123out = stateout(:,1:4);
ptpout = Quat2Eu(q0123out);
pqrout = stateout(:,5:7);
%{
%%%Make an Earth
[X,Y,Z] = sphere(100);
X = X*R/1000;
Y = Y*R/1000;
Z = Z*R/1000;

%%%Plot 3D orbit
fig = figure();
set(fig,'color','white')
plot3(xout,yout,zout,'b-','LineWidth',4)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
hold on
earth_sphere()
axis equal

%%%Plot Quaternions
fig2 = figure();
set(fig2,'color','white')
plot(tout,q0123out,'-','LineWidth',2);
grid on
xlabel('Time (sec)')
ylabel('Quaternions')
legend('q0','q1','q2','q3')

%%%Plot Euler Angles
fig3 = figure();
set(fig3,'color','white')
plot(tout,ptpout*180/pi,'-','LineWidth',2);
grid on
xlabel('Time (sec)')
ylabel('Euler Angles (deg)')
legend('Phi','Theta','Psi')

%%%Plot Angular Velocity
fig4 = figure();
set(fig4,'color','white')
plot(tout,pqrout,'-','LineWidth',2);
grid on
xlabel('Time (sec)')
ylabel('Angular Velocity (rad/s)')
legend('p','q','r')
%}

