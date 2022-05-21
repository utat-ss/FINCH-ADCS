function [sys,x0,str,ts] = sf_dyn(t,x,u,flag)
% S-function sf_aerodyn.M
% This S-function represents the nonlinear aircraft dynamics

% Copyright 1986-2007 The MathWorks, Inc. 


switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
    otherwise
        ctrlMsgUtils.error('Controllib:general:UnexpectedError',['Unhandled flag = ',num2str(flag)]);

end

% end sfuntmpl


%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 7;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 7;
sizes.NumInputs      = 7;
%sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% Initial conditions
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
x0 = [q0123_0;p0;q0;r0];
%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
%%%Moments of Inertia Cubesat
Ix = 0.9;
Iy = 0.9;
Iz = 0.3;
I = [Ix,0,0;0,Iy,0;0,0,Iz]; %%kg-m^2
invI = inv(I);

%%select states
q0123 = x(1:4);
p = x(5);
q = x(6);
r = x(7);
pqr = x(5:7);
%%Add process noise
pqr = pqr + [u(1);u(2);u(3)];

%%%Rotational Kinematics
%Derivative of Quaternions
PQRMAT = [0 -p -q -r;p 0 r -q;q -r 0 p;r q -p 0];
q0123dot = 0.5*PQRMAT*q0123 + [u(4);u(5);u(6);u(7)];

%%%Rotational Dynamics
%%%Total External Disturbance Moments
LMN = [0;0;0];
H = I*pqr;
pqrdot = invI*(LMN - cross(pqr,H));

%%%Return Derivatives State
sys = [q0123dot;pqrdot];

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

sys = x;

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 0.1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate