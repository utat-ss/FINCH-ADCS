function dstatedt = Satellite(t,state)
%%genericState = [q0123;p;q;r];
global m invI I

%%select states
q0123 = state(1:4);
p = state(5);
q = state(6);
r = state(7);
pqr = state(5:7);

%%%Rotational Kinematics
%Derivative of Quaternions
PQRMAT = [0 -p -q -r;p 0 r -q;q -r 0 p;r q -p 0];
q0123dot = 0.5*PQRMAT*q0123;

%%%Rotational Dynamics
%%%Total External Disturbance Moments
LMN = [0;0;0];
H = I*pqr;
pqrdot = invI*(LMN - cross(pqr,H));

%%%Return Derivatives State
dstatedt = [q0123dot;pqrdot];

%%%Discretization (check this)
%dt = 1; 
%dstatedt_dis = state + dstatedt*dt;
end