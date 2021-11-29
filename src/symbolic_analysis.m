clearvars
clc
% Symbolic analysis
syms x1 x2 z1 z2 a1 a2 u1 u2 B mb mf Cd rho As Ab In wc hc g

eqn1 = x2;
eqn2 = (u1*(cos(B)*sin(a1) + sin(B)*cos(a1)) + u2*(cos(B)*sin(a1) - sin(B)*cos(a1)) - (1/2)*Cd*rho*(As*cos(a1 - atan(z2/x2)) + Ab*sin(a1 - atan(z2/x2)))*x2*sqrt(x2^2 + z2^2))/(mb + mf);
eqn3 = z2;
eqn4 = (u1*(cos(B)*cos(a1) - sin(B)*sin(a1)) + u2*(cos(B)*cos(a1) + sin(B)*sin(a1)) - (1/2)*Cd*rho*(As*cos(a1 - atan(z2/x2)) + Ab*sin(a1 - atan(z2/x2)))*z2*sqrt(x2^2 + z2^2))/(mb + mf);
eqn5 = a2;
eqn6 = (1/In)*((u1 - u2)*cos(B)*wc + (u2 - u1)*sin(B)*hc);

J(1,1) = diff(eqn1,x1);
J(1,2) = diff(eqn1,x2);
J(1,3) = diff(eqn1,z1);
J(1,4) = diff(eqn1,z2);
J(1,5) = diff(eqn1,a1);
J(1,6) = diff(eqn1,a2);
R(1,1) = diff(eqn1,u1);
R(1,2) = diff(eqn1,u2);

J(2,1) = diff(eqn2,x1);
J(2,2) = diff(eqn2,x2);
J(2,3) = diff(eqn2,z1);
J(2,4) = diff(eqn2,z2);
J(2,5) = diff(eqn2,a1);
J(2,6) = diff(eqn2,a2);
R(2,1) = diff(eqn2,u1);
R(2,2) = diff(eqn2,u2);

J(3,1) = diff(eqn3,x1);
J(3,2) = diff(eqn3,x2);
J(3,3) = diff(eqn3,z1);
J(3,4) = diff(eqn3,z2);
J(3,5) = diff(eqn3,a1);
J(3,6) = diff(eqn3,a2);
R(3,1) = diff(eqn3,u1);
R(3,2) = diff(eqn3,u2);

J(4,1) = diff(eqn4,x1);
J(4,2) = diff(eqn4,x2);
J(4,3) = diff(eqn4,z1);
J(4,4) = diff(eqn4,z2);
J(4,5) = diff(eqn4,a1);
J(4,6) = diff(eqn4,a2);
R(4,1) = diff(eqn4,u1);
R(4,2) = diff(eqn4,u2);

J(5,1) = diff(eqn5,x1);
J(5,2) = diff(eqn5,x2);
J(5,3) = diff(eqn5,z1);
J(5,4) = diff(eqn5,z2);
J(5,5) = diff(eqn5,a1);
J(5,6) = diff(eqn5,a2);
R(5,1) = diff(eqn5,u1);
R(5,2) = diff(eqn5,u2);

J(6,1) = diff(eqn6,x1);
J(6,2) = diff(eqn6,x2);
J(6,3) = diff(eqn6,z1);
J(6,4) = diff(eqn6,z2);
J(6,5) = diff(eqn6,a1);
J(6,6) = diff(eqn6,a2);
R(6,1) = diff(eqn6,u1);
R(6,2) = diff(eqn6,u2);

J = subs(J,u1,0.5*g*(mb + mf)/cos(B));
J = subs(J,u2,0.5*g*(mb + mf)/cos(B));
J = subs(J,x1,0);
J = subs(J,z1,20);
J = subs(J,z2,0);
J = subs(J,a1,0);
J = subs(J,a2,0);
J = limit(J,x2,0);

R_sym = R 
R = subs(R,u1,0.5*g*(mb + mf)/cos(B));
R = subs(R,u2,0.5*g*(mb + mf)/cos(B));
R = subs(R,x1,0);
R = subs(R,z1,20);
R = subs(R,z2,0);
R = subs(R,a1,0);
R = subs(R,a2,0);
R = limit(R,x2,0);

J = subs(J,g,3.711);
R = subs(R,g,3.711);

R = subs(R,B,pi/4);
R = subs(R,wc,3.2/2);
R = subs(R,hc,0.9421);
R = subs(R,In,(1/12)*(1510*(3.2^2 + 2.5^2) + 390*(1^2 + 0.5^2)));
R = subs(R,mf,390);
R = subs(R,mb,1510);

double(J)
double(R)
% System Parameters
rho = 0.020;            % Mars surface atmosphere density           [kg/m3]
g = 3.711;              % Mars surface gravitational constant       [m/s2]
B = pi/4;               % Thruster mounting angle                   [rad]
Cd = 0.2;               % Drag Coefficient                          [-]

mf = 390;               % Propellant Mass                           [kg]
wf = 1;                 % Width of propellant housing               [m]
hf = 0.5;               % Height of propellant housing              [m]
df = 1;                 % Depth of propellant housing               [m]

mb = 1510;              % Mass of aircraft and payload              [kg]
wb = 3.2;               % Width of skycrane body                    [m]
hb = 2.5;               % Height of skycrane body                   [m]
db = 2.9;               % Depth of skycrane body                    [m]

wc = wb/2;              % Vehicle center of mass width dimension    [m]
hc = 0.9421;            % Vehicle center of mass height dimension   [m]

As = (hb*db) + (hf*df); % Side Area                                 [m2]
Ab = (wb*db) + (wf*df); % Bottom Area                               [m2]


% Moment of Inertia
In = (1/12)*(mb*(wb^2 + hb^2) + mf*(wf^2 + hf^2));


