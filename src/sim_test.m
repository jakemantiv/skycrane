clearvars; clc

% Simulation Testing
control = true;

% Get system parameters
skycrane_params

% Get non-linear functions
[F,H] = skycraneNL;

% Nominal conditions
x_nom  = 0;
z_nom  = 20;
a_nom  = 0;
dx_nom = 0;
dz_nom = 0;
da_nom = 0;
T1_nom = 0.5*g*(mb + mf)/cos(B);
T2_nom = T1_nom;

% Nominal Points
X_nom = [x_nom, dx_nom, z_nom, dz_nom, a_nom, da_nom]';
U_nom = [T1_nom, T2_nom]';

% Initial condition
delX0 = [-2; 1; 1; 2; -.1; -.02];
X0 = delX0 + X_nom;

% Simulation time
dT = 0.1;
time = 0:dT:50;
[X,Y,U] = truthModel(time,F,H,X0,X_nom,U_nom,control);

% Plot states
symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
title = 'Simulated System States';
make_plots(time,X,symbols,title,false);

% Plot Measurements
symbols = {'$\xi$','$z$','$\dot{\theta}$','$\ddot{\xi}$'};
title = 'Simulated System Measurements';
make_plots(time,Y,symbols,title,false);