clearvars; clc; close all

% Simulation Testing
control = true;
Tsim = 120;
dT = 0.1;
Nsim = 100;
qW = .0009;

% Simulation time
time = 0:dT:Tsim;
N = numel(time);

% Get system parameters
skycrane_params

% Get truth parameters
data = load("..\lib\skycrane_finalproj_KFdata.mat");
Qt = data.Qtrue;
Rt = data.Rtrue;

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

% Get non-linear functions
[Fnl,Hnl] = skycraneNL;
Y_nom = Hnl(X_nom,U_nom,zeros(4,1));

% Get linearized functions
[A,~,~,~,F,G,H,M,L] = skycraneLin(dT,X_nom,U_nom);
[p,n] = size(H);

% Initialize Filter
P0 = 1E4*eye(6);
W = qW*diag([1,1,1]);
Z = dT*[-A, L*W*L'; zeros(n), A'];
eZ = expm(Z);
Q = F*eZ(1:n,n+1:end);
R = Rt;
X_var = [2,.5,2,.5,deg2rad(1),deg2rad(.05)]';

% MonteCarlo Simulation Runs
ex = zeros(Nsim,N-1);
ey = zeros(Nsim,N-1);
for i = 1:Nsim
    % Initial condition
    delX0 = randn(n,1).*X_var;
    X0 = delX0 + X_nom;
    
    % Truth model Simulation
    [X,Y,U] = truthModel(time,U_nom,Fnl,Hnl,Qt,Rt,X0,X_nom,control);
    
    % Linearized Kalman Filter Estimate
    dY = Y - Y_nom;
    dU = U - U_nom;
    dX0 = zeros(n,1);
    [dXh,dYh,P,S,Sx] = KF(time,dY,dU,dX0,P0,F,G,H,M,Q,R);
    Xh = dXh + X_nom;
    Yh = dYh + Y_nom;
    
    % Calculate NEES and NIS at each time step
    ex(i,:) = NEES(X,Xh,P);
    ey(i,:) = NIS(Y,Yh,S);
    
    fprintf('Done with %d\n',i)
end

% NEES and NIS averaged over each time step
alpha = 0.05;
exb = mean(ex,1);
eyb = mean(ey,1);
rx = [chi2inv(alpha/2,Nsim*n), chi2inv(1 - alpha/2,Nsim*n)]'/Nsim;
ry = [chi2inv(alpha/2,Nsim*p), chi2inv(1 - alpha/2,Nsim*p)]'/Nsim;

% Plot NIS and NEES Statistics
figure
subplot(2,1,1);
plot(time(2:end),exb,'or',time,ones(2,numel(time)).*rx,'--r')
ylabel('NEES')
subplot(2,1,2);
plot(time(2:end),eyb,'or',time,ones(2,numel(time)).*ry,'--r')
ylabel('NIS')

% Plot options for states
state_opts = struct;
state_opts.symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
state_opts.title = 'Simulated System States';
state_opts.saveFigs = false;
state_opts.filename = '';
state_opts.legends = {'Truth Sim','Kalman Filter'};

% Plot states
make_plots(state_opts,time,X,Xh)

% Plot states and error covariance
state_opts.title = 'Simulated System State Errors and 2-Sigma Bounds';
state_opts.legends = {'Truth Sim','+2 sigma','-2 sigma'};
make_plots(state_opts,time,X-Xh,2*Sx,-2*Sx)


% Plot options for states
meas_opts = struct;
meas_opts.symbols = {'$\xi$','$z$','$\dot{\theta}$','$\ddot{\xi}$'};
meas_opts.title = 'Simulated System Measurements';
meas_opts.saveFigs = false;
meas_opts.filename = '';
meas_opts.legends = {'Truth Sim','Kalman Filter'};

% Plot Measurements
make_plots(meas_opts,time(2:end),Y,Yh)