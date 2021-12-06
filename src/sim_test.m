clearvars; clc; close all

% fix random seed
rng(0)

% where to save the data? 
dataFileName = 'mcTestData';

% Simulation Testing
control = true;
Tsim = 120;
dT = 0.1;
Nsim = 100;
qW = .011;

% Simulation time
time = 0:dT:Tsim;
N = numel(time);

% Get system parameters
skycrane_params

% Get truth parameters
data = load(['..',filesep,'lib', filesep, 'skycrane_finalproj_KFdata.mat']);
Qt = data.Qtrue;
Rt = data.Rtrue;

% Get non-linear functions
[Fnl,Hnl] = skycraneNL;

% Get linearized functions
[A,B,C,D,Gam,F,G,H,M,Om] = skycraneLin;

% Nominal conditions
x_nom  = 0;
z_nom  = 20;
a_nom  = 0;
dx_nom = 0;
dz_nom = 0;
da_nom = 0;
T1_nom = 0.5*g*(mb + mf)/cos(Be);
T2_nom = T1_nom;

% Nominal Points
n = 6; p = 4; m = 2;
Xnom = [x_nom, dx_nom, z_nom, dz_nom, a_nom, da_nom]';
Unom = [T1_nom, T2_nom]';
Ynom = Hnl(Xnom,Unom,zeros(p,1));

% Nominal Trajectory Functions
Xnom = @(t) Xnom;
Unom = @(t) Unom;
Ynom = @(t) Ynom;

% Initial Condition Variance
X_var = .05*[2,.5,2,.5,deg2rad(1),deg2rad(.2)]';

% Initialize Filter
P0 = diag([50,10,50,10,45,5]);
Q = @(dT,t) qW*diag([1,1,1]);
R = @(dt,t) Rt;

% MonteCarlo Simulation Runs
ex_lkf = zeros(Nsim,N-1);
ey_lkf = zeros(Nsim,N-1);
ex_ekf = zeros(Nsim,N-1);
ey_ekf = zeros(Nsim,N-1);
ex_ukf = zeros(Nsim,N-1);
ey_ukf = zeros(Nsim,N-1);
for i = 1:Nsim
    % Initial condition
    delX0 = randn(n,1).*X_var;
    X0 = delX0 + Xnom(0);
    
    % Truth model Simulation
    [X,Y,U] = truthModel(time,Fnl,Hnl,X0,Xnom,Unom,control);
    
    % Linearized Kalman Filter Estimate
    [Xh_lkf,Yh_lkf,P_lkf,S_lkf,Sx_lkf] = KF(time,Y,U,X0,P0,Xnom,Unom,Ynom,F,G,H,M,Om,Q,R);
    
    % Extended Kalman Filter
    [Xh_ekf,Yh_ekf,P_ekf,S_ekf,Sx_ekf] = EKF(time,Y,U,X0,P0,Xnom,Unom,Fnl,F,Hnl,H,Om,Q,R);
    
    % Unscented Kalman Filter
    [Xh_ukf,Yh_ukf,P_ukf,S_ukf,Sx_ukf] = UKF(time,Y,U,X0,0.1*P0,Fnl,Hnl,H,Om,Q,R);

    % Calculate NEES and NIS
    ex_lkf(i,:) = NEES(X,Xh_lkf,P_lkf);
    ey_lkf(i,:) = NIS(Y,Yh_lkf,S_lkf);
    ex_ekf(i,:) = NEES(X,Xh_ekf,P_ekf);
    ey_ekf(i,:) = NIS(Y,Yh_ekf,S_ekf);
    ex_ukf(i,:) = NEES(X,Xh_ukf,P_ukf);
    ey_ukf(i,:) = NIS(Y,Yh_ukf,S_ukf);
    
    % Progress Notification
    fprintf('Done with %d\n',i)
end

% NEES and NIS averaged over each time step
alpha = 0.05;
exb_lkf = mean(ex_lkf,1);
eyb_lkf = mean(ey_lkf,1);
exb_ekf = mean(ex_ekf,1);
eyb_ekf = mean(ey_ekf,1);
exb_ukf = mean(ex_ukf,1);
eyb_ukf = mean(ey_ukf,1);
rx = [chi2inv(alpha/2,Nsim*n), chi2inv(1 - alpha/2,Nsim*n)]'/Nsim;
ry = [chi2inv(alpha/2,Nsim*p), chi2inv(1 - alpha/2,Nsim*p)]'/Nsim;

% print some outputs of NEES/NIS tests to screen
% ratio of samples that were within ry bounds, should equal 1-alpha
nis_succss_ekf = sum(exb_ekf<rx(2) & exb_ekf>rx(1))./numel(exb_ekf);
nis_succss_lkf = sum(exb_lkf<rx(2) & exb_lkf>rx(1))./numel(exb_lkf);
nis_succss_ukf = sum(exb_ukf<rx(2) & exb_ukf>rx(1))./numel(exb_ukf);

% ratio of samples that were within ry bounds, should equal 1-alpha
nees_succss_ekf = sum(eyb_ekf<ry(2) & eyb_ekf>ry(1))./numel(eyb_ekf);
nees_succss_lkf = sum(eyb_lkf<ry(2) & eyb_lkf>ry(1))./numel(eyb_lkf);
nees_succss_ukf = sum(eyb_ukf<ry(2) & eyb_ukf>ry(1))./numel(eyb_ukf);

fprintf('LKF NIS Test:\n');
fprintf('Observed ratio between bounds: %f\n', nis_succss_lkf);
fprintf('Expected ratio between bounds: %f\n\n', 1-alpha);

fprintf('LKF NEES Test:\n');
fprintf('Observed ratio between bounds: %f\n', nees_succss_lkf);
fprintf('Expected ratio between bounds: %f\n\n', 1-alpha);

fprintf('EKF NEES Test:\n');
fprintf('Observed ratio between bounds: %f\n', nis_succss_ekf);
fprintf('Expected ratio between bounds: %f\n\n', 1-alpha);

fprintf('EKF NEES Test:\n');
fprintf('Observed ratio between bounds: %f\n', nees_succss_ekf);
fprintf('Expected ratio between bounds: %f\n\n', 1-alpha);

fprintf('UKF NEES Test:\n');
fprintf('Observed ratio between bounds: %f\n', nis_succss_ukf);
fprintf('Expected ratio between bounds: %f\n\n', 1-alpha);

fprintf('UKF NEES Test:\n');
fprintf('Observed ratio between bounds: %f\n', nees_succss_ukf);
fprintf('Expected ratio between bounds: %f\n\n', 1-alpha);

%% LKF plots
% Plot NIS and NEES Statistics
figure
subplot(2,1,1);
plot(time(2:end),exb_lkf,'or',time,ones(2,numel(time)).*rx,'--r')
ylabel('NEES')
subplot(2,1,2);
plot(time(2:end),eyb_lkf,'or',time,ones(2,numel(time)).*ry,'--r')
ylabel('NIS')

% Plot options for states
state_opts = struct;
state_opts.symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
state_opts.title = 'Simulated System States';
state_opts.saveFigs = false;
state_opts.filename = '';
state_opts.legends = {'Truth Sim','Linearized Kalman Filter'};

% Plot states
make_plots(state_opts,time,X,Xh_lkf)

% Plot states and error covariance
state_opts.title = 'Simulated System State Errors and 2-Sigma Bounds';
state_opts.legends = {'Truth Sim','+2 sigma','-2 sigma'};
make_plots(state_opts,time,X-Xh_lkf,2*Sx_lkf,-2*Sx_lkf)


% Plot options for measurements
meas_opts = struct;
meas_opts.symbols = {'$\xi$','$z$','$\dot{\theta}$','$\ddot{\xi}$'};
meas_opts.title = 'Simulated System Measurements';
meas_opts.saveFigs = false;
meas_opts.filename = '';
meas_opts.legends = {'Truth Sim','Linearized Kalman Filter'};

% Plot Measurements
make_plots(meas_opts,time,Y,Yh_lkf)
%% EKF plots
% Plot NIS and NEES Statistics
figure
subplot(2,1,1);
plot(time(2:end),exb_ekf,'or',time,ones(2,numel(time)).*rx,'--r')
ylabel('NEES')
subplot(2,1,2);
plot(time(2:end),eyb_ekf,'or',time,ones(2,numel(time)).*ry,'--r')
ylabel('NIS')

% Plot options for states
state_opts = struct;
state_opts.symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
state_opts.title = 'Simulated System States';
state_opts.saveFigs = false;
state_opts.filename = '';
state_opts.legends = {'Truth Sim','Extended Kalman Filter'};

% Plot states
make_plots(state_opts,time,X,Xh_ekf)

% Plot states and error covariance
state_opts.title = 'Simulated System State Errors and 2-Sigma Bounds';
state_opts.legends = {'Truth Sim','+2 sigma','-2 sigma'};
make_plots(state_opts,time,X-Xh_ekf,2*Sx_ekf,-2*Sx_ekf)


% Plot options for states
meas_opts = struct;
meas_opts.symbols = {'$\xi$','$z$','$\dot{\theta}$','$\ddot{\xi}$'};
meas_opts.title = 'Simulated System Measurements';
meas_opts.saveFigs = false;
meas_opts.filename = '';
meas_opts.legends = {'Truth Sim','Extended Kalman Filter'};

% Plot Measurements
make_plots(meas_opts,time,Y,Yh_ekf)

%% UKF plots
% Plot NIS and NEES Statistics
figure
subplot(2,1,1);
plot(time(2:end),exb_ukf,'or',time,ones(2,numel(time)).*rx,'--r')
ylabel('NEES')
subplot(2,1,2);
plot(time(2:end),eyb_ukf,'or',time,ones(2,numel(time)).*ry,'--r')
ylabel('NIS')

% Plot options for states
state_opts = struct;
state_opts.symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
state_opts.title = 'UKF Simulated System States';
state_opts.saveFigs = false;
state_opts.filename = '';
state_opts.legends = {'Truth Sim','Extended Kalman Filter'};

% Plot states
make_plots(state_opts,time,X,Xh_ukf)

% Plot states and error covariance
state_opts.title = 'UKF Simulated System State Errors and 2-Sigma Bounds';
state_opts.legends = {'Truth Sim','+2 sigma','-2 sigma'};
make_plots(state_opts,time,X-Xh_ukf,2*Sx_ukf,-2*Sx_ukf)


% Plot options for states
meas_opts = struct;
meas_opts.symbols = {'$\xi$','$z$','$\dot{\theta}$','$\ddot{\xi}$'};
meas_opts.title = 'UKF Simulated System Measurements';
meas_opts.saveFigs = false;
meas_opts.filename = '';
meas_opts.legends = {'Truth Sim','Extended Kalman Filter'};

% Plot Measurements
make_plots(meas_opts,time,Y,Yh_ukf)


%% Provided Data Analysis

% Use tuned filters to estimate states from provided observation and
% control data.

% Provided data
time = data.tvec;
Y = data.ydata;
U = data.uhist;

% Initial Conditions
X0 = Xnom(0);
P0 = diag([50,10,50,10,45,5]);

% Linearized Kalman Filter Estimate
[Xh_lkf,Yh_lkf,P_lkf,S_lkf,Sx_lkf] = KF(time,Y,U,X0,P0,Xnom,Unom,Ynom,F,G,H,M,Om,Q,R);

% Extended Kalman Filter
[Xh_ekf,Yh_ekf,P_ekf,S_ekf,Sx_ekf] = EKF(time,Y,U,X0,P0,Xnom,Unom,Fnl,F,Hnl,H,Om,Q,R);

% Unscented Kalman Filter
[Xh_ukf,Yh_ukf,P_ukf,S_ukf,Sx_ukf] = UKF(time,Y,U,X0,0.1*P0,Fnl,Hnl,H,Om,Q,R);

%% LKF Plots
% Plot options for states
state_opts = struct;
state_opts.symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
state_opts.title = 'Simulated System States';
state_opts.saveFigs = false;
state_opts.filename = '';
state_opts.legends = {'LKF','EKF','UKF'};

% Plot states
make_plots(state_opts,time,Xh_lkf,Xh_ekf,Xh_ukf)

% Plot options for 2-sigma
state_opts.title = '2-Sigma Bounds for Simulated System States';

% Plot 2-sigma
make_plots(state_opts,time,2*Sx_lkf,2*Sx_ekf,2*Sx_ukf)
