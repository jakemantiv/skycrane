clearvars; clc

% where to save the data? 
dataFileName = 'mcTestData';

% How many Monte Carlo Runs?
numMC = 10;

% fix random seed
rng(0)

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

% preallocate memory 
X = zeros(6,numel(time), numMC);
Y = zeros(4,numel(time), numMC);
U = zeros(2,numel(time), numMC);
X0 = zeros(6,numMC);
for i = 1:numMC
    % TODO: Draw the initial condition from a random distribution
    delX0 = [-2; 1; 1; 2; -.1; -.02];
    X0(:,i) = delX0 + X_nom;
    [X(:,:,i),Y(:,:,i),U(:,:,i)] = truthModel(time,F,H,X0(:,i),X_nom,U_nom,control);
end

% collect outputs
data.X = X; 
data.Y = Y; 
data.U = U;
data.X0 = X0;
data.time = time;

filePath = ['..', filesep, 'lib', filesep];
save([filePath, dataFileName], 'data');