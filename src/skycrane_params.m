function skycrane_params
%#ok<*NASGU> 

% System Parameters
rho = 0.020;            % Mars surface atmosphere density           [kg/m3]
g = 3.711;              % Mars surface gravitational constant       [m/s2]
Be = pi/4;              % Thruster mounting angle                   [rad]
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

vars = who;

for i = 1:numel(vars)
    evalin('caller',sprintf('%s = %f;',vars{i},eval(vars{i})));
end
