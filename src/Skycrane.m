% Skycrane
clearvars
clc; close all

% save plots? 
saveFigs = false;

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

% Nominal conditions
x_nom  = 0;
z_nom  = 20;
a_nom  = 0;
dx_nom = 0;
dz_nom = 0;
da_nom = 0;
T1_nom = 0.5*g*(mb + mf)/cos(B);
T2_nom = T1_nom;

% Assume mf is constant

% Moment of Inertia
In = (1/12)*(mb*(wb^2 + hb^2) + mf*(wf^2 + hf^2));

% Alpha function
al = @(X) atan2(X(4),X(2));

% Drag functions
Fdx = @(X) (1/2)*Cd*rho*(As*cos(X(5) - al(X)) + Ab*sin(X(5) - al(X)))*X(2)*sqrt(X(2)^2 + X(4)^2);
Fdz = @(X) (1/2)*Cd*rho*(As*cos(X(5) - al(X)) + Ab*sin(X(5) - al(X)))*X(4)*sqrt(X(2)^2 + X(4)^2);

% Non-Linear functions 
dx =  @(X,U) X(2);
dx2 = @(X,U) (U(1)*(cos(B)*sin(X(5)) + sin(B)*cos(X(5))) + U(2)*(cos(B)*sin(X(5)) - sin(B)*cos(X(5))) - Fdx(X))/(mb + mf);
dz =  @(X,U) X(4);
dz2 = @(X,U) (U(1)*(cos(B)*cos(X(5)) - sin(B)*sin(X(5))) + U(2)*(cos(B)*cos(X(5)) + sin(B)*sin(X(5))) - Fdz(X))/(mb + mf) - g;
da =  @(X,U) X(6);
da2 = @(X,U) (1/In)*((U(1) - U(2))*cos(B)*wc + (U(2) - U(1))*sin(B)*hc);
dX =  @(X,U) [dx(X,U), dx2(X,U), dz(X,U), dz2(X,U), da(X,U), da2(X,U)]';
Yf =  @(X,U) [X(1); X(3); X(6); dx2(X,U)];


% Jacobians
dxdx   = @(X,U) 0;
dxdx2  = @(X,U) 1;
dxdz   = @(X,U) 0;
dxdz2  = @(X,U) 0;
dxda   = @(X,U) 0;
dxda2  = @(X,U) 0;
dxdu1  = @(X,U) 0;
dxdu2  = @(X,U) 0;

dx2dx  = @(X,U) 0;
dx2dx2 = @(X,U) -Cd*rho/(2*(mb + mf))*((As*sin(X(5)) - Ab*cos(X(5)))*X(4) + 2*(As*cos(X(5)) + Ab*sin(X(5)))*X(2));
dx2dz  = @(X,U) 0;
dx2dz2 = @(X,U) -Cd*rho/(2*(mb + mf))*((As*sin(X(5)) - Ab*cos(X(5)))*X(2));
dx2da  = @(X,U) (U(1)*(cos(B)*cos(X(5)) - sin(B)*sin(X(5))) + U(2)*(cos(B)*cos(X(5)) + sin(B)*sin(X(5))) - (1/2)*Cd*rho*((As*cos(X(5)) + Ab*sin(X(5)))*X(2)*X(4) + (-As*sin(X(5)) + Ab*cos(X(5)))*X(2)^2))/(mb + mf);
dx2da2 = @(X,U) 0;
dx2du1 = @(X,U) (cos(B)*sin(X(5)) + sin(B)*cos(X(5)))/(mb + mf);
dx2du2 = @(X,U) (cos(B)*sin(X(5)) - sin(B)*cos(X(5)))/(mb + mf);

dzdx   = @(X,U) 0;
dzdx2  = @(X,U) 0;
dzdz   = @(X,U) 0;
dzdz2  = @(X,U) 1;
dzda   = @(X,U) 0;
dzda2  = @(X,U) 0;
dzdu1  = @(X,U) 0;
dzdu2  = @(X,U) 0;

dz2dx  = @(X,U) 0;
dz2dx2 = @(X,U) -Cd*rho/(2*(mb + mf))*((As*sin(X(5)) - Ab*cos(X(5)))*2*X(4) + (As*cos(X(5)) + Ab*sin(X(5)))*X(2));
dz2dz  = @(X,U) 0;
dz2dz2 = @(X,U) -Cd*rho/(2*(mb + mf))*((As*cos(X(5)) + Ab*sin(X(5)))*X(4));
dz2da  = @(X,U) (U(1)*(-cos(B)*sin(X(5)) - sin(B)*cos(X(5))) + U(2)*(-cos(B)*sin(X(5)) + sin(B)*cos(X(5))) - (1/2)*Cd*rho*((As*cos(X(5)) + Ab*sin(X(5)))*X(4)^2 + (-As*sin(X(5)) + Ab*cos(X(5)))*X(2)*X(4)))/(mb + mf);
dz2da2 = @(X,U) 0;
dz2du1 = @(X,U) (cos(B)*cos(X(5)) - sin(B)*sin(X(5)))/(mb + mf);
dz2du2 = @(X,U) (cos(B)*cos(X(5)) + sin(B)*sin(X(5)))/(mb + mf);

dadx   = @(X,U) 0;
dadx2  = @(X,U) 0;
dadz   = @(X,U) 0;
dadz2  = @(X,U) 0;
dada   = @(X,U) 0;
dada2  = @(X,U) 1;
dadu1  = @(X,U) 0;
dadu2  = @(X,U) 0;

da2dx  = @(X,U) 0;
da2dx2 = @(X,U) 0;
da2dz  = @(X,U) 0;
da2dz2 = @(X,U) 0;
da2da  = @(X,U) 0;
da2da2 = @(X,U) 0;
da2du1 = @(X,U) (1/In)*(cos(B)*wc - sin(B)*hc);
da2du2 = @(X,U) (1/In)*(sin(B)*hc - cos(B)*wc);

% Linearized System State Matrix
A_fcn = @(X,U) [dxdx(X,U)   dxdx2(X,U)  dxdz(X,U)   dxdz2(X,U)  dxda(X,U)   dxda2(X,U)
                dx2dx(X,U)  dx2dx2(X,U) dx2dz(X,U)  dx2dz2(X,U) dx2da(X,U)  dx2da2(X,U)
                dzdx(X,U)   dzdx2(X,U)  dzdz(X,U)   dzdz2(X,U)  dzda(X,U)   dzda2(X,U)
                dz2dx(X,U)  dz2dx2(X,U) dz2dz(X,U)  dz2dz2(X,U) dz2da(X,U)  dz2da2(X,U)
                dadx(X,U)   dadx2(X,U)  dadz(X,U)   dadz2(X,U)  dada(X,U)   dada2(X,U)
                da2dx(X,U)  da2dx2(X,U) da2dz(X,U)  da2dz2(X,U) da2da(X,U)  da2da2(X,U)];

% Linearized System Control Matrix
B_fcn = @(X,U) [dxdu1(X,U)  dxdu2(X,U)
                dx2du1(X,U) dx2du2(X,U)
                dzdu1(X,U)  dzdu2(X,U)
                dz2du1(X,U) dz2du2(X,U)
                dadu1(X,U)  dadu2(X,U)
                da2du1(X,U) da2du2(X,U)];

% Linearized Measurement Matrix
C_fcn = @(X,U) [1           0           0           0           0           0
                0           0           1           0           0           0
                0           0           0           0           0           1
                dx2dx(X,U)  dx2dx2(X,U) dx2dz(X,U)  dx2dz2(X,U) dx2da(X,U)  dx2da2(X,U)];

            
% Linearized Feedforward Matrix
D_fcn = @(X,U) [0           0
                0           0
                0           0
                dx2du1(X,U) dx2du2(X,U)];


% Nominal Points
X_nom = [x_nom, dx_nom, z_nom, dz_nom, a_nom, da_nom]';
U_nom = [T1_nom, T2_nom]';


% Linearize about nominal points
A = A_fcn(X_nom,U_nom);
B = B_fcn(X_nom,U_nom);
C = C_fcn(X_nom,U_nom);
D = D_fcn(X_nom,U_nom);


% Find discrete time state matrices
dt = 0.1;
[n,m] = size(B);
A_h = [A            B
       zeros(m,n)   zeros(m,m)];
eAh = expm(A_h*dt);
F = eAh(1:n,1:n);
G = eAh(1:n,n+1:end);
H = C;
M = D;


% Determine DT system observability
Obsv = [H; H*F; H*F^2; H*F^3; H*F^4; H*F^5];
is_obsv = rank(Obsv) == 6 

% Determine DT system stability
evals = eig(F)

% Determine DT system controllability
Co = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B];
is_ctrb = rank(Co) == 6

% Simulate linearized system about nominal point, with small perturbation
N = 500;
time = (0:N)*dt;
delX0 = [5, -.01, 10, -0.15, deg2rad(.02), -deg2rad(.005)]';
delX0 = [0; 0.2; 0; 0; 0; 0.001]
delX = delX0.*ones(n,N);
U0 = [0,0]';
for k = 1:N
    delX(:,k+1) = F*delX(:,k) + G*U0;
end
X = delX + X_nom;
Y = H*X + M*U0;


% Simulate non-linear system about same nominal point, with same small
% perturbation
sol = ode45(@(t,X) dX(X,U_nom),[0,N*dt],X_nom + delX0);
X_nl = deval(sol,time);
Y_nl = zeros(size(Y));
for k = 1:N+1
    Y_nl(:,k) = Yf(X_nl(:,k),U_nom);
end

% Plot states
symbols = {'$\xi$','$\dot{\xi}$','$z$','$\dot{z}$','$\theta$','$\dot{\theta}$'};
figure
for i = 1:n
    ax(i) = subplot(n,1,i);
    plot(ax(i),time,X(i,:),'-',time,X_nl(i,:),'--', 'LineWidth', 1.5);
    grid(ax(i),'on');grid(ax(i),'minor')
    ylabel(ax(i),symbols{i},'Interpreter','latex', 'FontSize', 20)
    legend('Nonlinear', 'Linearized', 'location', 'best');
end
xlabel('Time (s)');
linkaxes(ax,'x')
sgtitle('Simulated System States');
if saveFigs
    saveas(gcf,'../figs/simulated_states.png');
end

% Plot measurements
p = 4;
symbols = {'$\xi$','$z$','$\dot{\theta}$','$\ddot{\xi}$'};
figure
for i = 1:p
    ax(i) = subplot(p,1,i);
    plot(ax(i),time,Y(i,:),'-',time,Y_nl(i,:),'--', 'LineWidth', 1.5);
    grid(ax(i),'on');grid(ax(i),'minor')
    ylabel(ax(i),symbols{i},'Interpreter','latex', 'FontSize', 20)

end

xlabel('Time (s)');
linkaxes(ax,'x')
sgtitle('Simulated System Measurements');
if saveFigs
    saveas(gcf,'../figs/simulated_measurements.png');
end
% Plot x/z trajectories
figure
plot(X(1,:),X(3,:),'-',X_nl(1,:),X_nl(3,:),'--', 'LineWidth', 1.5)
grid on; grid minor
if saveFigs
    saveas(gcf,'../figs/xz_traj.png');
end

% 