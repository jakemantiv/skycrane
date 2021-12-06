function [A,B,C,D,L,F,G,H,M,O] = skycraneLin
% Get system parameters
skycrane_params

% Moment of Inertia
In = (1/12)*(mb*(wb^2 + hb^2) + mf*(wf^2 + hf^2));

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
dx2da  = @(X,U) (U(1)*(cos(Be)*cos(X(5)) - sin(Be)*sin(X(5))) + U(2)*(cos(Be)*cos(X(5)) + sin(Be)*sin(X(5))) - (1/2)*Cd*rho*((As*cos(X(5)) + Ab*sin(X(5)))*X(2)*X(4) + (-As*sin(X(5)) + Ab*cos(X(5)))*X(2)^2))/(mb + mf);
dx2da2 = @(X,U) 0;
dx2du1 = @(X,U) (cos(Be)*sin(X(5)) + sin(Be)*cos(X(5)))/(mb + mf);
dx2du2 = @(X,U) (cos(Be)*sin(X(5)) - sin(Be)*cos(X(5)))/(mb + mf);

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
dz2da  = @(X,U) (U(1)*(-cos(Be)*sin(X(5)) - sin(Be)*cos(X(5))) + U(2)*(-cos(Be)*sin(X(5)) + sin(Be)*cos(X(5))) - (1/2)*Cd*rho*((As*cos(X(5)) + Ab*sin(X(5)))*X(4)^2 + (-As*sin(X(5)) + Ab*cos(X(5)))*X(2)*X(4)))/(mb + mf);
dz2da2 = @(X,U) 0;
dz2du1 = @(X,U) (cos(Be)*cos(X(5)) - sin(Be)*sin(X(5)))/(mb + mf);
dz2du2 = @(X,U) (cos(Be)*cos(X(5)) + sin(Be)*sin(X(5)))/(mb + mf);

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
da2du1 = @(X,U) (1/In)*(cos(Be)*wc - sin(Be)*hc);
da2du2 = @(X,U) (1/In)*(sin(Be)*hc - cos(Be)*wc);

% Linearized System State Matrix
A = @(X,U) [dxdx(X,U)   dxdx2(X,U)  dxdz(X,U)   dxdz2(X,U)  dxda(X,U)   dxda2(X,U)
            dx2dx(X,U)  dx2dx2(X,U) dx2dz(X,U)  dx2dz2(X,U) dx2da(X,U)  dx2da2(X,U)
            dzdx(X,U)   dzdx2(X,U)  dzdz(X,U)   dzdz2(X,U)  dzda(X,U)   dzda2(X,U)
            dz2dx(X,U)  dz2dx2(X,U) dz2dz(X,U)  dz2dz2(X,U) dz2da(X,U)  dz2da2(X,U)
            dadx(X,U)   dadx2(X,U)  dadz(X,U)   dadz2(X,U)  dada(X,U)   dada2(X,U)
            da2dx(X,U)  da2dx2(X,U) da2dz(X,U)  da2dz2(X,U) da2da(X,U)  da2da2(X,U)];

% Linearized System Control Matrix
B = @(X,U) [dxdu1(X,U)  dxdu2(X,U)
            dx2du1(X,U) dx2du2(X,U)
            dzdu1(X,U)  dzdu2(X,U)
            dz2du1(X,U) dz2du2(X,U)
            dadu1(X,U)  dadu2(X,U)
            da2du1(X,U) da2du2(X,U)];

% Linearized Measurement Matrix
C = @(X,U) [1           0           0           0           0           0
            0           0           1           0           0           0
            0           0           0           0           0           1
            dx2dx(X,U)  dx2dx2(X,U) dx2dz(X,U)  dx2dz2(X,U) dx2da(X,U)  dx2da2(X,U)];

            
% Linearized Feedforward Matrix
D = @(X,U) [0           0
            0           0
            0           0
            dx2du1(X,U) dx2du2(X,U)];

% Noise Translation Matrix
L = [0  0   0
    1   0   0
    0   0   0
    0   1   0
    0   0   0
    0   0   1];

% Discrete Time Estimate Functions
n = 6; m = 2;
Ah = @(X,U) [A(X,U), B(X,U); zeros(m,n), zeros(m,m)];
eAh = @(dt,X,U) expm(dt*Ah(X,U));
Sf = substruct('()',{[1:n];[1:n]}); %#ok<*NBRAK>
Sg = substruct('()',{[1:n];[(n+1):(n+m)]});
Sq = substruct('()',{[1:n];[(n+1):(2*n)]});
F = @(dt,t,X,U) subsref(eAh(dt,X,U),Sf);
G = @(dt,t,X,U) subsref(eAh(dt,X,U),Sg);
H = @(dt,t,X,U) C(X,U);
M = @(dt,t,X,U) D(X,U);
O = @(dt,t,X,U) dt*L;
% Z = @(dt,X,U,W) dt*[-A(X,U), L*W*L'; zeros(n,n), A(X,U)'];
% eZ = @(dt,X,U,W) expm(Z(dt,X,U,W));
% O = @(dt,t,X,U,W) F(dt,t,X,U)*subsref(eZ(dt,X,U,W),Sq);
