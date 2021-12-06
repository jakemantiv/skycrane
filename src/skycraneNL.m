function [F,H] = skycraneNL

% Get System Parameters
skycrane_params

% TODO These anonymous functions are making MC runs very slow
% may speed things up if we un-anonymous these guys 

% Moment of Inertia
In = (1/12)*(mb*(wb^2 + hb^2) + mf*(wf^2 + hf^2));

% Alpha function
al = @(X) atan2(X(4),X(2));

% Drag functions
Fdx = @(X) (1/2)*Cd*rho*(As*cos(X(5) - al(X)) + Ab*sin(X(5) - al(X)))*X(2)*sqrt(X(2)^2 + X(4)^2);
Fdz = @(X) (1/2)*Cd*rho*(As*cos(X(5) - al(X)) + Ab*sin(X(5) - al(X)))*X(4)*sqrt(X(2)^2 + X(4)^2);

% Non-Linear functions 
dx =  @(X,U,W) X(2);
dx2 = @(X,U,W) (U(1)*(cos(Be)*sin(X(5)) + sin(Be)*cos(X(5))) + U(2)*(cos(Be)*sin(X(5)) - sin(Be)*cos(X(5))) - Fdx(X))/(mb + mf) + W(1);
dz =  @(X,U,W) X(4);
dz2 = @(X,U,W) (U(1)*(cos(Be)*cos(X(5)) - sin(Be)*sin(X(5))) + U(2)*(cos(Be)*cos(X(5)) + sin(Be)*sin(X(5))) - Fdz(X))/(mb + mf) - g + W(2);
da =  @(X,U,W) X(6);
da2 = @(X,U,W) (1/In)*((U(1) - U(2))*cos(Be)*wc + (U(2) - U(1))*sin(Be)*hc) + W(3);

F =  @(X,U,W) [ dx(X,U,W)
                dx2(X,U,W)
                dz(X,U,W)
                dz2(X,U,W)
                da(X,U,W)
                da2(X,U,W)];
            
H =  @(X,U,V) [ X(1) + V(1)
                X(3) + V(2)
                X(6) + V(3)
                dx2(X,U,zeros(3,1)) + V(4)];