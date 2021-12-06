% Kalman Filter Code

function [Xh,Yh,P,S,Sx] = EKF(t,Y,U,X0,P0,Fnl,Flin,G,Hnl,Hlin,M,Q,R)
[p,n] = size(Hlin);

% Initialize Filter
Xp = X0;
Pp = P0;

% Initialize Output Vectors
Xh = X0.*ones(size(X0,1),numel(t));
Yh = NaN(size(Y));
P = P0.*ones(size(X0,1),size(X0,1),numel(t));
S = NaN(size(Y,1),size(Y,1),numel(t));
Sx = sqrt(diag(P0)).*ones(size(X0,1),numel(t));
Omega = eye(size(Q,1)); % what is omega?
% Loop through time input 
for k = 1:numel(t)-1
    % Perform Time Update
    %Xm = Fnl(Xp,U(:,k),[0;0;0]);
    Xm = Xp + 0.1*Fnl(Xp,U(:,k),[0;0;0]); % euler integration 
    
    % debug code, ode45 for more exact solution, too slow
    %sol = ode45(@(~,X) Fnl(X,U(:,k),[0; 0; 0]), [0,0.1], Xp)
    
    Pm = Flin*Pp*Flin' + Omega*Q*Omega;
    
    % Perform Measurement Update
    Ym = Hnl(Xm, U(:,k+1) , [0; 0; 0; 0]);
    K = Pm*Hlin'/(Hlin*Pm*Hlin' + R);
    Xp = Xm + K*(Y(:,k+1) - Ym);
    Pp = (eye(n) - K*Hlin)*Pm;
    
    % Update output Vectors
    Xh(:,k+1) = Xp;
    P(:,:,k+1) = Pp;
    S(:,:,k+1) = Hlin*Pm*Hlin' + R;
    Sx(:,k+1) = sqrt(diag(Pp));
    
    % Measurement
    Yh(:,k+1) = Hlin*Xm + M*U(:,k+1);
end

end
