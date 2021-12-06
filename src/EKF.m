% Kalman Filter Code

function [Xh,Yh,P,S,Sx] = EKF(t,Y,U,X0,P0,Xnom,Unom,Fnl,Flin,Hnl,Hlin,Om,Q,R)
n = size(X0,1);

% Initialize Filter
Xp = X0;
Pp = P0;

% Initialize Output Vectors
Xh = X0.*ones(size(X0,1),numel(t));
Yh = NaN(size(Y));
P = P0.*ones(size(X0,1),size(X0,1),numel(t));
S = NaN(size(Y,1),size(Y,1),numel(t));
Sx = sqrt(diag(P0)).*ones(size(X0,1),numel(t));

% Loop through time input 
for k = 1:numel(t)-1
    % Time values
    tk = t(k);
    tp = t(k+1);
    dT = tp - tk;
    
    % -------- Time Update Section --------
    
    % Matrices for time update
    Fk = Flin(dT, tk, Xp, Unom(tk));
    Omk = Om(dT, tk, Xp, Unom(tk));
    Qk = Q(dT, tk);
    
    %Xm = Fnl(Xp,U(:,k),[0;0;0]);
    Xm = Xp + dT*Fnl(Xp,U(:,k),[0;0;0]); % euler integration 
    
    % debug code, ode45 for more exact solution, too slow
%     sol = ode45(@(~,X) Fnl(X,U(:,k),[0; 0; 0]), [0,dT], Xp);
%     Xm = sol.y(:,end);
    
    Pm = Fk*Pp*Fk' + Omk*Qk*Omk';
    
    
    % -------- Measurement Update Section --------
    
    % Matrices for measurement update
    Hk = Hlin(dT, tp, Xm, Unom(tp));
    Rk = R(dT, tp);
    
    % Perform Measurement Update
    Ym = Hnl(Xm, U(:,k+1) , [0; 0; 0; 0]);
    K = Pm*Hk'/(Hk*Pm*Hk' + Rk);
    Xp = Xm + K*(Y(:,k+1) - Ym);
    Pp = (eye(n) - K*Hk)*Pm;
    
    % Update output Vectors
    Xh(:,k+1) = Xp;
    P(:,:,k+1) = Pp;
    S(:,:,k+1) = Hk*Pm*Hk' + Rk;
    Sx(:,k+1) = sqrt(diag(Pp));
    
    % Measurement
    Yh(:,k+1) = Ym;
end

end
