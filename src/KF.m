% Linearized Kalman Filter Code
function [Xh,Yh,P,S,Sx] = KF(t,Y,U,X0,P0,Xnom,Unom,Ynom,F,G,H,M,Om,Q,R)
% Size components for preallocations
p = size(Y,1);
n = size(X0,1);
N = numel(t);
I = eye(n);

% Control perturbation
dU = U - Unom(t);

% Measurement Perturbation
dY = Y - Ynom(t);

% Initialize Filter
dXp = X0 - Xnom(0);
Pp = P0;

% Initialize Output Vectors
dXh = dXp.*ones(n,N);
dYh = NaN(p,N);
P = Pp.*ones(n,n,N);
S = NaN(p,p,N);
Sx = sqrt(diag(Pp)).*ones(n,N);

% Loop through time input
for k = 1:N-1
    % Time values
    tk = t(k);
    tp = t(k+1);
    dT = tp - tk;
    
    % -------- Time Update Section --------
    
    % Matrices for time update
    Fk = F(dT, tk, Xnom(tk), Unom(tk));
    Gk = G(dT, tk, Xnom(tk), Unom(tk));
    Omk = Om(dT, tk, Xnom(tk), Unom(tk));
    Qk = Q(dT, tk);
    
    % Control at time = tk
    dUk = dU(:,k);
    
    % Perform Time Update
    dXm = Fk*dXp + Gk*dUk;
    Pm = Fk*Pp*Fk' + Omk*Qk*Omk';
    
    
    % -------- Measurement Update Section --------
    
    % Matrices for measurement update
    Hk = H(dT, tp, Xnom(tp), Unom(tp));
    Mk = M(dT, tp, Xnom(tp), Unom(tp));
    Rk = R(dT, tp);
    
    % Measurement perturbation at time = tk+1
    dYk = dY(:,k+1);
    
    % Perform Measurement Update
    Kk = Pm*Hk'*inv(Hk*Pm*Hk' + Rk);
    dXp = dXm + Kk*(dYk - Hk*dXm - Mk*dUk);
    Pp = (I - Kk*Hk)*Pm;
    
    % Update output Vectors
    dXh(:,k+1) = dXp;
    dYh(:,k+1) = Hk*dXm + Mk*dUk;
    P(:,:,k+1) = Pp;
    S(:,:,k+1) = Hk*Pm*Hk' + Rk;
    Sx(:,k+1) = sqrt(diag(Pp));
end

% Evaluate actual state for output
Xh = Xnom(t) + dXh;
Yh = Ynom(t) + dYh;
end
