% Unscented Kalman Filter
function [Xh,Yh,P,S,Sx] = UKF(t,Y,U,X0,P0,Fnl,Hnl,Hlin,Om,Q,R)
n = size(X0,1);
p = size(Y,1);

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
    Omk = Om(dT, tk, Xp, U(:,k));
    Qk = Q(dT, tk);
    
    % Sigma-Point Selection
    Sp = chol(n*Pp,'upper');
    Xsp = zeros(n,2*n);
    for i = 1:n
        Xsm = Xp + Sp(i,:)';
        Xsp(:,i) = Xsm + dT*Fnl(Xsm,U(:,k),[0;0;0]); % euler integration
        
        Xsm = Xp - Sp(i,:)';
        Xsp(:,n+i) = Xsm + dT*Fnl(Xsm,U(:,k),[0;0;0]); % euler integration
    end
    Xm = 1/(2*n)*sum(Xsp,2);
    
    Ps = zeros(n);
    for i = 1:(2*n)
        Ps = Ps + (Xsp(:,i) - Xm)*(Xsp(:,i) - Xm)';
    end
    Pm = 1/(2*n)*Ps + Omk*Qk*Omk';
    
    
    % -------- Measurement Update Section --------
    
    % Matrices for measurement update
    Hk = Hlin(dT, tp, Xm, U(:,k+1));
    Rk = R(dT, tp);
    
    % Sigma-Point Selection
    Sp = chol(n*Pm,'upper');
    Ysp = zeros(p,2*n);
    Xsp = zeros(n,2*n);
    for i = 1:n
        Xsp(:,i) = Xm + Sp(i,:)';
        Ysp(:,i) = Hnl(Xsp(:,i),U(:,k+1),zeros(p,1));
        
        Xsp(:,n+i) = Xm - Sp(i,:)';
        Ysp(:,n+i) = Hnl(Xsp(:,n+i),U(:,k+1),zeros(p,1));
    end
    
    % Perform Measurement Update
    Ym = 1/(2*n)*sum(Ysp,2);
    
    Pys = zeros(p);
    Pxys = zeros(n,p);
    for i = 1:(2*n)
        Pys = Pys + (Ysp(:,i) - Ym)*(Ysp(:,i) - Ym)';
        Pxys = Pxys + (Xsp(:,i) - Xm)*(Ysp(:,i) - Ym)';
    end
    Py = 1/(2*n)*Pys + Rk;
    Pxy = 1/(2*n)*Pxys;
    
    K = Pxy*inv(Py);
    Xp = Xm + K*(Y(:,k+1) - Ym);
    Pp = Pm - K*Py*K';
    
    % Update output Vectors
    Xh(:,k+1) = Xp;
    P(:,:,k+1) = Pp;
    S(:,:,k+1) = Hk*Pm*Hk' + Rk;
    Sx(:,k+1) = sqrt(diag(Pp));
    
    % Measurement
    Yh(:,k+1) = Ym;
end

end