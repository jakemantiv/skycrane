% Kalman Filter Code

function [Xh,Yh,P,S,Sx] = KF(t,Y,U,X0,P0,F,G,H,M,Q,R)
[p,n] = size(H);

% Remove control feedthrough inputs from measurements
Y = Y - M*U;

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
    % Perform Time Update
    Xm = F*Xp + G*U(:,k);
    Pm = F*Pp*F' + Q;
    
    % Perform Measurement Update
    K = Pm*H'*inv(H*Pm*H' + R);
    Xp = Xm + K*(Y(:,k+1) - H*Xm);
    Pp = (eye(n) - K*H)*Pm;
    
    % Update output Vectors
    Xh(:,k+1) = Xp;
    P(:,:,k+1) = Pp;
    S(:,:,k+1) = H*Pm*H' + R;
    Sx(:,k+1) = sqrt(diag(Pp));
    
    % Measurement
    Yh(:,k+1) = H*Xm + M*U(:,k+1);
end

end
