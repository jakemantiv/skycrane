function ey = NIS(Y,Yh,S)

% Calculate estimation errors
err = Y - Yh;

% Initialize NEES statistic
N = size(Y,2);
ey = NaN(1,N-1);

% Loop through each time step
for k = 1:N-1
    ey(k) = err(:,k+1)'*inv(S(:,:,k+1))*err(:,k+1);
end