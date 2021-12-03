function ey = NIS(Y,Yh,S)

% Calculate estimation errors
err = Y - Yh;

% Initialize NEES statistic
N = size(Y,2);
ey = NaN(1,N);

% Loop through each time step
for k = 1:N
    ey(k) = err(:,k)'*inv(S(:,:,k))*err(:,k);
end