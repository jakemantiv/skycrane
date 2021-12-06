% Function to evaluate NIS statistic on entire block of results
function eym = NIS(Y,Yh,S)

% Size parameters
Nsim = size(Y,3);
N = size(Y,2)-1;

% Calculate estimation errors
err = Y - Yh;

% Initialize NIS statistic
eym = NaN(1,N);
ey = NaN(1,Nsim);

% Loop through each time step
for k = 1:N
    % Loop through each simulation
    for s = 1:Nsim
        ey(s) = err(:,k+1,s)'*inv(S(:,:,k+1,s))*err(:,k+1,s);
    end
    eym(k) = mean(ey);
end