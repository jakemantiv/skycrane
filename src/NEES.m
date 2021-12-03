function ex = NEES(X,Xh,P)

% Calculate estimation errors
err = X - Xh;

% Initialize NEES statistic
N = size(X,2)-1;
ex = NaN(1,N);

% Loop through each time step
for k = 1:N
    ex(k) = err(:,k+1)'*inv(P(:,:,k+1))*err(:,k+1);
end