function yk = awgn_draw(my,R)
% Fuction to sample Gaussian random vector
[p,N] = size(my);
Sv = chol(R,'lower');
qk = randn(p,N);
yk = my + Sv*qk;
end