% Truth Model

function [X,Y,U] = truthModel(t,F,H,X0,Xnom,Unom,control)

% Load data matrix
data = load(['..',filesep,'lib', filesep, 'skycrane_finalproj_KFdata.mat']);

% Feedback control gain matrix
Klin = data.Klin;

% True Noise Matrices
Q = data.Qtrue;
R = data.Rtrue;

% Feedback control law
if control
    Ucl = @(X,Xn,Un) Un - Klin*(X - Xn);
else
    Ucl = @(X,Xn,Un) Un;
end

% Initialize Outputs
N = numel(t);
n = numel(X0);
m = size(Unom,1);
p = 4;
X = X0.*ones(n,N);
U = NaN(m,N);

% Nominal Matrices
if size(Xnom,2) == 1
    Xnom = Xnom.*ones(n,N);
end
if size(Unom,2) == 1
    Unom = Unom.*ones(m,N);
end

% Simulation loop
for k = 1:N-1
    % Time step
    dT = t(k+1) - t(k);
    
    % Get Control based on prior step
    U(:,k) = Ucl(X(:,k),Xnom(:,k),Unom(:,k));
    
    % Process Noise Sample
    w = awgn_draw(zeros(size(Q,1),1),Q);
    
    % Evaluate Time update 
    sol = ode45(@(~,X) F(X,U(:,k),w), [0,dT], X(:,k));
    X(:,k+1) = sol.y(:,end);
end

% Measurement Noise Sample
v = awgn_draw(zeros(p,N),R);

% Measurements
Y = NaN(p,N);
for k = 1:N
    Y(:,k) = H(X(:,k),U(:,k),v(:,k));
end
