% Truth Model

function [X,Y,U] = truthModel(t,F,H,X0,Xnom,Unom,control)

% Load data matrix
data = load(['..',filesep,'lib', filesep, 'skycrane_finalproj_KFdata.mat']);

% True Noise Parameters
Q = data.Qtrue;
R = data.Rtrue;

% Feedback control gain matrix
Klin = data.Klin;

% Feedback control law
if control
    Ucl = @(t,X,U) U - Klin*(X - Xnom(t));
else
    Ucl = @(t,X,U) U;
end

% Initialize Outputs
N = numel(t);
n = numel(X0);
m = size(Unom(0),1);
p = size(R,1);
X = X0.*ones(n,N);
U = NaN(m,N);

% Simulation loop
for k = 1:N-1
    % Time step
    dT = t(k+1) - t(k);
    
    % Get Control based on prior step
    U(:,k) = Ucl(t(k),X(:,k),Unom(t(k)));
    
    % Process Noise Sample
    w = awgn_draw(zeros(size(Q,1),1),Q);
    
    % Evaluate Time update 
    sol = ode45(@(~,X) F(X,U(:,k),w), [0,dT], X(:,k));
    X(:,k+1) = sol.y(:,end);
end
U(:,k+1) = Ucl(t(k+1),X(:,k+1),U(:,k+1));

% Measurement Noise Sample
v = awgn_draw(zeros(p,N-1),R);

% Measurements
Y = NaN(p,N);
for k = 1:N-1
    Y(:,k+1) = H(X(:,k+1),U(:,k+1),v(:,k));
end
