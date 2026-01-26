close all
clear
clc
format compact

% Function that takes in a parameter ctr for which controller to use, and
% returns state and input history.
% ctr = 1: fault-tolerant
% ctr = 2: no fault-tolerance
% ctr = 3; memoryless

M = 2;      % Number of systems
T = 10;     % Time horizon

% System dynamics
% Continuous-time
% A{1} = [-0.9967 0 0.6176; 0 -0.5057 0; -0.0939 0 -0.2127];
% B{1} = [0 -4.2423 4.2423 1.4871; 1.6532 -1.2735 -1.2735 0.0024; 0 -0.2805 0.2805 -0.8823];
% Discrete-time
A{1} = [0.355 0 0.3428; 0 0.6031 0; -0.0521 0 0.7901];
B{1} = [0 -2.72 2.72 0.7376; 1.298 -0.9996 -0.9996 0.0019; 0 -0.1153 0.1153 -0.8362];
n = size(A{1}, 1);
p = size(B{1}, 2);
C{1} = eye(n);
m = size(C{1}, 1);
% Switched dynamics --- partial sensor failure
A{2} = A{1};
B{2} = B{1};
C{2} = diag([1 0 0]);

% Noise bounds
wbar = 1;
vbar = 1;

% Setting up subwords. This is specific to the case of only one switch.
k = 4;                  % Subword length (memory over switching signal)
lengthP = k+1;          % Number of possible subwords (M^k in the most general case)
P = cell(1, lengthP);   % Initializing set of all possible subwords
% lengthP is thus the number of possible controllers after subword-based
% clustering.

% Create set of possible subwords
for i = 1:lengthP
    P{i} = [ones(1, lengthP-i) M*ones(1, i-1)];
end

% We now use cvx to code the optimization using SLS.
% This is the state-feedback case. Working on output-feedback as well.

cvx_begin quiet
% variable Phixx(n*(T+1), n*(T+1), lengthP) lower_triangular;
% variable Phixy(n*(T+1), m*(T+1), lengthP) lower_triangular;
% variable Phiux(p*(T+1), n*(T+1), lengthP) lower_triangular;
% variable Phiuy(p*(T+1), m*(T+1), lengthP) lower_triangular;
variable Phix(n*(T+1), n*(T+1), lengthP)
variable Phiu(p*(T+1), n*(T+1), lengthP)
variable gam

% Setting up objective: is this correct under subword memory?
for i = 1:lengthP
    N(i) = norm(Phix(:,:,i), Inf);
end

% Setting up next subwords
P_next = cell(1, lengthP);      % Set of possible next subwords
P_next{1} = {1, 2};   % If p = 11, p+ = 11 or 12.
for i = 2:lengthP-1
    P_next{i} = {i+1};     % Intermediate subwords are just shifted.
end
P_next{end} = {lengthP};       % If p = 22, p+ = 22.

minimize gam
gam >= max(N)

% Identity constraint
for t = 1:T-1
    for j = 1:lengthP
        Phix(((t+1)*n+1):(t+2)*n, (t*n+1):(t+1)*n, j) == eye(n)
    end
end

% Affine constraints from SLS imposed iteratively over subwords.
for tau = 0:T
for t = tau:T-1
    for i = 1:lengthP
    for j = P_next{i}
        j = j{:};
        Phix(((t+1)*n+1):(t+2)*n, (tau*n+1):(tau+1)*n, j) == A{P{i}(end)}*Phix((t*n+1):(t+1)*n, (tau*n+1):(tau+1)*n, i) ...
          + B{P{i}(end)}*Phiu((t*p+1):(t+1)*p, (tau*n+1):(tau+1)*n, i)
    end
    end
end
end

cvx_end

K = cell(1, lengthP);
for i = 1:lengthP
    K{i} = Phiu(:,:,i)/Phix(:,:,i);
end

%% Below this is incomplete. Simulation to be modified.
x0 = -1+2*rand(n,1);
x = [x0 zeros(n,T)];
u = zeros(p,T);
% y = zeros(m*(T+1),1);

[~,tf] = max(N);
t_fault = tf+1;         % Fault time step (in Matlab array index)
sigma = [ones(1, tf) M*ones(1, T-tf+1)];

for t = 1:T
    w = 2*randi(2,n,1)-3;
    v = 2*randi(2,m,1)-3;

    sig = [ones(1,k-t) sigma(t-k+1:t)];             % Last k letters of switching signal
    K_final = K{cellfun(@(x) isequal(x,sig), P)};   % Use controller corresponding to last k letters

    % y((m*(t-1)+1):m*t) = C{L_fault(t)}*x(:,t) + v;
    u(:,t) = K_final(p*(t-1)+1:p*t,:)*vec(x);
    x(:,t+1) = A{L_fault(t)}*x(:,t) + B{L_fault(t)}*u(:,t) + w;
end