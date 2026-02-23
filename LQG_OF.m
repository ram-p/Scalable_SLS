% LQG problem using Scalable SLS and switching signal memory. Continues to
% result in a cost of 1.7321 always. Not fully explainable at the moment.
% Perhaps we should run trajectories and compare between different memory
% sizes.

close all
clear
clc
format compact

M = 2;      % Number of systems
T = 10;     % Time horizon

% System dynamics in discrete-time
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

% Setting up next subwords
P_next = cell(1, lengthP);      % Set of possible next subwords
P_next{1} = {1, 2};   % If p = 1111, p+ = 1111 or 1112.
for i = 2:lengthP-1
    P_next{i} = {i+1};     % Intermediate subwords are just shifted.
end
P_next{end} = {lengthP};       % If p = 2222, p+ = 2222.

% We now use cvx to code the optimization using SLS.
% This is the state-feedback case. Working on output-feedback as well.

cvx_begin quiet
variable Phixxs(n, n, T+1, T+1, lengthP)
variable Phixys(n, m, T+1, T+1, lengthP)
variable Phiuxs(p, n, T+1, T+1, lengthP)
variable Phiuys(p, m, T+1, T+1, lengthP)
variable gam

% Populating cells for ease of use
Phixx = repmat({zeros(n,n)}, T+1, T+1, lengthP);
Phixy = repmat({zeros(n,m)}, T+1, T+1, lengthP);
Phiux = repmat({zeros(p,n)}, T+1, T+1, lengthP);
Phiuy = repmat({zeros(p,m)}, T+1, T+1, lengthP);

for i = 1:T+1
    for j = max(1,i-k+1):i
        for l = 1:lengthP
            Phixx{i,j,l} = Phixxs(:,:,i,j,l);
            Phixy{i,j,l} = Phixys(:,:,i,j,l);
            Phiux{i,j,l} = Phiuxs(:,:,i,j,l);
            Phiuy{i,j,l} = Phiuys(:,:,i,j,l);
        end
    end
end

minimize gam

% Setting up objective
N = zeros(T+1, lengthP);
for i = 1:T+1
for j = 1:T+1
for l = 1:lengthP
    gam >= norm([Phixx{i,j,l}  Phixy{i,j,l};  Phiux{i,j,l}  Phiuy{i,j,l}], "fro");
end
end
end

% Identity constraint
for t = 1:T
    for j = 1:lengthP
        Phixx{t+1, t, j} == eye(n)
        Phixy{t+1, t, j} == zeros(n,m)
        Phiux{t+1, t, j} == zeros(p,n)
    end
end

% Affine constraints from SLS imposed iteratively over subwords.
for t = 1:T
for tau = max(1,t-k+1):t-1
    for i = 1:lengthP
    % Remove the following loop? Since this is an issue with system
    % responses, we're only concerned with how system responses change from
    % t to t+1, and maybe not the fact that p also changes to p_next?
    for j = P_next{i}
        J = j{:};
        Phixx{t+1, tau, J} == A{P{i}(end)}*Phixx{t, tau, i} + B{P{i}(end)}*Phiux{t, tau, i}
        Phixy{t+1, tau, J} == A{P{i}(end)}*Phixy{t, tau, i} + B{P{i}(end)}*Phiuy{t, tau, i}
        Phixx{t+1, tau, J} == Phixx{t, tau, i}*A{P{i}(end)} + Phixy{t, tau, i}*C{P{i}(end)}
        Phiux{t+1, tau, J} == Phiux{t, tau, i}*A{P{i}(end)} + Phiuy{t, tau, i}*C{P{i}(end)}
    end
    end
end
end

cvx_end

% Recovering controller K
K = cell(1, lengthP);
for i = 1:lengthP
    K{i} = cell2mat(Phiuy(:,:,i)) - cell2mat(Phiux(:,:,i))*(cell2mat(Phixx(:,:,i))\cell2mat(Phixy(:,:,i)));
end

% Simulation
x0 = -1+2*rand(n,1);
x = [x0 zeros(n,T)];
u = zeros(p,T);
y = zeros(m*(T+1),1);

tf = 6;
sigma = [ones(1, tf) M*ones(1, T-tf+1)];

for t = 1:T
    w = -1 + 2*randn(n,1);
    v = -1 + 2*randn(m,1);

    sig = [ones(1,k-t) sigma(max(1,t-k+1):t)];             % Last k letters of switching signal
    K_final = K{cellfun(@(x) isequal(x,sig), P)};   % Use controller corresponding to last k letters
    K_final = K_final(:,(n+1):end);

    % y((m*(t-1)+1):m*t) = C{L_fault(t)}*x(:,t) + v;
    vecx1 = x(:);
    vecx = vecx1(1:T*n);
    u(:,t) = K_final(p*(t-1)+1:p*t,:)*vecx;
    x(:,t+1) = A{sigma(t)}*x(:,t) + B{sigma(t)}*u(:,t) + w;
end