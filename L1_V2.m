% The currently working version. An LQG (+safety?) version is currently in
% development. L1_V1.m and vibecode.m are not useful.

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
k = 8;                  % Subword length (memory over switching signal)
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
P_next{1} = {1, 2};   % If p = 11, p+ = 11 or 12.
for i = 2:lengthP-1
    P_next{i} = {i+1};     % Intermediate subwords are just shifted.
end
P_next{end} = {lengthP};       % If p = 22, p+ = 22.

% We now use cvx to code the optimization using SLS.
% This is the state-feedback case. Working on output-feedback as well.

cvx_begin quiet
variable Phixs(n, n, T+1, T+1, lengthP)
variable Phius(p, n, T+1, T+1, lengthP)
variable gam

% Populating cells for ease of use
% Phix = cell(T+1, T+1, lengthP);
% Phiu = cell(T+1, T+1, lengthP);
Phix = repmat({zeros(n,n)}, T+1, T+1, lengthP);
Phiu = repmat({zeros(p,n)}, T+1, T+1, lengthP);
Phixnorm = cell(T+1, lengthP);      % For matrix norm

% for i = 1:T+1
%     for j = 1:i
%         for l = 1:lengthP
%             expression Phixs(n,n)
%             expression Phius(p,n)
%             Phix{i,j,l} = Phixs;
%             Phiu{i,j,l} = Phius;
%             Phixnorm{i,l} = [Phix{i,:,l}];
%         end
%     end
% end

for i = 1:T+1
    for j = max(1,i-k+1):i
        for l = 1:lengthP
            Phix{i,j,l} = Phixs(:,:,i,j,l);
            Phiu{i,j,l} = Phius(:,:,i,j,l);
            Phixnorm{i,l} = [Phix{i,:,l}];
        end
    end
end

minimize gam

% Setting up objective: is this correct under subword memory?
N = zeros(T+1, lengthP);
for i = 1:T+1
for l = 1:lengthP
    gam >= 2*norm(Phixnorm{i,l}, Inf);
end
end

% Identity constraint
for t = 1:T
    for j = 1:lengthP
        Phix{t+1, t, j} == eye(n)
    end
end

% Affine constraints from SLS imposed iteratively over subwords.
for t = 1:T
for tau = max(1,t-k+1):t-1
    for i = 1:lengthP
    for j = P_next{i}
        J = j{:};
        Phix{t+1, tau, J} == A{P{i}(end)}*Phix{t, tau, i} + B{P{i}(end)}*Phiu{t, tau, i}
    end
    end
end
end

cvx_end

% Recovering controller K
K = cell(1, lengthP);
for i = 1:lengthP
    K{i} = cell2mat(Phiu(:,:,i))/cell2mat(Phix(:,:,i));
end

% Simulation
x0 = -1+2*rand(n,1);
x = [x0 zeros(n,T)];
u = zeros(p,T);
% y = zeros(m*(T+1),1);

tf = 6;
% t_fault = tf+1;         % Fault time step (in Matlab array index)
sigma = [ones(1, tf) M*ones(1, T-tf+1)];

for t = 1:T
    w = -1 + 2*rand(n,1);
    v = -1 + 2*rand(m,1);

    sig = [ones(1,k-t) sigma(max(1,t-k+1):t)];             % Last k letters of switching signal
    K_final = K{cellfun(@(x) isequal(x,sig), P)};   % Use controller corresponding to last k letters
    K_final = K_final(:,(n+1):end);

    % y((m*(t-1)+1):m*t) = C{L_fault(t)}*x(:,t) + v;
    vecx1 = x(:);
    vecx = vecx1(1:T*n);
    u(:,t) = K_final(p*(t-1)+1:p*t,:)*vecx;
    x(:,t+1) = A{sigma(t)}*x(:,t) + B{sigma(t)}*u(:,t) + w;
end