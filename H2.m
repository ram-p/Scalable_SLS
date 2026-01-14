function [x, u, K_final, t_fault] = H2(ctr, T)

% Function that takes in a parameter ctr for which controller to use, and
% returns state and input history.
% ctr = 1: fault-tolerant
% ctr = 2: no fault-tolerance
% ctr = 3; memoryless

M = 2;      % Number of systems

% System dynamics % Continuous-time
% A{1} = [-0.9967 0 0.6176; 0 -0.5057 0; -0.0939 0 -0.2127];
% B{1} = [0 -4.2423 4.2423 1.4871; 1.6532 -1.2735 -1.2735 0.0024; 0 -0.2805 0.2805 -0.8823];
% Discrete-time
A{1} = [0.355 0 0.3428; 0 0.6031 0; -0.0521 0 0.7901];
B{1} = [0 -2.72 2.72 0.7376; 1.298 -0.9996 -0.9996 0.0019; 0 -0.1153 0.1153 -0.8362];
n = size(A{1}, 1);
p = size(B{1}, 2);
C{1} = eye(n);
m = size(C{1}, 1);
% Switched dynamics --- drift term added
A{2} = A{1}-1.5*eye(n);
B{2} = B{1};
C{2} = C{1};

% Noise covariances
Sw = eye(n*(T+1));
Sv = eye(m*(T+1));
Cov = blkdiag(Sw, Sv);

% Cost matrices
Q = eye(n);
R = 2*eye(p);
Qblock = kron(eye(T+1), Q);
Rblock = kron(eye(T+1), R);
Cost = blkdiag(sqrtm(Qblock), sqrtm(Rblock));

% lengthL = T+2;              % Number of possible switching signals, including switches at start/end.
% pi_sigma = [0.15 0.07*ones(1,lengthL-2) 0.15]; % Just hard-coded for now.

% Continue to assume non-deterministic switching and not stochastic
% switching, unlike the L-CSS paper.

% Setting up subwords. This is specific to the case of only one switch.
k = 2;                  % Subword length (memory over switching signal)
lengthP = k+1;          % Number of possible subwords (M^k in the most general case)
P = cell(1, lengthP);   % Initializing set of all possible subwords
% lengthP is thus the number of possible controllers after subword-based
% clustering.

% Create set of possible subwords
for i = 1:lengthP
    P{i} = [ones(1, lengthP-i) M*ones(1, i-1)];
end

% Count for prefix constraints - compare number of ones with last signal
for i = 1:lengthP
    ct(i) = nnz(P{i} == P{end});
end

%% Everything below this is incomplete, and must be modified for subword clustering.
% We now use cvx to code the optimization using SLS.

cvx_begin quiet
variable Phixx(n*(T+1), n*(T+1), lengthL) lower_triangular;
variable Phixy(n*(T+1), m*(T+1), lengthL) lower_triangular;
variable Phiux(p*(T+1), n*(T+1), lengthL) lower_triangular;
variable Phiuy(p*(T+1), m*(T+1), lengthL) lower_triangular;

for i = 1:lengthL
    N(i) = square_pos(norm(Cost*[Phixx(:,:,i) Phixy(:,:,i); Phiux(:,:,i) Phiuy(:,:,i)]*Cov, 'fro'));
end

minimize sum(pi_sigma.*N)

for i = 1:lengthL
    % Affine constraints from SLS, also block diagonality.
    [eye(n*(T+1))-Z*Ablock{i} -Z*Bblock{i}]*[Phixx(:,:,i) Phixy(:,:,i); Phiux(:,:,i) Phiuy(:,:,i)] ...
        == [eye(n*(T+1)) zeros(n*(T+1),m*(T+1))]
    [Phixx(:,:,i) Phixy(:,:,i); Phiux(:,:,i) Phiuy(:,:,i)]*[eye(n*(T+1))-Z*Ablock{i}; -Cblock{i}] ...
        == [eye(n*(T+1)); zeros(p*(T+1),n*(T+1))]

    % Prefix constraints due to switching
    % All are equated to the last Phi^sigma since that is how the count
    % is set up.
    Phixx(1:ct(i)*n, 1:ct(i)*n, i) == Phixx(1:ct(i)*n, 1:ct(i)*n, end)
    Phixy(1:ct(i)*n, 1:ct(i)*m, i) == Phixy(1:ct(i)*n, 1:ct(i)*m, end)
    Phiux(1:ct(i)*p, 1:ct(i)*n, i) == Phiux(1:ct(i)*p, 1:ct(i)*n, end)
    Phiuy(1:ct(i)*p, 1:ct(i)*m, i) == Phiuy(1:ct(i)*p, 1:ct(i)*m, end)
end
cvx_end

K = cell(1, lengthL);
for i = 1:lengthL
    K{i} = Phiuy(:,:,i) - Phiux(:,:,i)*(Phixx(:,:,i)\Phixy(:,:,i));
end

x0 = randn(n,1);
x = [x0 zeros(n,T)];
u = zeros(p,T+1);
y = zeros(m*(T+1),1);

% t_fault = randi(T+1);   % Fault time step (in Matlab array index)
t_fault = 5;
L_fault = L{t_fault};   % Faulty switching signal

% Note that for convenience in code, we can define K_FTC ahead of time due
% to prefix constraints.
K_FTC = K{t_fault};             % Fault tolerant controller for that signal
K_NF = K{end};                  % Non-fault tolerant controller
% Controller with no memory
K_diag = zeros(size(K_FTC));
for i = 1:size(K_FTC,1)
    for j = 1:size(K_FTC,2)
        if i == j
            K_diag(i,j) = K_FTC(i,j);
        end
    end
end

K_final = (ctr==1)*K_FTC + (ctr==2)*K_NF + (ctr==3)*K_diag;

for t = 1:T
    w = randn(n,1);
    v = randn(m,1);

    y((m*(t-1)+1):m*t) = C{L_fault(t)}*x(:,t) + v;
    u(:,t) = K_final(p*(t-1)+1:p*t,:)*y;
    % if t > t_fault
    %     u(:,t) = zeros(p,1);
    % end
    x(:,t+1) = A{L_fault(t)}*x(:,t) + B{L_fault(t)}*u(:,t) + w;
end

end