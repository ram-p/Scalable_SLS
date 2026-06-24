% LQG problem using Scalable Output-Feedback SLS, with memory over
% measurements and switching signal, as well as Markov parameter
% clustering. Run gamvk.m to see impacts of memory on cost and runtime.

% NOTE ——— Code does not currently work. Switching to next cluster is
% incorporated, but program still results in infeasibility.

function [gam, cvx_cputime] = LQG_MP(k, T, Nc)

% k denotes memory over both switching signal and measurements, T denotes
% the full time horizon, and Nc denotes the number of clusters.

M = 2;      % Number of systems

% System dynamics in discrete-time: double integrator, unit sampling time
A{1} = blkdiag([1 1; 0 1], [1 1; 0 1]);
B{1} = [0 0; 1 0; 0 0; 0 1];
n = size(A{1}, 1);
m = size(B{1}, 2);
C{1} = [1 0 0 0; 0 0 1 0];
p = size(C{1}, 1);
% Switched dynamics --- changed integrator gains 
A{2} = blkdiag([1 4; 0 1], [1 4; 0 1]);
B{2} = [0 0; 4 0; 0 0; 0 4];
C{2} = C{1};

% % System dynamics in discrete-time: ADMIRE fighter-jet model
% A{1} = [0.355 0 0.3428; 0 0.6031 0; -0.0521 0 0.7901];
% B{1} = [0 -2.72 2.72 0.7376; 1.298 -0.9996 -0.9996 0.0019; 0 -0.1153 0.1153 -0.8362];
% n = size(A{1}, 1);
% m = size(B{1}, 2);
% C{1} = eye(n);
% p = size(C{1}, 1);
% % Switched dynamics --- partial sensor failure
% A{2} = A{1};
% B{2} = B{1};
% C{2} = diag([1 0 0]);

% Noise bounds
wbar = 1;
vbar = 1;

% Matrix of Markov parameters
X = MP(A,B,C,M,k);
Nseq = M^k;      % Number of possible sequences of Markov parameters

% Setting up clustering
% Nc = 4;                           % Number of clusters = Number of controllers
[idx, centres] = kmeans(X, Nc);     % k-means clustering on Markov parameters, with Nc clusters.

% % Create set of possible subwords
% for i = 1:lengthP
%     P{i} = [ones(1, lengthP-i) M*ones(1, i-1)];
% end

% % Setting up next subwords
% P_next = cell(1, lengthP);      % Set of possible next subwords
% P_next{1} = {1, 2};   % If p = 1111, p+ = 1111 or 1112.
% for i = 2:lengthP-1
%     P_next{i} = {i+1};     % Intermediate subwords are just shifted.
% end
% P_next{end} = {lengthP};       % If p = 2222, p+ = 2222.

% All possible sequences of modes
seqs = dec2base(0:Nseq-1, M) - '0' + 1;

% Using cvx to code the optimization using SLS.
cvx_begin
% cvx_solver sedumi
variable Phixxs(n, n, T+1, T+1, Nc)
variable Phixys(n, p, T+1, T+1, Nc)
variable Phiuxs(m, n, T+1, T+1, Nc)
variable Phiuys(m, p, T+1, T+1, Nc)
variable gam

% Populating cells for ease of use
Phixx = repmat({zeros(n,n)}, T+1, T+1, Nc);
Phixy = repmat({zeros(n,p)}, T+1, T+1, Nc);
Phiux = repmat({zeros(m,n)}, T+1, T+1, Nc);
Phiuy = repmat({zeros(m,p)}, T+1, T+1, Nc);

for i = 1:T+1
    for j = max(1,i-k+1):i
    % for j = 1:i
        for l = 1:Nc
            Phixx{i,j,l} = Phixxs(:,:,i,j,l);
            Phixy{i,j,l} = Phixys(:,:,i,j,l);
            Phiux{i,j,l} = Phiuxs(:,:,i,j,l);
            Phiuy{i,j,l} = Phiuys(:,:,i,j,l);
        end
    end
end

minimize gam

% Setting up objective
for l = 1:Nc
    gam >= norm([Phixx{:,:,l}  Phixy{:,:,l}; Phiux{:,:,l}  Phiuy{:,:,l}], "fro")
end

% Identity constraint
for t = 1:T
    for j = 1:Nc
        Phixx{t+1, t, j} == eye(n)
        Phixy{t+1, t, j} == zeros(n,p)
        Phiux{t+1, t, j} == zeros(m,n)
    end
end

% Affine constraints from SLS imposed iteratively.
for i = 1:Nc
% Finding the sigmas that produced a path in the given cluster i
    % Finding the rows of X that correspond to the current cluster
    rowsX = ismember(idx, i, 'rows');
    % Finding all sigmas corresponding to those rows of X
    sig = seqs(rowsX,:);
    % Finding the unique elements of each sig.
    sig_u = unique(sig, 'rows');
    % Evaluating all possible "next" sigmas to incorporate inter-cluster
    % dependency.
    sig_next_pre = sig_u(:, 2:end);     % Prefix for next sigma uses 2:end most recent letters.
    sig_next = [];
    for l = 1:size(sig_next_pre, 1)
        for model = 1:M
            % All possible next sigmas by appending each letter to each
            % prefix.
            sig_next = [sig_next; [sig_next_pre(l,:) model]];
        end
    end
    sig_next = unique(sig_next, 'rows');

    % Finding all possible next clusters
    next_clusters = unique(idx(ismember(seqs, sig_next, 'rows')), 'rows')';

    for t = 1:T
    for tau = max(1,t-k+1):t-1
        for j = next_clusters

        % SLS constraints: remember: len(sig) families of constraints since
        % the constraints are written for EACH possible switching signal
        % found above. TO CONSIDER: Next cluster on the left-hand side, and
        % how to do that.
        for m_ind = 1:size(sig_u, 1)
            mode = sig_u(m_ind, t-tau+1);
            % Indicates which mode the system is in. Note t-tau is between
            % 1 and k-1, since tau can take a maximum of k-1 values: t-k+1 <= tau <= t-1.
            Phixx{t+1, tau, j} == A{mode}*Phixx{t, tau, i} + B{mode}*Phiux{t, tau, i}
            Phixy{t+1, tau, j} == A{mode}*Phixy{t, tau, i} + B{mode}*Phiuy{t, tau, i}
            Phixx{t+1, tau, j} == Phixx{t, tau, i}*A{mode} + Phixy{t, tau, i}*C{mode}
            Phiux{t+1, tau, j} == Phiux{t, tau, i}*A{mode} + Phiuy{t, tau, i}*C{mode}
        end
        end
    end
    end
end

cvx_end

end

% % Recovering controller K
% K = cell(1, Nc);
% for i = 1:Nc
%     K{i} = cell2mat(Phiuy(:,:,i)) - cell2mat(Phiux(:,:,i))*(cell2mat(Phixx(:,:,i))\cell2mat(Phixy(:,:,i)));
% end
% 
% % Simulation: initializing
% x0 = -1+2*rand(n,1);
% x = [x0 zeros(n,T)];
% u = zeros(m,T);
% y = zeros(p*(T+1),1);
% 
% % Time of fault tf and switching signal
% tf = 6;
% sigma = [ones(1, tf) M*ones(1, T-tf+1)];
% 
% for t = 1:T
%     % Noise vectors
%     w = -1 + 2*randn(n,1);
%     v = -1 + 2*randn(p,1);
% 
%     sig = [ones(1,k-t) sigma(max(1,t-k+1):t)];          % Last k letters of switching signal
%     K_final = K{cellfun(@(x) isequal(x,sig), P)};       % Use controller corresponding to last k letters
% 
%     y((p*(t-1)+1):p*t) = C{sigma(t)}*x(:,t) + v;        % y = Cx + v
%     u(:,t) = K_final(m*(t-1)+1:m*t,:)*y;                % u = Ky
%     x(:,t+1) = A{sigma(t)}*x(:,t) + B{sigma(t)}*u(:,t) + w;
% end