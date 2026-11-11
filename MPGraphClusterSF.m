% The implementation of scalable, state-feedback SLS through clustering and
% limited memory. SLS constraints are imposed by looping over the 'quotient
% graph', i.e., the graph between clusters.

close all
clear
clc
format compact

T = 10;         % Time horizon
M = 2;          % Number of models
Nc = 4;         % Number of clusters = Number of controllers
k = 3;          % Memory

% % System dynamics in discrete-time: double integrator, unit sampling time
% A{1} = blkdiag([1 1; 0 1], [1 1; 0 1]);
% B{1} = [0 0; 1 0; 0 0; 0 1];
% n = size(A{1}, 1);
% m = size(B{1}, 2);
% C{1} = [1 0 0 0; 0 0 1 0];
% p = size(C{1}, 1);
% % Switched dynamics --- changed integrator gains 
% A{2} = blkdiag([1 4; 0 1], [1 4; 0 1]);
% % A{2} = A{1};
% B{2} = [0 0; 4 0; 0 0; 0 4];
% % B{2} = B{1};
% C{2} = C{1};

% System dynamics in discrete-time: ADMIRE fighter-jet model
A{1} = [0.355 0 0.3428; 0 0.6031 0; -0.0521 0 0.7901];
B{1} = [0 -2.72 2.72 0.7376; 1.298 -0.9996 -0.9996 0.0019; 0 -0.1153 0.1153 -0.8362];
n = size(A{1}, 1);
m = size(B{1}, 2);
C{1} = eye(n);
p = size(C{1}, 1);
% Switched dynamics --- actuator failure
A{2} = A{1};
B{2} = [0 0 2.72 0.7376; 1.298 0 -0.9996 0.0019; 0 0 0.1153 -0.8362];
% C{2} = diag([1 0 0]);
C{2} = C{1};

X = MP(A,B,C,M,k);

Nseq = M^k;     % Total number of sequences
% Nc = Nseq;

% Every row is one history
seqs = dec2base(0:Nseq-1,M) - '0' + 1;

% Graph will contain Nseq*M edges
HistoryEdges = struct();
edge = 1;

for node = 1:Nseq
    history = seqs(node,:);
    prefix = history(2:end);

    for mode = 1:M
        nextHistory = [prefix mode];

        % Find which row this history is
        nextNode = find(ismember(seqs,nextHistory,'rows'));

        HistoryEdges(edge).from = node;
        HistoryEdges(edge).to   = nextNode;
        HistoryEdges(edge).mode = history(end);
        HistoryEdges(edge).history_from = history;
        HistoryEdges(edge).history_to   = nextHistory;

        edge = edge + 1;
    end
end

% Clustering with k-means
rng(1);
% historyCluster = (1:Nseq)';
[historyCluster, centres] = kmeans(X,Nc);
% Each row of X corresponds to the particular row of seqs, using MP.m.
% Then, each row of historyCluster does the same, and corresponds to the
% particular row in seqs!

QuotientEdges = struct();
cedge = 1;

for e = 1:length(HistoryEdges)

    % History indices
    fromHist = HistoryEdges(e).from;
    toHist   = HistoryEdges(e).to;

    % Corresponding clusters
    fromCluster = historyCluster(fromHist);
    toCluster   = historyCluster(toHist);

    % Store cluster transition
    QuotientEdges(cedge).from = fromCluster;
    QuotientEdges(cedge).to   = toCluster;
    QuotientEdges(cedge).mode = HistoryEdges(e).mode;

    HistoryEdges(e).fromCluster = fromCluster;
    HistoryEdges(e).toCluster   = toCluster;

    % Store the history edge that generated it
    QuotientEdges(cedge).history_from_index = fromHist;
    QuotientEdges(cedge).history_to_index   = toHist;
    QuotientEdges(cedge).history_from = seqs(fromHist,:);
    QuotientEdges(cedge).history_to   = seqs(toHist,:);

    cedge = cedge + 1;

end

for e = 1:length(QuotientEdges)

    fromStr = sprintf('%d', QuotientEdges(e).history_from);
    toStr   = sprintf('%d', QuotientEdges(e).history_to);

    fprintf('%s (C%d) --(%d)--> %s (C%d)\n', ...
        fromStr, ...
        QuotientEdges(e).from, ...
        QuotientEdges(e).mode, ...
        toStr, ...
        QuotientEdges(e).to);

end

cvx_begin

variable Phixs(n, n, T+1, T+1, Nc)
variable Phius(m, n, T+1, T+1, Nc)
variable gam

% Populating cells for ease of use
Phix = repmat({zeros(n,n)}, T+1, T+1, Nc);
Phiu = repmat({zeros(m,n)}, T+1, T+1, Nc);

for tau = 1:T
    for t = tau+1:T+1
    % for tau = 1:t
    % With memory over measurements, only i-k+1:i. If no, then 1:i.
        for l = 1:Nc
            Phix{t,tau,l} = Phixs(:,:,t,tau,l);
            Phiu{t,tau,l} = Phius(:,:,t,tau,l);
        end
    end
end

minimize gam

% Setting up objective
for l = 1:Nc
    gam >= norm([Phix{:,:,l};  Phiu{:,:,l}], "fro")
end

% Identity constraint
for t = 1:T
    for l = 1:Nc
        Phix{t+1, t, l} == eye(n)
    end
end

% Affine constraints from SLS imposed iteratively
for tau = 1:T
    for t = tau+1:T
        for e = 1:length(QuotientEdges)

            i = QuotientEdges(e).from;
            j = QuotientEdges(e).to;
            mode = QuotientEdges(e).mode;

            Phix{t+1, tau, j} == A{mode}*Phix{t, tau, i} + B{mode}*Phiu{t, tau, i}

        end
    end
end

cvx_end