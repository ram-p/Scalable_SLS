% The actual cvx implementation of the test code in MPGraphCluster.m. This
% imposes SLS constraints by looping over the 'quotient graph', i.e., the
% graph between clusters.

close all
clear
clc
format compact

M = 2;          % Number of models
k = 3;          % Memory
Nc = 4;         % Number of clusters = Number of controllers
T = 10;         % Time horizon

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

X = MP(A,B,C,M,k);

Nseq = M^k;     % Total number of sequences

% Every row is one history
seqs = dec2base(0:Nseq-1,M) - '0' + 1;

% Graph will contain Nseq*M edges
HistoryEdges = struct();
edge = 1;

for node = 1:Nseq
    history = seqs(node,:);
    prefix = history(1:end-1);

    for mode = 1:M
        nextHistory = [mode prefix];

        % Find which row this history is
        nextNode = find(ismember(seqs,nextHistory,'rows'));

        HistoryEdges(edge).from = node;
        HistoryEdges(edge).to   = nextNode;
        HistoryEdges(edge).mode = mode;
        HistoryEdges(edge).history_from = history;
        HistoryEdges(edge).history_to   = nextHistory;

        edge = edge + 1;
    end
end

% Clustering with k-means
rng(1);
[historyCluster, centres] = kmeans(X,Nc);

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

for e = 1:length(HistoryEdges)

    fromStr = sprintf('%d', HistoryEdges(e).history_from);
    toStr   = sprintf('%d', HistoryEdges(e).history_to);

    fprintf('%s (C%d) --(%d)--> %s (C%d)\n', ...
        fromStr, ...
        HistoryEdges(e).fromCluster, ...
        HistoryEdges(e).mode, ...
        toStr, ...
        HistoryEdges(e).toCluster);

end

cvx_begin

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

for t = 1:T+1
    for tau = max(1,t-k+1):t-1
    % for tau = 1:t
    % With memory over measurements, only i-k+1:i. If no, then 1:i.
        for l = 1:Nc
            Phixx{t,tau,l} = Phixxs(:,:,t,tau,l);
            Phixy{t,tau,l} = Phixys(:,:,t,tau,l);
            Phiux{t,tau,l} = Phiuxs(:,:,t,tau,l);
            Phiuy{t,tau,l} = Phiuys(:,:,t,tau,l);
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
    for l = 1:Nc
        Phixx{t+1, t, l} == eye(n)
        Phixy{t+1, t, l} == zeros(n,p)
        Phiux{t+1, t, l} == zeros(m,n)
    end
end

% Affine constraints from SLS imposed iteratively
for t = 1:T
    for tau = max(1,t-k+1):t-1
        for e = 1:length(QuotientEdges)

            i = QuotientEdges(e).from;
            j = QuotientEdges(e).to;
            mode = QuotientEdges(e).mode;

            Phixx{t+1, tau, j} == A{mode}*Phixx{t, tau, i} + B{mode}*Phiux{t, tau, i}
            Phixy{t+1, tau, j} == A{mode}*Phixy{t, tau, i} + B{mode}*Phiuy{t, tau, i}

            Phixx{t, tau+1, j} == Phixx{t, tau, i}*A{mode} + Phixy{t, tau, i}*C{mode}
            Phiux{t, tau+1, j} == Phiux{t, tau, i}*A{mode} + Phiuy{t, tau, i}*C{mode}

        end
    end
end

cvx_end