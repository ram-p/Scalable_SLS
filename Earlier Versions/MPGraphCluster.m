% A test code with specific M, k, Nc.

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
    prefix = history(2:end);

    for mode = 1:M
        nextHistory = [prefix mode];

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
        historyCluster(HistoryEdges(e).from), ...
        HistoryEdges(e).mode, ...
        toStr, ...
        historyCluster(HistoryEdges(e).to));

end

for t = 1:T

    fprintf('\n========== t = %d ==========\n', t);

    for e = 1:length(QuotientEdges)

        i = QuotientEdges(e).from;
        j = QuotientEdges(e).to;
        mode = QuotientEdges(e).mode;

        fprintf('C%d --(%d)--> C%d\n', i, mode, j);

    end

end

% For nondeterministic transitions.
EdgeTable = zeros(length(QuotientEdges),3);

for e = 1:length(QuotientEdges)

    EdgeTable(e,:) = [ ...
        QuotientEdges(e).from,...
        QuotientEdges(e).mode,...
        QuotientEdges(e).to];

end

% EdgeTable = unique(EdgeTable,'rows');

fprintf('\nAre there nondeterminstic transitions?\n')

for c = 1:Nc
    for mode = 1:M

        % Find all quotient edges with this source cluster and mode
        rows = EdgeTable( ...
            EdgeTable(:,1)==c & EdgeTable(:,2)==mode, :);

        if isempty(rows)
            fprintf('Cluster %d --(%d)--> no transition\n', c, mode);
            continue
        end

        % Unique successor clusters
        successors = unique(rows(:,3));

        fprintf('Cluster %d -- mode %d --> ', c, mode);
        fprintf('C%d ', successors);
        fprintf('\n');

    end
end