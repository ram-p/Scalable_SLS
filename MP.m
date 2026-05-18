% Function to form matrix of all possible Markov parameters (M^k), where M
% is number of modes and k is memory. Once this matrix is formed, it is
% easy to perform clustering with Nc clusters and assign the most recent
% subword to a given cluster.

function X = MP(A, B, C, M, k)

    m = size(B{1},2);
    p = size(C{1},1);

    % Total number of sequences
    Nseq = M^k;

    % Generate all sequences (each row is one sequence)
    % Entries are in {1,...,M}
    seqs = dec2base(0:Nseq-1, M) - '0' + 1;  % size: Nseq x k

    % Preallocate data matrix
    X = zeros(Nseq, p*m*k);

    for s = 1:Nseq
        sigma = seqs(s,:);  % length-k sequence

        % Store Markov parameters for this sequence
        vec = zeros(1, p*m*k);

        for ell = 1:k
            % ell = 1 corresponds to p_t
            % ell = k corresponds to p_{t-k+1}

            idx_start = (ell-1)*p*m + 1;
            idx_end   = ell*p*m;

            % Build A product
            Aprod = eye(size(A{1}));
            for j = 2:ell
                Aprod = A{sigma(j-1)} * Aprod;
            end

            % Select C and B
            Cmat = C{sigma(1)};  % always sigma_t
            Bmat = B{sigma(ell)};

            % Compute Markov parameter
            P = Cmat * Aprod * Bmat;

            % Vectorize (column-wise)
            vec(idx_start:idx_end) = P(:)';
        end

        % Store in data matrix
        X(s,:) = vec;
    end
end