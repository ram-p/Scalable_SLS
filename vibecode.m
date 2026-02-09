n = 3;      % block size
B = 4;      % number of blocks
N = n * B;

cvx_begin
    % Start from a constant zero matrix
    X = zeros(N,N);

    % Fill lower-triangular blocks only
    for i = 1:B
        for j = 1:i
            % Declare block variable
            variable Xij(n,n)

            % Insert into global matrix
            rows = (i-1)*n + (1:n);
            cols = (j-1)*n + (1:n);
            X(rows, cols) = Xij;
        end
    end

    % Example constraints
    % (optional, just to show usage)
    for i = 1:B
        rows = (i-1)*n + (1:n);
        X(rows, rows) >= eye(n);
    end

    minimize( norm(X, 'fro') )
cvx_end