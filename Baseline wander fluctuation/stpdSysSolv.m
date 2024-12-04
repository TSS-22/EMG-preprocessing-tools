function b = stpdSysSolv(A, b)
% Based on algorithm 4.3.6, page 181
% GOLUB, Gene H. et VAN LOAN, Charles F. Matrix computations. JHU press, 2013.
% Requires 8n flops
%
% Solve the system Ax=b (useful to avoid computing the inverse when dividing matrices)
% Return the solution x (replace b to limit memory usage)


% Check for the minimum size and that the arrays are not empty
n = length(alpha);

if n >= 3 && length(beta) >= 2 && length(b) == n
    % alpha is the diagonal values
    % beta is the supra-diagonal values
    alpha = diag(A);
    beta = diag(A,1);

    % Forward elimination
    for i = 1:n-1
        temp = beta(i);
        beta(i) = temp / alpha(i);
        alpha(i+1) = alpha(i+1) - temp * beta(i);
    end

    % Modify b during elimination
    for i = 1:n-1
        b(i+1) = b(i+1) - beta(i) * b(i);
    end

    % Back substitution
    b(n) = b(n) / alpha(n);
    for i = n-1:-1:1
        b(i) = (b(i) - beta(i) * b(i+1)) / alpha(i);
    end
end

end
