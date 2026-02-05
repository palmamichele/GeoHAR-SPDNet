function Y = vl_logeloss(X, c, dzdy)
% VL_LOGELOSS
% Log-Euclidean Metric/Loss function for SPD matrices
% Supports full Fréchet derivative in backward pass


batchSize = length(c);
n = size(X{1}, 1);  % assuming square SPD matrices

if nargin < 3
    % ------------------------------
    % Forward pass: compute loss
    % ------------------------------
    Y = 0;
    for i = 1:batchSize
        diff = logm(X{i}) - logm(c{i});
        logEuclideanLoss = sum(diff(:).^2) / (n^2);
        Y = Y + logEuclideanLoss;
    end
    Y = Y / batchSize;
else
    % ------------------------------
    % Backward pass: compute gradient w.r.t X
    % ------------------------------
    Y = cell(1, batchSize);  % Initialize cell array for gradients
    for i = 1:batchSize
        Xi = X{i};
        Ci = c{i};

        % Compute difference in tangent space
        Ei = logm(Xi) - logm(Ci);

        % Eigen-decomposition of Xi
        [U, Lambda] = eig(Xi);   % Xi = U*Lambda*U'
        lambda = diag(Lambda);

        % Transform difference into eigenbasis
        G = U' * Ei * U;

        % Compute K matrix for Fréchet derivative
        D = length(lambda);
        K = zeros(D,D);
        for p = 1:D
            for q = 1:D
                if p == q
                    K(p,q) = 1 / lambda(p);
                else
                    K(p,q) = (log(lambda(p)) - log(lambda(q))) / (lambda(p) - lambda(q));
                end
            end
        end

        % Apply Fréchet derivative in eigenbasis
        gradXi = U * (K .* G) * U';

        % Scale by 2/n^2
        gradXi = (2/n^2) * gradXi;

        % Multiply by upstream gradient dzdy
        Y{i} = gradXi * dzdy;
    end
end

end
