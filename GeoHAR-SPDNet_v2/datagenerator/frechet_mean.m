function M = frechet_mean(varargin)
    % Input:
    % varargin - SPD matrices passed as individual arguments
    % varargin{end} can contain the computation method as a string: 'procrustes' or 'log-euclidean'
    % Output:
    % M - the Fréchet mean of the SPD matrices
    

    

    % Default method
    if ischar(varargin{end}) || isstring(varargin{end})
        method = varargin{end}; % Extract the method
        varargin = varargin(1:end-1); % Exclude the method from matrices
    else
        method = 'log-euclidean'; % Default method
    end



    
    n = numel(varargin); % Number of matrices passed
    d = size(varargin{1}, 1); % Dimensionality of the SPD matrices



    M=zeros(d);
    switch lower(method)
        case 'log-euclidean'
            % Compute mean in the Log-Euclidean metric
            sumLog = zeros(d);
            for i = 1:n
                sumLog = sumLog + logm(makespd(varargin{i}));
            end
            M = expm(sumLog / n);

        case 'procrustes'
            % Compute mean using the Procrustes method matrix sqrt based
            tol = 1e-6; % Convergence tolerance
            maxIter = 100; % Maximum iterations
            L = cell(1, n);
            for i=1:n
                Si = makespd(varargin{i}); % Ensure the matrix is SPD
                L{i} = msqrt(Si);
            end
            L_hat = L{1}; 
            iter=0; 

            while true
                sumProcrustes = zeros(d);
                for i = 1:n
                    Li = L{i};
                    [U, ~, V] = svd(L_hat' * Li, 'econ');
                    Ri= U * V';
                    sumProcrustes = sumProcrustes + (Li*Ri);
                end

                L_hat_next = sumProcrustes / n;
                iter = iter + 1;
        
                if norm(L_hat_next - L_hat, 'fro') < tol || iter >= maxIter
                    L_hat = L_hat_next;
                    break;
                end
                L_hat=L_hat_next;
                
            end
            M = (L_hat * L_hat');
            M = (M+M')/2;

        otherwise
            error('Unsupported method. Choose "log-euclidean" or "procrustes".');
    end
end

function L = msqrt(A)
% L = sqrt(A) returns symmetric square root L such that A = L*L'
% uses eigen-decomposition and clips tiny negative eigenvalues to zero.
    [V,D,~] = svd(A);
    d = diag(D);
    L = V * diag(sqrt(d)) * V';
    L = (L + L')/2; %ensure symmetry
end

