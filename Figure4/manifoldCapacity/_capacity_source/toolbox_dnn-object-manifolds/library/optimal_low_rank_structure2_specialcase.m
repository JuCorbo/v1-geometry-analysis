function [Vopt, Xopt, Kopt, residual_centers_norm, mean_square_corrcoef, mean_abs_corrcoef, mean_square_corr, mean_abs_corr, results] = ...
    optimal_low_rank_structure2_specialcase(X, MAX_K, verbose, minSquare, N_REPEATS)
    % Problem dimensions
    [N, P] = size(X);

    % Handle special case when P = 2
    if P == 2
        warning('P = 2: Handling limited manifold geometry.');

        % Initialize results
        residual_centers_norm = nan(MAX_K + 1, P);  % Match expected dimensions
        residual_centers_norm(1, :) = sqrt(sum(X.^2, 1));  % Norm of original columns

        % Compute covariance matrix
        C = X' * X;  % Covariance of input data

        % Compute centroid distance
        centroid_distance = norm(X(:, 1) - X(:, 2));  % Distance between the two columns (centroids)

        % Compute overlap
        norm1 = sqrt(C(1, 1));
        norm2 = sqrt(C(2, 2));
        if norm1 == 0 || norm2 == 0
            warning('Centroid norm is zero. Overlap cannot be computed.');
            overlap = NaN;
        else
            overlap = abs(C(1, 2)) / (norm1 * norm2);  % Normalized overlap
        end

        % incorrect computation of radii!
        % Compute intra-class variability (manifold radii) for each class
        % Rm_class1 = mean(sqrt(sum((X(:, 1) - mean(X(:, 1), 1)).^2, 1)));  % Class 1 (first column)
        % Rm_class2 = mean(sqrt(sum((X(:, 2) - mean(X(:, 2), 1)).^2, 1)));  % Class 2 (second column)
        % % Average manifold radii
        % Rm = (Rm_class1 + Rm_class2) / 2;

        % Debugging output
        % fprintf('Centroid Distance: %1.3f, Overlap: %1.3f, Rm: %1.3f\n', centroid_distance, overlap, Rm);
        % 
        % Create a struct to save the metrics
        results = struct();
        results.centroid_distance = centroid_distance;
        results.overlap = overlap;
        % results.Rm_class1 = Rm_class1;
        % results.Rm_class2 = Rm_class2;
        % results.Rm = Rm;
        
        % Save the results
        %save(sprintf('%s_pairwise_metrics.mat', prefix), 'results');
        
        % Set other outputs as NaN for consistency
        Vopt = [];
        Xopt = X;
        Kopt = 0;
        mean_square_corr = nan(MAX_K + 1, 1);
        mean_abs_corr = nan(MAX_K + 1, 1);
        mean_square_corrcoef = nan(MAX_K + 1, 1);
        mean_abs_corrcoef = nan(MAX_K + 1, 1);

        return;
    end

    % Default parameter values
    if nargin < 2
        MAX_K = ceil(P / 2);
    end
    if nargin < 3
        verbose = 1;
    end
    if nargin < 4
        minSquare = true;  % T: minimize square correlations, F: minimize abs correlations
    end
    if nargin < 5
        N_REPEATS = 1;  % Number of optimization repeats, taking the most stable one
    end

    early_termination = true;
    % Results variables
    mean_square_corr = nan(MAX_K + 1, 1);
    mean_abs_corr = nan(MAX_K + 1, 1);
    mean_square_corrcoef = nan(MAX_K + 1, 1);
    mean_abs_corrcoef = nan(MAX_K + 1, 1);
    residual_centers_norm = nan(MAX_K + 1, P);

    % Add optimization package
    path_to_add = [pwd, '/../FOptM/'];
    if ~exist(path_to_add, 'dir')
        path_to_add = [pwd, '/FOptM/'];
    end
    path(path, path_to_add);

    % Optimization options
    MAX_ITER = 10000;
    opts = struct;
    opts.record = 0;  % no printout
    opts.mxitr = MAX_ITER;  % max number of iterations
    opts.gtol = 1e-6;  % stop control for the projected gradient
    opts.xtol = 1e-6;  % stop control for ||X_k - X_{k-1}||
    opts.ftol = 1e-8;  % stop control for |F_k - F_{k-1}| / (1 + |F_{k-1}|)

    % Move from representation from dimension of N to P-1, the rank
    if N > P - 1
        [Q, ~] = qr(X, 0);
        Q = Q(:, 1:P - 1);  % [N, P-1]
        Xq = Q' * X;
    else
        Xq = X;
        Q = eye(N);
    end

    % Initialize cost to infinite
    bestCost = inf;
    Vopt = [];
    Kopt = 0;

    % Main optimization loop
    for ik = 1:MAX_K + 1
        k = ik - 1;

        % Early termination
        if early_termination && (k > Kopt + 3)
            fprintf('Early termination Kopt=%d\n', Kopt);
            break;
        end

        tic;

        if k == 0
            V = [];
            Xk = Xq;
        else
            % Optimization for k > 0
            best_stability = 0;
            best_V = [];
            for in = 1:N_REPEATS
                s = randn(P, 1);
                V0 = [Xq * s, V];
                [V0, ~] = qr(V0, 0);
                [V1, output] = OptStiefelGBB(V0, @square_corrcoeff_full_cost, opts, Xq');
                assert(output.itr < MAX_ITER, 'Max iterations reached at k=%d', k);
                Xk = Xq - V1 * (V1' * Xq);
                stability = min(sqrt(sum(Xk.^2, 1)) ./ sqrt(sum(Xq.^2, 1)));
                if stability > best_stability
                    best_stability = stability;
                    best_V = V1;
                end
            end
            V = best_V;
            Xk = Xq - V * (V' * Xq);
        end

        % Save output variables
        Xk_norm = sqrt(sum(Xk.^2, 1));
        residual_centers_norm(ik, :) = Xk_norm;

        Ck = Xk' * Xk;
        square_offdiagonal_corr = (Ck - diag(diag(Ck))).^2;
        mean_square_corr(ik) = sum(square_offdiagonal_corr(:)) / (P - 1) / P;

        abs_offdiagonal_corr = abs(Ck - diag(diag(Ck)));
        mean_abs_corr(ik) = sum(abs_offdiagonal_corr(:)) / (P - 1) / P;

        % Update best results
        currentCost = minSquare * mean_square_corr(ik) + ~minSquare * mean_abs_corr(ik);
        if currentCost < bestCost
            bestCost = currentCost;
            Vopt = Q * V;
            Xopt = Q * Xk;
            Kopt = k;
        end
    end
end
