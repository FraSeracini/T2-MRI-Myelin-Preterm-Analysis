function [S0_map, V1_map, V2_map, residuals, mean_residual, RSS] = estimateT2_NNLS_2comp(images, TEs, mask)
% estimateT2_NNLS_2comp estimates T2 compartment volume fractions using non-negative least squares (NNLS).
%
% Syntax:
%   [S0_map, V1_map, V2_map, residuals, mean_residual, RSS] = ...
%       estimateT2_NNLS_2comp(images, TEs, mask)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] of multi-echo MRI data.
%   TEs    - Vector of echo times (in milliseconds), one for each echo.
%   mask   - 3D binary or probabilistic brain mask to restrict voxel-wise fitting.
%
% Outputs:
%   S0_map        - Estimated total signal amplitude per voxel (sum of components).
%   V1_map        - Estimated normalized volume fraction of the first component (T2 = 60 ms).
%   V2_map        - Estimated normalized volume fraction of the second component (T2 = 2000 ms).
%   residuals     - 4D matrix of squared residuals between actual and predicted signals.
%   mean_residual - Mean absolute residual across all voxels.
%   RSS           - 3D matrix of residual sum of squares per voxel.
%
% Description:
%   This function fits a two-component exponential decay model to multi-echo MRI
%   signals using non-negative least squares (NNLS):
%
%       S(TE) = S0 * [v1 * exp(-TE/T2_1) + v2 * exp(-TE/T2_2)]
%
%   where T2_1 and T2_2 are fixed to 60 ms and 2000 ms, respectively (representing
%   intra-/extra-cellular and CSF compartments).
%
%   NNLS is used to ensure the solution is physically meaningful (non-negative volumes).
%   The estimated fractions are normalized by the sum (S0), and residuals are computed
%   for each voxel across all echo times.
%
% Notes:
%   - Voxels with zero or invalid signal are excluded.
%   - Residuals are set to zero for invalid voxels.
%   - Parallel processing (`parfor`) is used for efficient voxel-wise computation.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize residuals and RSS maps
residuals = zeros(size(images));
RSS = zeros(96, 96, 55);

% Fixed T2 values (ms) - two components: intra/extra-cellular and CSF
T2_values = [60, 2000];  % intra/extra-cellular water, cerebrospinal fluid (CSF)

% Matrix A [TE x compartments]
A_base = exp(-TEs(:) ./ T2_values);  % size: [num_echoes x 2]

% Output maps initialization
S0_map = zeros(rows, cols, slices, 'double');
V1_map = zeros(rows, cols, slices, 'double');
V2_map = zeros(rows, cols, slices, 'double');

% Parallel loop over all voxels
parfor i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0) && mask(i, j, k) > 0
                % NNLS fitting: solves A * v = y
                v = lsqnonneg(A_base, signal);

                % Estimate S0 as the sum of all components
                S0 = sum(v);
                if S0 > 0
                    v_norm = v / S0;
                else
                    v_norm = [0; 0];
                end

                % Store the results
                S0_map(i,j,k) = S0;
                V1_map(i,j,k) = v_norm(1);
                V2_map(i,j,k) = v_norm(2);
            end
        end
    end
end

% Remove NaN and Inf values from the outputs
S0_map(isnan(S0_map) | isinf(S0_map)) = 0;
V1_map(isnan(V1_map) | isinf(V1_map)) = 0;
V2_map(isnan(V2_map) | isinf(V2_map)) = 0;

% Compute residuals (squared differences)
for t = 1:num_echoes
    A_t = exp(-TEs(t) ./ T2_values);  % 1x2 vector
    predicted = S0_map .* ( ...
        V1_map * A_t(1) + ...
        V2_map * A_t(2));
    actual = images(:,:,:,t);
    residuals(:,:,:,t) = (actual - predicted).^2;
end

% Invalidate non-meaningful voxels (e.g., background)
invalid_mask = (S0_map <= 0) | isnan(S0_map);
for t = 1:num_echoes
    temp = residuals(:,:,:,t);
    temp(invalid_mask) = 0;
    residuals(:,:,:,t) = temp;
end

% Compute the Residual Sum of Squares (RSS) per voxel
for i = 1 : 96
    for j = 1 : 96
        for k = 1 : 55
            RSS(i, j, k) = sum(residuals(i, j, k, :));
        end
    end
end

% Compute the average residual error across all voxels
mean_residual = mean(abs(RSS(:)));
end