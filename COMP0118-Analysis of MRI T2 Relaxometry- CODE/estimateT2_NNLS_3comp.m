function [S0_map, V1_map, V2_map, V3_map, mean_residual, RSS] = estimateT2_NNLS_3comp(images, TEs, mask)
% estimateT2_NNLS_3comp estimates compartment volume fractions using a 3-component NNLS model with fixed T2 values.
%
% Syntax:
%   [S0_map, V1_map, V2_map, V3_map, mean_residual, RSS] = ...
%       estimateT2_NNLS_3comp(images, TEs, mask)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] containing multi-echo MRI signal data.
%   TEs    - Vector of echo times (in milliseconds), one per echo image.
%   mask   - 3D binary or probabilistic mask specifying the voxels to include in the fitting.
%
% Outputs:
%   S0_map        - Estimated total signal amplitude per voxel (sum of components).
%   V1_map        - Normalized volume fraction of the short T2 component (myelin-like, T2 = 20 ms).
%   V2_map        - Normalized volume fraction of the medium T2 component (intra/extra-cellular, T2 = 80 ms).
%   V3_map        - Normalized volume fraction of the long T2 component (CSF-like, T2 = 2000 ms).
%   mean_residual - Mean absolute residual error across all voxels.
%   RSS           - 3D matrix of residual sum of squares for each voxel.
%
% Description:
%   This function fits a three-compartment T2 relaxation model to multi-echo MRI data using
%   non-negative least squares (NNLS):
%
%       S(TE) = S0 * [v1 * exp(-TE/T2_1) + v2 * exp(-TE/T2_2) + v3 * exp(-TE/T2_3)]
%
%   where T2_1, T2_2, and T2_3 are fixed to 20 ms, 80 ms, and 2000 ms, respectively, representing
%   myelin water, intra/extra-cellular water, and CSF compartments. The output volume fractions
%   are normalized by the total signal (S0).
%
%   Only voxels with positive signal across all echo times and within the provided mask are considered.
%   NaN or Inf values are cleaned from the outputs, and residuals are computed between predicted and actual signals.
%
% Notes:
%   - Uses `lsqnonneg` to ensure non-negativity of estimated component contributions.
%   - Residuals are calculated voxel-wise and used to compute the overall mean residual error.
%   - Parallel computation (`parfor`) is used to accelerate voxel-wise fitting.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize residuals and RSS maps
residuals = zeros(size(images));
RSS = zeros(96, 96, 55);

% Fixed T2 values (ms) for the three compartments
T2_values = [20, 80, 2000];  % myelin water, intra/extra-cellular water, CSF

% Design matrix A [TE x compartments]
A_base = exp(-TEs(:) ./ T2_values);  % size: [num_echoes x 3]

% Initialize output maps
S0_map = zeros(rows, cols, slices, 'double');
V1_map = zeros(rows, cols, slices, 'double');
V2_map = zeros(rows, cols, slices, 'double');
V3_map = zeros(rows, cols, slices, 'double');

% Fit the model voxel-wise
parfor i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0) && mask(i, j, k) > 0
                % Perform NNLS fit: solve A * v = y
                v = lsqnonneg(A_base, signal);

                % Estimate S0 as the sum of all components
                S0 = sum(v);
                if S0 > 0
                    v_norm = v / S0;
                else
                    v_norm = [0; 0; 0];
                end

                % Store the results
                S0_map(i,j,k) = S0;
                V1_map(i,j,k) = v_norm(1);
                V2_map(i,j,k) = v_norm(2);
                V3_map(i,j,k) = v_norm(3);
            end
        end
    end
end

% Clean NaN and Inf values from the maps
S0_map(isnan(S0_map) | isinf(S0_map)) = 0;
V1_map(isnan(V1_map) | isinf(V1_map)) = 0;
V2_map(isnan(V2_map) | isinf(V2_map)) = 0;
V3_map(isnan(V3_map) | isinf(V3_map)) = 0;

% Compute residuals (squared difference between actual and predicted signals)
for t = 1:num_echoes
    A_t = exp(-TEs(t) ./ T2_values);  % 1x3 vector
    predicted = S0_map .* ( ...
        V1_map * A_t(1) + ...
        V2_map * A_t(2) + ...
        V3_map * A_t(3));
    actual = images(:,:,:,t);
    residuals(:,:,:,t) = (actual - predicted).^2;
end

% Invalidate non-meaningful voxels
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

% Compute the mean residual across all voxels
mean_residual = mean(abs(RSS(:)));

end