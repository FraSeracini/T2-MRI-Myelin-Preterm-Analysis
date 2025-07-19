function [S0_map, V_maps, mean_residual, RSS] = estimateT2_NNLS_10comp(images, TEs, mask)
% estimateT2_NNLS_10comp estimates multi-compartment T2 volume fractions using NNLS with 10 fixed T2 values.
%
% Syntax:
%   [S0_map, V_maps, mean_residual, RSS] = estimateT2_NNLS_10comp(images, TEs, mask)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] containing multi-echo MRI data.
%   TEs    - Vector of echo times (in milliseconds), one for each echo image.
%   mask   - 3D binary or probabilistic mask indicating which voxels to process.
%
% Outputs:
%   S0_map        - 3D matrix of the total signal amplitude per voxel (sum of components).
%   V_maps        - 4D matrix [rows, cols, slices, 10] of normalized volume fractions for each T2 compartment.
%   mean_residual - Mean absolute residual across all voxels.
%   RSS           - 3D matrix of residual sum of squares per voxel.
%
% Description:
%   This function estimates the T2 distribution using a 10-compartment model, with fixed T2 values:
%       [10, 20, 30, 50, 80, 120, 300, 500, 1000, 2000] ms.
%
%   Each voxel is modeled as:
%       S(TE) = S0 * sum_i (v_i * exp(-TE / T2_i))
%
%   where v_i are the non-negative volume fractions of each compartment, normalized by the total signal S0.
%
%   The NNLS algorithm (`lsqnonneg`) ensures non-negative estimates of each component. Residuals between
%   predicted and actual signals are used to compute the residual sum of squares (RSS) and mean error.
%
% Notes:
%   - Only voxels with strictly positive signal across all echo times and included in the mask are processed.
%   - Residuals are computed voxel-wise across all echo times.
%   - NaN and Inf values are cleaned from outputs.
%   - Parallelization (`parfor`) is used to speed up voxel-wise fitting.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize
residuals = zeros(size(images));
RSS = zeros(96, 96, 55);  % Adjust size if image dimensions vary

% Fixed T2 values (ms) - 10 compartments covering a wide range of T2 values
T2_values = [10, 20, 30, 50, 80, 120, 300, 500, 1000, 2000];

% Design matrix A [TE x compartments]
A_base = exp(-TEs(:) ./ T2_values);  % size: [num_echoes x 10]

% Output maps initialization
S0_map = zeros(rows, cols, slices, 'double');
V_maps = zeros(rows, cols, slices, length(T2_values), 'double');

% Voxel-wise fitting
parfor i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0) && mask(i, j, k) > 0
                % NNLS fitting: solve A * v = y
                v = lsqnonneg(A_base, signal);

                % Estimate S0 as the sum of all components
                S0 = sum(v);
                if S0 > 0
                    v_norm = v / S0;
                else
                    v_norm = zeros(length(T2_values), 1);
                end

                % Save the results
                S0_map(i,j,k) = S0;
                for c = 1:10
                    V_maps(i,j,k,c) = v_norm(c);
                end
            end
        end
    end
end

% Clean NaN and Inf values
S0_map(isnan(S0_map) | isinf(S0_map)) = 0;
V_maps(isnan(V_maps) | isinf(V_maps)) = 0;

% Compute residuals (squared difference between predicted and actual signal)
for t = 1:num_echoes
    A_t = exp(-TEs(t) ./ T2_values);  % 1x10 vector
    predicted = S0_map .* sum(V_maps .* reshape(A_t, [1 1 1 length(T2_values)]), 4);
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
