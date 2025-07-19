function [T2_1_map, T2_2_map, S0_map, V1_map, residuals, mean_residual, RSS] = estimateT2_fixedT2s(images, TEs, mask)
% estimateT2_fixedT2s estimates S0 and volume fraction (v1) using a bi-exponential model with fixed T2 values.
%
% Syntax:
%   [T2_1_map, T2_2_map, S0_map, V1_map, residuals, mean_residual, RSS] = ...
%       estimateT2_fixedT2s(images, TEs, mask)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] containing multi-echo MRI data.
%   TEs    - Vector of echo times (in milliseconds), matching the 4th dimension of `images`.
%   mask   - 3D binary or probabilistic mask defining voxels to include in the fitting.
%
% Outputs:
%   T2_1_map      - 3D matrix with constant value equal to the fixed short T2 component (60 ms).
%   T2_2_map      - 3D matrix with constant value equal to the fixed long T2 component (2000 ms).
%   S0_map        - Estimated initial signal intensity per voxel.
%   V1_map        - Estimated volume fraction of the short T2 component per voxel.
%   residuals     - 4D matrix of squared residuals between observed and predicted signal.
%   mean_residual - Mean absolute residual across all voxels.
%   RSS           - Total residual sum of squares over the entire volume.
%
% Description:
%   This function fits a simplified bi-exponential model of the form:
%       S(TE) = S0 * [v1 * exp(-TE / T2_1) + (1 - v1) * exp(-TE / T2_2)]
%   where T2_1 and T2_2 are fixed (default: 60 ms and 2000 ms, respectively).
%   The model estimates only the parameters S0 and v1 using non-linear least squares.
%
%   The fitting is done voxel-wise for all valid voxels (positive signal, within the mask).
%   Invalid or non-physical values (e.g., NaN, Inf, negative) are set to zero.
%
% Notes:
%   - Uses `lsqcurvefit` with bounds on S0 and v1.
%   - Parallel processing (parfor) is used for speed.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize
residuals = zeros(size(images));
RSS = zeros(96, 96, 55);

% Fixed T2 values (short and long components)
T2_1_fixed = 60;
T2_2_fixed = 2000;

% Ensure TEs is double precision
TEs = double(TEs);

% Initialize output maps
T2_1_map = T2_1_fixed * ones(rows, cols, slices, 'double');
T2_2_map = T2_2_fixed * ones(rows, cols, slices, 'double');
S0_map   = zeros(rows, cols, slices, 'double');
V1_map   = zeros(rows, cols, slices, 'double');

% Fit for all the voxels
parfor i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0)
                if(mask(i, j, k) > 0)

                    % Parameters to estimate: [S0, v1]
                    p0 = [max(signal), 0.5];
                    lb = [0, 0];
                    ub = [inf, 1];

                    % Fixed T2 bicomponent model
                    model_fun = @(p, TE) p(1) * (p(2) * exp(-TE / T2_1_fixed) + (1 - p(2)) * exp(-TE / T2_2_fixed));

                    options = optimset('Display', 'off', 'MaxFunEvals', 100);
                    try
                        params = lsqcurvefit(model_fun, p0, TEs(:), signal(:), lb, ub, options);

                        S0_map(i, j, k) = params(1);
                        V1_map(i, j, k) = params(2);
                    
                    catch
                        continue
                    end
                end
            end
        end
    end
end

% Clean up NaNs and Infs
S0_map(isinf(S0_map) | isnan(S0_map) | S0_map < 0) = 0;
V1_map(isinf(V1_map) | isnan(V1_map) | V1_map < 0) = 0;

% Compute predicted signal and residuals
for t = 1:num_echoes
    predicted = S0_map .* (V1_map .* exp(-TEs(t) ./ T2_1_fixed) + (1 - V1_map) .* exp(-TEs(t) ./ T2_2_fixed));
    actual = images(:,:,:,t);
    residuals(:,:,:,t) = (actual - predicted).^2;
end

% Eliminate out of scale values
invalid_mask = (S0_map <= 0) | isnan(S0_map);
for t = 1:num_echoes
    temp = residuals(:,:,:,t);
    temp(invalid_mask) = 0;
    residuals(:,:,:,t) = temp;
end

% Compute the sum of Square differences
for i = 1 : 96
    for j = 1 : 96
        for k = 1 : 55
            RSS(i, j, k) = sum(residuals(i, j, k, :));
        end
    end
end

% Compute the average error
mean_residual = mean(abs(RSS(:)));

return;
end
