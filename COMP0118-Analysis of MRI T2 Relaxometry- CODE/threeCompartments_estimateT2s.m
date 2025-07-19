function [S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map, mean_residual, RSS] = threeCompartments_estimateT2s(images, TEs, mask)
% threeCompartments_estimateT2s estimates voxel-wise T2 parameters using a three-compartment exponential model.
%
% Syntax:
%   [S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map, mean_residual, RSS] = ...
%       threeCompartments_estimateT2s(images, TEs, mask)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] containing multi-echo MRI signal data.
%   TEs    - Vector of echo times (in milliseconds) corresponding to each echo.
%   mask   - 3D binary or probabilistic brain mask specifying voxels for fitting.
%
% Outputs:
%   S0_map        - Estimated initial signal amplitude per voxel.
%   V1_map        - Volume fraction of the short T2 compartment (e.g., myelin water).
%   V2_map        - Volume fraction of the intermediate T2 compartment (e.g., intra/extracellular).
%   V3_map        - Volume fraction of the long T2 compartment (e.g., CSF-like).
%   T2_1_map      - Estimated T2 map (in ms) of the short relaxation component.
%   T2_2_map      - Estimated T2 map (in ms) of the intermediate component.
%   T2_3_map      - Estimated T2 map (in ms) of the long component.
%   mean_residual - Mean absolute residual error across all voxels.
%   RSS           - 3D matrix of residual sum of squares (RSS) per voxel.
%
% Description:
%   This function fits a three-compartment mono-exponential decay model to multi-echo MRI data
%   using non-linear least squares (`lsqcurvefit`) with constraints:
%
%       S(TE) = S0 * (v1 * exp(-TE / T2_1) + v2 * exp(-TE / T2_2) + v3 * exp(-TE / T2_3))
%       with v1 + v2 + v3 = 1
%
%   The T2 relaxation times and corresponding volume fractions (v1, v2, v3) are estimated voxel-wise
%   for all valid voxels (mask > 0 and positive signal values).
%
% Notes:
%   - All volume fractions are constrained to be non-negative and sum to 1.
%   - T2 values are constrained to biologically meaningful ranges (e.g., T2_1 ∈ [5, 60] ms).
%   - Invalid values (NaN, Inf, or S0 ≤ 0) are zeroed out in the output.
%   - Residuals are computed per echo, per voxel and used to compute RSS and mean error.
%   - Parallel computation (`parfor`) is used to accelerate processing.
%

warning('off', 'MATLAB:nearlySingularMatrix');
[rows, cols, slices, num_echoes] = size(images);

% Initialize
residuals = zeros(size(images));
RSS = zeros(96, 96, 55);


% Output maps
S0_map   = zeros(rows, cols, slices, 'double');
V1_map   = zeros(rows, cols, slices, 'double');
V2_map   = zeros(rows, cols, slices, 'double');
V3_map   = zeros(rows, cols, slices, 'double');
T2_1_map = zeros(rows, cols, slices, 'double');
T2_2_map = zeros(rows, cols, slices, 'double');
T2_3_map = zeros(rows, cols, slices, 'double');

TEs = double(TEs);

% Fit the model for all the voxels
parfor i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0) && mask(i, j, k) > 0

                % x = [S0, v1, v2, v3, T2_1, T2_2, T2_3]
                x0 = [max(signal), 0.3, 0.5, 0.2, 20, 80, 2000];
                lb = [0,     0,   0,   0,   5,   65,  1800];
                ub = [Inf,   1,   1,   1, 60, 100, 2500];

                Aeq = [0 1 1 1 0 0 0];
                beq = 1;

                model_fun = @(x, TE) x(1) * ( ...
                    x(2)*exp(-TE/x(5)) + ...
                    x(3)*exp(-TE/x(6)) + ...
                    x(4)*exp(-TE/x(7)) );

                options = optimoptions('lsqcurvefit', 'Display', 'off');

                try
                    params = lsqcurvefit(model_fun, x0, TEs(:), signal(:), lb, ub, [], [], Aeq, beq, [], options);

                    S0_map(i, j, k)   = params(1);
                    V1_map(i, j, k)   = params(2);
                    V2_map(i, j, k)   = params(3);
                    V3_map(i, j, k)   = params(4);
                    T2_1_map(i, j, k) = params(5);
                    T2_2_map(i, j, k) = params(6);
                    T2_3_map(i, j, k) = params(7);

                catch
                    continue
                end
            end
        end
    end
end

% Clean invalid
S0_map(isinf(S0_map) | isnan(S0_map)) = 0;
V1_map(isinf(V1_map) | isnan(V1_map)) = 0;
V2_map(isinf(V2_map) | isnan(V2_map)) = 0;
V3_map(isinf(V3_map) | isnan(V3_map)) = 0;

T2_1_map(isinf(T2_1_map) | isnan(T2_1_map)) = 0;
T2_2_map(isinf(T2_2_map) | isnan(T2_2_map)) = 0;
T2_3_map(isinf(T2_3_map) | isnan(T2_3_map)) = 0;

% Residuals for all voxels and TEs
for t = 1:num_echoes
    predicted = S0_map .* ( ...
        V1_map .* exp(-TEs(t) ./ T2_1_map) + ...
        V2_map .* exp(-TEs(t) ./ T2_2_map) + ...
        V3_map .* exp(-TEs(t) ./ T2_3_map));
    actual = images(:,:,:,t);
    residuals(:,:,:,t) = (actual - predicted).^2;
end

invalid_mask = (S0_map <= 0) | isnan(S0_map);
for t = 1:num_echoes
    temp = residuals(:,:,:,t);
    temp(invalid_mask) = 0;
    residuals(:,:,:,t) = temp;
end

% Residuals for all voxels
for i = 1 : rows
    for j = 1 : cols
        for k = 1 : slices
            RSS(i, j, k) = sum(residuals(i, j, k, :));
        end
    end
end

% Average error
mean_residual = mean(abs(RSS(:)));

return;
end
