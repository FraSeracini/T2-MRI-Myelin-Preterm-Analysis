function [T2_1_map, T2_2_map, S0_map, V1_map, residuals, mean_residual, RSS] = estimateT2_multipoint_NLLS_bicomponent(images, TEs, mask)
% estimateT2_multipoint_NLLS_bicomponent estimates bi-exponential T2 parameters using non-linear least squares (NLLS).
%
% Syntax:
%   [T2_1_map, T2_2_map, S0_map, V1_map, residuals, mean_residual, RSS] = ...
%       estimateT2_multipoint_NLLS_bicomponent(images, TEs, mask)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] of multi-echo MRI data.
%   TEs    - Vector of echo times (in milliseconds), one per echo image.
%   mask   - 3D binary or probabilistic mask defining voxels to include in the fitting.
%
% Outputs:
%   T2_1_map      - Estimated T2 map (in ms) for the short T2 component.
%   T2_2_map      - Estimated T2 map (in ms) for the long T2 component.
%   S0_map        - Estimated initial signal amplitude (S0) map.
%   V1_map        - Volume fraction map of the short T2 component (between 0 and 1).
%   residuals     - 4D matrix of squared residuals for each voxel and echo time.
%   mean_residual - Mean residual error across all voxels (scalar).
%   RSS           - Residual Sum of Squares (3D map), computed across all echo times.
%
% Description:
%   This function fits a bi-exponential signal decay model to multi-echo MRI data:
%       S(TE) = S0 * [v1 * exp(-TE/T2_1) + (1 - v1) * exp(-TE/T2_2)]
%   using non-linear least squares (NLLS) optimization (`lsqcurvefit`) with
%   lower and upper bounds on the parameters.
%
%   For each voxel within the specified mask and with all positive signal values,
%   the algorithm estimates four parameters:
%       - S0  : initial signal amplitude
%       - v1  : volume fraction of the short T2 component
%       - T2_1: short T2 relaxation time
%       - T2_2: long T2 relaxation time
%
%   Residuals are computed between the observed signal and the model prediction,
%   and the residual sum of squares (RSS) is calculated for each voxel.
%
% Notes:
%   - Voxels with non-positive or invalid S0 values are excluded.
%   - Output maps are filtered to remove NaN, Inf, or negative values.
%   - The initial parameter guess and bounds can be adjusted for different datasets.
%
[rows, cols, slices, num_echoes] = size(images);

% Initialize
residuals = zeros(size(images));
RSS = zeros(96, 96, 55);


% Ensure TEs is double precision
TEs = double(TEs);

% Initialize output maps as double
T2_1_map = zeros(rows, cols, slices, 'double');
T2_2_map = zeros(rows, cols, slices, 'double');
S0_map   = zeros(rows, cols, slices, 'double');
V1_map   = zeros(rows, cols, slices, 'double');

% Fit for all the voxels
parfor i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0)
                if(mask(i, j, k) > 0)

                    % Starting point: [S0, v1, T2_1, T2_2]
                    p0 = [max(signal), 0.5, 20, 70];
                    lb = [0, 0, 10, 60];
                    ub = [inf, 1, 60, 2500];

                    % Non-linear model: S(TE) = S0 * [v1 * exp(-TE/T2_1) + (1 - v1) * exp(-TE/T2_2)]
                    objective_fun = @(p, TE) p(1) * (p(2) * exp(-TE / p(3)) + (1 - p(2)) * exp(-TE / p(4)));

                    % Fit with non-linear least squares
                    options = optimset('Display', 'off', 'MaxFunEvals', 100);

                    try
                        params = lsqcurvefit(objective_fun, p0, TEs(:), signal(:), lb, ub, options);

                        S0_map(i, j, k)   = params(1);
                        V1_map(i, j, k)   = params(2);
                        T2_1_map(i, j, k) = params(3);
                        T2_2_map(i, j, k) = params(4);

                    catch
                        continue
                    end
                end
            end
        end
    end
end

% Clean up NaNs and Infs
T2_1_map(isinf(T2_1_map) | isnan(T2_1_map) | T2_1_map < 0) = 0;
T2_2_map(isinf(T2_2_map) | isnan(T2_2_map) | T2_2_map < 0) = 0;
S0_map(isinf(S0_map) | isnan(S0_map) | S0_map < 0) = 0;
V1_map(isinf(V1_map) | isnan(V1_map) | V1_map < 0) = 0;

for t = 1:num_echoes
    predicted = S0_map .* (V1_map.*exp(-TEs(t) ./ T2_1_map) +(1 - V1_map).*exp(-TEs(t) ./ T2_2_map)) ;
    actual = images(:,:,:,t);
    residuals(:,:,:,t) = (actual - predicted).^2;
end

% Invalida voxels non significativi
invalid_mask = (S0_map <= 0) | isnan(S0_map);
for t = 1:num_echoes
    temp = residuals(:,:,:,t);
    temp(invalid_mask) = 0;
    residuals(:,:,:,t) = temp;
end

% Calculate the Sum of Square differences for all the voxels
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