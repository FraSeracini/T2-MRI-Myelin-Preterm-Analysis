function [S0_map, V1_map, V2_map, V3_map, V4_map, T2_1_map, T2_2_map, T2_3_map, T2_4_map, mean_residual, RSS] = ...
    fourCompartments_estimateT2s(images, TEs, mask)

% fourCompartments_estimateT2s estimates T2 relaxation parameters using a four-compartment model with free T2s.
%
% Syntax:
%   [S0_map, V1_map, V2_map, V3_map, V4_map, T2_1_map, T2_2_map, T2_3_map, T2_4_map, mean_residual, RSS] = ...
%       fourCompartments_estimateT2s(images, TEs, mask)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] containing multi-echo MRI data.
%   TEs    - Vector of echo times (in milliseconds), one for each echo image.
%   mask   - 3D binary or probabilistic mask defining the voxels to include in the fitting.
%
% Outputs:
%   S0_map        - Estimated signal amplitude per voxel.
%   V1_map        - Volume fraction map of the first T2 component.
%   V2_map        - Volume fraction map of the second T2 component.
%   V3_map        - Volume fraction map of the third T2 component.
%   V4_map        - Volume fraction map of the fourth T2 component.
%   T2_1_map      - Estimated T2 map (in ms) of the first compartment.
%   T2_2_map      - Estimated T2 map (in ms) of the second compartment.
%   T2_3_map      - Estimated T2 map (in ms) of the third compartment.
%   T2_4_map      - Estimated T2 map (in ms) of the fourth compartment.
%   mean_residual - Mean absolute residual error across all voxels.
%   RSS           - 3D matrix of residual sum of squares for each voxel.
%
% Description:
%   This function fits a four-compartment signal decay model to multi-echo MRI data using
%   non-linear least squares (`lsqcurvefit`). The model has the form:
%
%       S(TE) = S0 * (v1 * exp(-TE/T2_1) + v2 * exp(-TE/T2_2) + v3 * exp(-TE/T2_3) + v4 * exp(-TE/T2_4))
%
%   where the T2s and volume fractions (v1...v4) are freely estimated for each voxel, under the constraint:
%
%       v1 + v2 + v3 + v4 = 1
%
%   This allows the model to represent a range of tissue compartments (e.g., myelin water, intra/extracellular,
%   slow-exchanging components, and CSF). Initial guesses and bounds are provided to guide the optimization.
%
% Notes:
%   - Only voxels with all-positive signal across all TEs and within the mask are processed.
%   - Voxels with invalid or missing results are zeroed out in the output maps.
%   - Residuals are calculated for each echo and voxel to compute voxelwise and global errors.
%   - Parallel processing (`parfor`) is used to accelerate fitting across voxels.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize
residuals = zeros(size(images));
RSS = zeros(96, 96, 55);

% Output maps
S0_map   = zeros(rows, cols, slices, 'double');
V1_map   = zeros(rows, cols, slices, 'double');
V2_map   = zeros(rows, cols, slices, 'double');
V3_map   = zeros(rows, cols, slices, 'double');
V4_map   = zeros(rows, cols, slices, 'double');

T2_1_map = zeros(rows, cols, slices, 'double');
T2_2_map = zeros(rows, cols, slices, 'double');
T2_3_map = zeros(rows, cols, slices, 'double');
T2_4_map = zeros(rows, cols, slices, 'double');

TEs = double(TEs);

% Fit the model for all the voxels
parfor i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0) && mask(i, j, k) > 0

                % x = [S0, v1, v2, v3, v4, T2_1, T2_2, T2_3, T2_4]
                x0 = [max(signal), 0.2, 0.3, 0.3, 0.2, 15, 50, 300, 2000];
                lb = [0, 0, 0, 0, 0,  5,  10, 100, 1500];
                ub = [Inf, 1, 1, 1, 1, 40, 100, 1000, 3000];

                % Vincolo: v1 + v2 + v3 + v4 = 1
                Aeq = [0 1 1 1 1 0 0 0 0];
                beq = 1;

                model_fun = @(x, TE) x(1) * ( ...
                    x(2)*exp(-TE/x(6)) + ...
                    x(3)*exp(-TE/x(7)) + ...
                    x(4)*exp(-TE/x(8)) + ...
                    x(5)*exp(-TE/x(9)) );

                options = optimoptions('lsqcurvefit', 'Display', 'off');

                try
                    params = lsqcurvefit(model_fun, x0, TEs(:), signal(:), lb, ub, [], [], Aeq, beq, [], options);

                    S0_map(i, j, k)   = params(1);
                    V1_map(i, j, k)   = params(2);
                    V2_map(i, j, k)   = params(3);
                    V3_map(i, j, k)   = params(4);
                    V4_map(i, j, k)   = params(5);
                    T2_1_map(i, j, k) = params(6);
                    T2_2_map(i, j, k) = params(7);
                    T2_3_map(i, j, k) = params(8);
                    T2_4_map(i, j, k) = params(9);
                catch
                    continue
                end
            end
        end
    end
end

% Clean invalid
fields = {'S0_map','V1_map','V2_map','V3_map','V4_map',...
          'T2_1_map','T2_2_map','T2_3_map','T2_4_map'};
for f = fields
    m = eval(f{1});
    m(isnan(m) | isinf(m)) = 0;
    eval([f{1} ' = m;']);
end

% Residuals for all the voxels and TEs
for t = 1:num_echoes
    predicted = S0_map .* ( ...
        V1_map .* exp(-TEs(t) ./ T2_1_map) + ...
        V2_map .* exp(-TEs(t) ./ T2_2_map) + ...
        V3_map .* exp(-TEs(t) ./ T2_3_map) + ...
        V4_map .* exp(-TEs(t) ./ T2_4_map));
    actual = images(:,:,:,t);
    residuals(:,:,:,t) = (actual - predicted).^2;
end

% Clean non-valid values
invalid_mask = (S0_map <= 0) | isnan(S0_map);
for t = 1:num_echoes
    temp = residuals(:,:,:,t);
    temp(invalid_mask) = 0;
    residuals(:,:,:,t) = temp;
end

% Residuals for all the voxels
for i = 1 : 96
    for j = 1 : 96
        for k = 1 : 55
            RSS(i, j, k) = sum(residuals(i, j, k, :));
        end
    end
end

% Average error
mean_residual = mean(abs(RSS(:)));

return;
end
