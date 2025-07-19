function [T2_map, S0_map] = estimateT2_multipoint_weighted(images, TEs)
% estimateT2_multipoint_weighted estimates T2 relaxation times using weighted linear least squares.
%
% Syntax:
%   [T2_map, S0_map] = estimateT2_multipoint_weighted(images, TEs)
%
% Inputs:
%   images - 4D matrix of size [rows, cols, slices, echoes], containing multi-echo MRI data.
%   TEs    - Vector of echo times (in ms), corresponding to the 4th dimension of `images`.
%
% Outputs:
%   T2_map - Estimated voxel-wise T2 relaxation time map (in ms).
%   S0_map - Estimated voxel-wise initial signal intensity map (S0).
%
% Description:
%   This function fits the mono-exponential decay model:
%       S(TE) = S0 * exp(-TE / T2)
%   to the logarithm of the signal at each voxel using weighted linear least squares,
%   where the weights are proportional to the square of the signal intensity.
%
%   The fitting is performed only for voxels with strictly positive signal values across all echo times.
%   Voxels resulting in non-physical values (e.g., negative, Inf, or NaN) are set to zero.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize
T2_map = zeros(rows, cols, slices);
S0_map = zeros(rows, cols, slices);

% Fit for all the voxels
for i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = squeeze(images(i, j, k, :));

            if all(signal > 0)
                log_signal = log(signal);

                X = [ones(num_echoes, 1), -TEs(:)];

                % Compute weights (higher weight for stronger signal)
                W = diag(signal.^2);

                % Compute Weighted Least Squares solution
                b = (X' * W * X) \ (X' * W * log_signal);

                % Extract S0 and T2 from regression coefficients
                S0_map(i, j, k) = exp(b(1));
                T2_map(i, j, k) = 1 / b(2);

            end
        end
    end
end

T2_map(isinf(T2_map) | isnan(T2_map) | T2_map < 0) = 0;
S0_map(isinf(S0_map) | isnan(S0_map) | S0_map < 0) = 0;

return;
end
