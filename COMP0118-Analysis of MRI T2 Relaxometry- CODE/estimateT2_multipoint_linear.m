function [T2_map, S0_map] = estimateT2_multipoint_linear(images, TEs)
% estimateT2_multipoint_linear estimates T2 relaxation times using linearized least squares.
%
% Syntax:
%   [T2_map, S0_map] = estimateT2_multipoint_linear(images, TEs)
%
% Inputs:
%   images - 4D matrix of size [rows, cols, slices, echoes], representing multi-echo MRI data.
%   TEs    - Vector of echo times (in ms), one for each image volume in the 4th dimension.
%
% Outputs:
%   T2_map - Estimated voxel-wise T2 relaxation time map (in ms), same spatial size as input images.
%   S0_map - Estimated voxel-wise initial signal (S0) map (same size as T2_map).
%
% Description:
%   This function estimates the T2 relaxation time and S0 value for each voxel
%   by applying a linear least squares fit to the logarithm of the signal across
%   multiple echo times. The model used is:
%       log(S(TE)) = log(S0) - TE / T2
%   which becomes a linear equation when log-transforming the signal.
%
%   Only voxels with positive signal values at all echo times are included in the fit.
%   The function handles non-physical results by zeroing negative, NaN, or Inf values.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize
T2_map = zeros(rows, cols, slices);
S0_map = zeros(rows, cols, slices);

% fit for every voxel
for i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = squeeze(images(i, j, k, :));

            if all(signal > 0)

                % Logarithmic linearization
                log_signal = log(signal);

                X = [ones(num_echoes, 1), -TEs(:)];
                b = X \ log_signal;

                S0_map(i, j, k) = exp(b(1));
                T2_map(i, j, k) = 1/b(2);

            end
        end
    end
end

% set the value of negative/ NaN/ Inf as 0
T2_map(isinf(T2_map) | isnan(T2_map) | T2_map < 0) = 0;
S0_map(isinf(S0_map) | isnan(S0_map) | S0_map < 0) = 0;

return;
end
