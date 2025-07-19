function [T2_map, S0_map] = estimateT2_multipoint_NNLS(images, TEs)
% estimateT2_multipoint_NNLS estimates T2 relaxation times using non-negative least squares (NNLS).
%
% Syntax:
%   [T2_map, S0_map] = estimateT2_multipoint_NNLS(images, TEs)
%
% Inputs:
%   images - 4D matrix of size [rows, cols, slices, echoes], containing multi-echo MRI signal intensities.
%   TEs    - Vector of echo times (in milliseconds), one for each echo in the 4th dimension of `images`.
%
% Outputs:
%   T2_map - 3D matrix of estimated T2 relaxation times (in ms) per voxel.
%   S0_map - 3D matrix of estimated initial signal intensities (S0) per voxel.
%
% Description:
%   This function estimates voxel-wise T2 relaxation times based on a linearized
%   mono-exponential decay model:
%       S(TE) = S0 * exp(-TE / T2)
%   by solving the log-transformed system using non-negative least squares (NNLS).
%
%   The use of NNLS ensures that the fitted parameters remain non-negative,
%   which is physically meaningful in the context of MRI signals.
%
%   Only voxels with strictly positive signal values across all echo times are considered.
%   Invalid or non-physical results (NaN, Inf, negative) are replaced with zeros.
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize
T2_map = zeros(rows, cols, slices, 'double'); % Ensure double precision
S0_map = zeros(rows, cols, slices, 'double'); % Ensure double precision

% Fit for every voxel
for i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = squeeze(images(i, j, k, :));

            if all(signal > 0)

                log_signal = log(double(signal));  % Convert to double before log

                X = double([ones(num_echoes, 1), -TEs(:)]);

                b = lsqnonneg(X, log_signal);

                S0_map(i, j, k) = exp(b(1));  % S0 = exp(intercept)
                T2_map(i, j, k) = 1 / max(b(2), eps);  % T2 = 1 / slope

            end
        end
    end
end

% Handle invalid values
T2_map(isinf(T2_map) | isnan(T2_map) | T2_map < 0) = 0;
S0_map(isinf(S0_map) | isnan(S0_map) | S0_map < 0) = 0;

return;
end
