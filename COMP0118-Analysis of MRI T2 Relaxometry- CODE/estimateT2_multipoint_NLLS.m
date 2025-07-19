function [T2_map, S0_map] = estimateT2_multipoint_NLLS(images, TEs)
% estimateT2_multipoint_NLLS estimates T2 relaxation times using non-linear least squares (NLLS).
%
% Syntax:
%   [T2_map, S0_map] = estimateT2_multipoint_NLLS(images, TEs)
%
% Inputs:
%   images - 4D matrix of size [rows, cols, slices, echoes], containing multi-echo MRI signal intensities.
%   TEs    - Vector of echo times (in milliseconds), corresponding to the 4th dimension of `images`.
%
% Outputs:
%   T2_map - 3D matrix of estimated T2 relaxation times (in ms) for each voxel.
%   S0_map - 3D matrix of estimated initial signal intensities (S0) for each voxel.
%
% Description:
%   This function estimates voxel-wise T2 relaxation times by fitting the
%   mono-exponential signal decay model using non-linear least squares:
%       S(TE) = S0 * exp(-TE / T2)
%
%   For each voxel with strictly positive signal values across all echo times,
%   a two-parameter model [S0, T2] is fit using `fminsearch` to minimize the
%   squared error between observed and predicted signals.
%
%   Invalid results (e.g., NaN, Inf, or negative estimates) are set to zero.
%   The initial guess for the optimization is based on the maximum signal
%   intensity and a default T2 value of 50 ms.
%


[rows, cols, slices, num_echoes] = size(images);

% Ensure TEs is double precision
TEs = double(TEs);

% Initialize output maps as double
T2_map = zeros(rows, cols, slices, 'double');
S0_map = zeros(rows, cols, slices, 'double');

% Initial parameter guess: [S0, T2]
initial_guess = double([max(images(:)), 50]);  % Ensure double precision

for i = 1:rows
    for j = 1:cols
        for k = 1:slices
            signal = squeeze(images(i, j, k, :));

            if all(signal > 0)
                signal = double(signal);

                % Define the objective function
                objective_fun = @(b) sum((signal(:) - (b(1) * exp(-TEs(:) / b(2)))).^2);

                % Solve non-linear least squares using fminsearch
                options = optimset('Display', 'off');
                params = fminsearch(objective_fun, initial_guess, options);

                S0_map(i, j, k) = max(params(1), 0);
                T2_map(i, j, k) = max(params(2), 0);

            end
        end
    end
end

T2_map(isinf(T2_map) | isnan(T2_map) | T2_map < 0) = 0;
S0_map(isinf(S0_map) | isnan(S0_map) | S0_map < 0) = 0;

return;
end
