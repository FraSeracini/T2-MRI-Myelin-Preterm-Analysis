function residuals = calculate_residuals(images, TEs, T2_map, S0_map)
% calculate_residuals computes voxel-wise squared residuals between actual and predicted signals for the mono-compartment model.
%
% Syntax:
%   residuals = calculate_residuals(images, TEs, T2_map, S0_map)
%
% Inputs:
%   images  - 4D matrix [rows, cols, slices, echoes] of original MRI signal intensities at different echo times.
%   TEs     - Vector of echo times (ms), one per echo image (must match the 4th dimension of `images`).
%   T2_map  - 3D matrix of estimated T2 relaxation times (ms) for each voxel.
%   S0_map  - 3D matrix of estimated initial signal intensities (S0) for each voxel.
%
% Outputs:
%   residuals - 4D matrix of squared residuals with the same size as `images`,
%               computed as (actual - predicted)^2 for each echo time.
%
% Description:
%   This function simulates the MRI signal at each echo time using the mono-exponential model:
%       S(TE) = S0 * exp(-TE / T2)
%   It then computes the squared residual between the predicted and actual signal for every voxel and echo.
%   Voxels with non-positive or invalid S0 values are excluded (residuals set to zero).
%

[rows, cols, slices, num_echoes] = size(images);

% Initialize
residuals = zeros(size(images));
TEs = TEs(:);

% compute the residuals of each voxel
for t = 1:num_echoes
    predicted = S0_map .* exp(-TEs(t) ./ T2_map);
    actual = images(:,:,:,t);
    residuals(:,:,:,t) = (actual - predicted).^2;
end

invalid_mask = (S0_map <= 0) | isnan(S0_map);
for t = 1:num_echoes
    temp = residuals(:,:,:,t);
    temp(invalid_mask) = 0;
    residuals(:,:,:,t) = temp;
end

return;
end