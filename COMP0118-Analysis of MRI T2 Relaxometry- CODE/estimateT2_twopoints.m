function [T2_map, mean_residual] = estimateT2_twopoints(image1, image2, TE1, TE2)
% estimateT2_twopoints estimates the T2 relaxation time map from two images.
%
% Syntax:
%   [T2_map, mean_residual] = estimateT2_twopoints(image1, image2, TE1, TE2)
%
% Inputs:
%   image1 - First MR image acquired at echo time TE1 (3D or 2D matrix)
%   image2 - Second MR image acquired at echo time TE2 (same dimensions as image1)
%   TE1    - Echo time corresponding to image1 (in ms)
%   TE2    - Echo time corresponding to image2 (in ms)
%
% Outputs:
%   T2_map        - Estimated T2 relaxation time map (same size as input images)
%   mean_residual - Mean squared residual between actual and predicted signal
%
% Description:
%   This function estimates the voxel-wise T2 relaxation times assuming
%   mono-exponential signal decay:
%       S(TE) = S0 * exp(-TE / T2)
%   It uses two images acquired at different echo times to solve for T2.
%   Only voxels with positive signal in both images are considered valid.
%
%   The function also computes the mean residual error between the
%   actual signal at TE2 and the predicted signal at TE2 using the T2 estimate.
%

% Ensure input dimensions match
if ~isequal(size(image1), size(image2))
    error('Input images must have the same dimensions');
end

% Initialize
T2_map = zeros(size(image1));

% Define valid voxels
valid_mask = (image1 > 0) & (image2 > 0);

% Estimate T2
T2_map(valid_mask) = (TE2 - TE1) ./ log(image1(valid_mask) ./ image2(valid_mask));

% Remove invalid or negative T2 values
invalid_T2 = isnan(T2_map) | isinf(T2_map) | (T2_map < 0);
T2_map(invalid_T2) = 0;
valid_mask(invalid_T2) = false;  % Also exclude them from residual calc

% Reconstruct signal2 from signal1 and estimated T2
S1 = image1(valid_mask);
T2 = T2_map(valid_mask);
predicted_S2 = S1 .* exp(-(TE2 - TE1) ./ T2);
actual_S2 = image2(valid_mask);

% Calculate Residuals
voxel_residuals = (actual_S2 - predicted_S2).^2;

mean_residual = mean(abs(voxel_residuals));
end
