function [residuals, mean_residual, RSS] = evaluate_model(T2_map, S0_map, images, TEs, running_time, seg)
% evaluate_model computes the residuals and error metrics for a fitted T2 model.
%
% Syntax:
%   [residuals, mean_residual, RSS] = evaluate_model(T2_map, S0_map, images, TEs, running_time, seg)
%
% Inputs:
%   T2_map      - 3D matrix of estimated T2 relaxation times (ms).
%   S0_map      - 3D matrix of estimated initial signal intensities (S0).
%   images      - 4D matrix of original multi-echo MRI data [X, Y, Z, TEs].
%   TEs         - Vector of echo times (ms) corresponding to the 4th dimension of `images`.
%   running_time - (Optional) Execution time of the model, not used in calculation but may be logged.
%   seg         - Segmentation map (not used directly in this function but included for completeness).
%
% Outputs:
%   residuals     - 4D matrix of squared differences between actual and predicted signals.
%   mean_residual - Mean absolute residual across all voxels (scalar).
%   RSS           - 3D matrix of residual sum of squares (RSS) per voxel, summed over echo times.
%
% Description:
%   This function calculates the residuals between the actual signal intensities
%   and the signal predicted using the estimated T2 and S0 values. It computes:
%     - voxel-wise squared residuals for each echo time;
%     - the residual sum of squares (RSS) per voxel;
%     - and the overall mean absolute residual across the volume.
%

% Initialize
RSS = zeros(96, 96, 55);

% 1. compute the residuals between the origin signal and the one
% simulated from S0
residuals = calculate_residuals(images, TEs, T2_map, S0_map);
for i = 1 : 96
    for j = 1 : 96
        for k = 1 : 55
            RSS(i, j, k) = sum(residuals(i, j, k, :));
        end
    end
end

mean_residual = mean(abs(RSS(:)));

return;
end
