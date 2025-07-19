function Analyze_all_tissues(T2_1_map, T2_2_map, S0_map, V1_map, seg)
% Analyze_all_tissues computes summary statistics (mean and 95% confidence intervals)
% for multiple MRI-derived parameters across major brain tissue types.
%
% Syntax:
%   Analyze_all_tissues(T2_1_map, T2_2_map, S0_map, V1_map, seg)
%
% Inputs:
%   T2_1_map - 3D matrix representing the T2 short component (in ms) per voxel.
%   T2_2_map - 3D matrix representing the T2 long component (in ms) per voxel.
%   S0_map   - 3D matrix representing the proton density (S0) per voxel.
%   V1_map   - 3D matrix representing the volume fraction of the short T2 component.
%   seg      - Struct containing the segmentation output, with field:
%              seg.img - 4D matrix [x, y, z, class] with probabilistic tissue maps:
%                        class 2: CSF, class 3: GM, class 4: WM.
%
% Outputs:
%   This function prints to the console a formatted table of results with the following columns:
%     - Tissue: White Matter, Grey Matter, CSF
%     - Parameter: T2 short, T2 long, S0, V1
%     - Mean: average value within the tissue
%     - 95% CI: lower and upper bounds of the confidence interval
%     - n: number of voxels included in the estimate
%
% Description:
%   The function extracts binary masks for White Matter, Grey Matter, and CSF
%   by thresholding the corresponding segmentation probability maps at > 0.99.
%   For each tissue type, it evaluates the following quantitative MRI parameters:
%     - T2_1 (short T2 relaxation time)
%     - T2_2 (long T2 relaxation time)
%     - S0   (signal amplitude)
%     - V1   (short T2 component volume fraction)
%
%   Parameter statistics are computed using the helper function:
%       calculate_parameter_estimate(param_map, mask)
%   which returns the mean, 95% confidence interval, and voxel count within the mask.
%
% Notes:
%   - Voxels are included only if their tissue probability > 0.99.
%   - Results are printed in a formatted table for easy reading.
%   - No variables are returned or saved — this function is display-only.
%




fprintf('\n===== Parameter Estimates with 95%% Confidence Intervals =====\n');
fprintf('Tissue        | Parameter     | Mean    | 95%% CI           | n\n');
fprintf('------------------------------------------------------------------\n');
%% 1. Estimate on white matter
wm_mask = seg.img(:, :, :, 4);
wm_voxels = wm_mask(:, :, :) > 0.99;

% estimate on T2 short
[mean_T2_wm,CI_lower_wm,CI_upper_wm,n_samples_wm] = calculate_parameter_estimate(T2_1_map,wm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'White Matter', 'T2 short (ms)', mean_T2_wm, CI_lower_wm, CI_upper_wm, n_samples_wm);

% estimate on T2 long
[mean_T2_wm,CI_lower_wm,CI_upper_wm,n_samples_wm] = calculate_parameter_estimate(T2_2_map,wm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'White Matter', 'T2 long (ms)', mean_T2_wm, CI_lower_wm, CI_upper_wm, n_samples_wm);

% estimate on S0
[mean_S0_wm,CI_lower_wm,CI_upper_wm,n_samples_wm] = calculate_parameter_estimate(S0_map,wm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'White Matter', 'S0', mean_S0_wm, CI_lower_wm, CI_upper_wm, n_samples_wm);

% estimate on V1
[mean_V1_wm,CI_lower_wm,CI_upper_wm,n_samples_wm] = calculate_parameter_estimate(V1_map,wm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'White Matter', 'V1', mean_V1_wm, CI_lower_wm, CI_upper_wm, n_samples_wm);

%% 2. Estimate on grey matter
gm_mask = seg.img(:, :, :, 3);
gm_voxels = gm_mask(:, :, :) > 0.99;

% estimate on T2 short
[mean_T2_gm,CI_lower_gm,CI_upper_gm,n_samples_gm] = calculate_parameter_estimate(T2_1_map,gm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'Grey Matter', 'T2 short (ms)', mean_T2_gm, CI_lower_gm, CI_upper_gm, n_samples_gm);

% estimate on T2 long
[mean_T2_gm,CI_lower_gm,CI_upper_gm,n_samples_gm] = calculate_parameter_estimate(T2_2_map,gm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'Grey Matter', 'T2 long (ms)', mean_T2_gm, CI_lower_gm, CI_upper_gm, n_samples_gm);

% estimate on S0
[mean_S0_gm,CI_lower_gm,CI_upper_gm,n_samples_gm] = calculate_parameter_estimate(S0_map,gm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'Grey Matter', 'S0', mean_S0_gm, CI_lower_gm, CI_upper_gm, n_samples_gm);

% estimate on V1
[mean_V1_gm,CI_lower_gm,CI_upper_gm,n_samples_gm] = calculate_parameter_estimate(V1_map,gm_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'Grey Matter', 'V1', mean_V1_gm, CI_lower_gm, CI_upper_gm, n_samples_gm);
%% 3. Estimate on CSF
csf_mask = seg.img(:, :, :, 2);
csf_voxels = csf_mask(:, :, :) > 0.99;

% estimate on T2 short
[mean_T2_csf,CI_lower_csf,CI_upper_csf,n_samples_csf] = calculate_parameter_estimate(T2_1_map,csf_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'CSF', 'T2 short (ms)', mean_T2_csf, CI_lower_csf, CI_upper_csf, n_samples_csf);

% estimate on T2 long
[mean_T2_csf,CI_lower_csf,CI_upper_csf,n_samples_csf] = calculate_parameter_estimate(T2_2_map,csf_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'CSF', 'T2 long (ms)', mean_T2_csf, CI_lower_csf, CI_upper_csf, n_samples_csf);

% estimate on S0
[mean_S0_csf,CI_lower_csf,CI_upper_csf,n_samples_csf] = calculate_parameter_estimate(S0_map,csf_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'CSF', 'S0', mean_S0_csf, CI_lower_csf, CI_upper_csf, n_samples_csf);

% estimate on V1
[mean_V1_csf,CI_lower_csf,CI_upper_csf,n_samples_csf] = calculate_parameter_estimate(V1_map,csf_voxels);
fprintf('%-14s| %-14s| %-8.3f| %.3f - %-8.3f| %d\n', ...
        'CSF', 'V1', mean_V1_csf, CI_lower_csf, CI_upper_csf, n_samples_csf);
end

