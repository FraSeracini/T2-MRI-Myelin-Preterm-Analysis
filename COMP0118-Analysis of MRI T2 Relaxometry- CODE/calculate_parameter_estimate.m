function [mean_value,CI_lower,CI_upper,n_samples] = calculate_parameter_estimate(map,mask)
% calculate_parameter_estimate computes the mean and 95% confidence interval
% for a given quantitative map within a binary mask.
%
% Syntax:
%   [mean_value, CI_lower, CI_upper, n_samples] = calculate_parameter_estimate(map, mask)
%
% Inputs:
%   map  - 3D matrix (or flattened vector) containing quantitative parameter values
%          (e.g., T2 map, S0 map, V1 map).
%
%   mask - 3D logical or binary mask of the same size as map, indicating
%          which voxels to include in the estimate (true = include).
%
% Outputs:
%   mean_value - Mean value of the parameter within the valid masked voxels.
%   CI_lower   - Lower bound of the 95% confidence interval for the mean.
%   CI_upper   - Upper bound of the 95% confidence interval for the mean.
%   n_samples  - Total number of valid voxels used for the estimation.
%
% Description:
%   This function extracts all values from the input map where the mask is true,
%   filters out invalid entries (e.g., non-finite or non-positive values),
%   and computes:
%     - the mean value
%     - the standard error of the mean (SEM)
%     - the 95% confidence interval based on the Student's t-distribution
%
%   The confidence interval is computed as:
%       CI = mean ± t * (std / sqrt(n))
%   where `t` is the critical value of the t-distribution at 97.5% confidence
%   for (n - 1) degrees of freedom.
%
% Notes:
%   - Only finite and strictly positive values are included in the estimate.
%   - If no valid voxels remain after filtering, the output may return NaNs.
%   - The function assumes parametric data and normality for CI estimation.
%

% estimate on the parameters
para_values = map(mask);
para_values = para_values(isfinite(para_values) & para_values > 0);
mean_value = mean(para_values);
std_value = std(para_values);
n_samples = length(para_values);

% compute the confidence intervals
SE = std_value / sqrt(n_samples);
t_critical = tinv(0.975, n_samples-1);
CI_lower = mean_value - t_critical * SE;
CI_upper = mean_value + t_critical * SE;
end

