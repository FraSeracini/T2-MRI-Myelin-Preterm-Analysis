function AIC_C = Compute_AIC_Corrected(RSS, n_parameters, num_echoes)
% Compute_AIC_Corrected calculates the corrected Akaike Information Criterion (AICc) for small sample sizes.
%
% Syntax:
%   AIC_C = Compute_AIC_Corrected(RSS, n_parameters, num_echoes)
%
% Inputs:
%   RSS          - Residual Sum of Squares from the model fit.
%   n_parameters - Number of parameters in the model (excluding the noise term).
%   num_echoes   - Number of data points (echo times) used in the fit.
%
% Outputs:
%   AIC_C - Corrected Akaike Information Criterion (AICc) value.
%
% Description:
%   This function computes the corrected version of the Akaike Information Criterion (AICc),
%   which adjusts the AIC for small sample sizes to avoid overfitting:
%
%       AIC = num_echoes * log(RSS / num_echoes) + 2 * (n_parameters + 1)
%       AICc = AIC + [2 * (k) * (k + 1)] / (n - k - 1)
%
%   where:
%     - k = number of parameters including the variance term (n_parameters + 1),
%     - n = number of data points (echoes),
%     - RSS = residual sum of squares from model fitting.
%
%   AICc is particularly useful when the number of observations is small relative to the number of parameters.
%

n_parameters = n_parameters + 1;
AIC = Compute_AIC_T2(n_parameters, num_echoes, RSS);
AIC_C = AIC + (2 * n_parameters * (n_parameters + 1)) / (num_echoes - n_parameters -1);

end