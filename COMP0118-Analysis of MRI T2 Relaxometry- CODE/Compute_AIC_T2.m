function AIC = Compute_AIC_T2(n_parameters, num_echoes, RSS)
% Compute_AIC_T2 calculates the Akaike Information Criterion (AIC) for T2 model fitting.
%
% Syntax:
%   AIC = Compute_AIC_T2(n_parameters, num_echoes, RSS)
%
% Inputs:
%   n_parameters - Number of model parameters used in the T2 fitting.
%   num_echoes   - Number of data points (echo times) used in the fitting.
%   RSS          - Residual Sum of Squares from the model fit.
%
% Outputs:
%   AIC - Computed Akaike Information Criterion value (scalar).
%
% Description:
%   This function computes the AIC, a model selection metric that balances model goodness-of-fit
%   (via the Residual Sum of Squares, RSS) with model complexity (number of parameters):
%
%       AIC = num_echoes * log(RSS / num_echoes) + 2 * n_parameters
%
%   A lower AIC value indicates a more optimal balance between accuracy and simplicity.
%   This criterion can be used to compare different T2 models (e.g., mono-, bi-, or multi-exponential).
%

% Compute AIC
AIC = num_echoes * log(RSS / num_echoes) + 2 * n_parameters;

end