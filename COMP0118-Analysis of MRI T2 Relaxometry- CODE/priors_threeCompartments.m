function [S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map, mean_residual, RSS] = priors_threeCompartments(images, TEs, mask, seg)
% priors_threeCompartments fits a three-compartment T2 relaxation model to 
% multi-echo MRI data using prior segmentation information.
%
% Syntax:
%   [S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map, mean_residual, RSS] = ...
%       priors_threeCompartments(images, TEs, mask, seg)
%
% Inputs:
%   images - 4D matrix [rows, cols, slices, echoes] of multi-echo MRI data.
%   TEs    - Vector of echo times (in milliseconds), one per echo.
%   mask   - 3D binary mask defining the voxels to include in the fitting.
%   seg    - 4D segmentation probability map [x, y, z, class], where:
%            seg(:,:,:,2) = CSF
%            seg(:,:,:,3) = Grey Matter
%            seg(:,:,:,4) = White Matter
%            seg(:,:,:,5) = Myelin (if available)
%
% Outputs:
%   S0_map    - Estimated proton density (S0) map.
%   V1_map    - Volume fraction map of the fast-relaxing (myelin-like) T2 component.
%   V2_map    - Volume fraction map of the intermediate T2 component (WM/GM).
%   V3_map    - Volume fraction map of the long T2 component (CSF).
%   T2_1_map  - Estimated short T2 component map (in ms).
%   T2_2_map  - Estimated intermediate T2 component map (in ms).
%   T2_3_map  - Estimated long T2 component map (in ms).
%   mean_residual - Mean residual error across all fitted voxels.
%   RSS       - Residual sum of squares (3D map), computed per voxel.
%
% Description:
%   This function fits a signal decay model of the form:
%       S(TE) = S0 * [v1 * exp(-TE/T2_1) + v2 * exp(-TE/T2_2) + v3 * exp(-TE/T2_3)]
%   for each voxel in the input mask using non-linear least squares (`lsqcurvefit`).
%
%   The volume fractions v1, v2, and v3 are initialized using the segmentation priors
%   (e.g., myelin, WM/GM, CSF) and normalized to sum to 1. These priors guide the fitting
%   process and are constrained via equality constraints.
%
%   The initial guesses and bounds for all parameters are set based on typical tissue
%   values (e.g., T2_myelin ≈ 20 ms, T2_WM/GM ≈ 80 ms, T2_CSF ≈ 2000 ms).
%   Residuals are computed at each TE, and the voxel-wise residual sum of squares (RSS)
%   is calculated and returned.
%
% Notes:
%   - Fitting is performed slice-wise (currently only slice 28).
%   - Voxels with non-positive signal or outside the mask are ignored.
%   - Output maps are cleaned of NaN and Inf values at the end.
%   - Optimization includes both bound constraints and a linear equality constraint: v1 + v2 + v3 = 1.
%   - Residuals are squared differences between actual and predicted signals for each TE.
%

[rows, cols, slices, num_echoes] = size(images);
slices = 1;
residuals = zeros(rows, cols, length(TEs), 'double');
RSS = zeros(rows, cols);

% Output maps
S0_map   = zeros(rows, cols, 'double');
V1_map   = zeros(rows, cols, 'double');
V2_map   = zeros(rows, cols, 'double');
V3_map   = zeros(rows, cols, 'double');
T2_1_map = zeros(rows, cols, 'double');
T2_2_map = zeros(rows, cols, 'double');
T2_3_map = zeros(rows, cols, 'double');

TEs = double(TEs);
delta = 0.2; % margine per i vincoli sui volumi frazionari

parfor i = 1:rows
    for j = 1:cols
        for k = 28:28
            signal = double(squeeze(images(i, j, k, :)));

            if all(signal > 0) && mask(i, j, k) > 0

                % Normalizza i canali di segmentazione in [0,1]
                seg_csf = seg(i,j,k,2);
                seg_gm  = seg(i,j,k,3);
                seg_wm  = seg(i,j,k,4);
                seg_mye = seg(i,j,k,5);  % possibile comparto mielinico

                % % Segnali utili
                % wm = mSeg(i, j, k, 4) + mSeg(i, j, k, 5);
                % % wm = mSeg(i, j, k, 4)
                % gm = mSeg(i, j, k, 3);
                % csf = mSeg(i, j, k, 2);
                %
                % % Definizione della myelin water come percentuale della WM (fissiamo 25%)
                % myelin_fraction = 0.25;
                % v1 = myelin_fraction * wm;
                % v2 = gm + (1 - myelin_fraction) * wm;
                % v3 = csf;


                options = optimoptions('lsqcurvefit', ...
                    'Display', 'off', ...
                    'MaxFunctionEvaluations', 3000, ...
                    'MaxIterations', 2500, ...
                    'FunctionTolerance', 1e-6, ...
                    'StepTolerance', 1e-6);

                % if seg_csf > 0.99
                % 
                %     % Inizializzazione
                %     x0 = double([max(signal), 2000]);  % [S0, T2]
                %     lb = double([0, 1800]);             % lower bounds
                %     ub = double([Inf, 2500]);         % upper bounds
                % 
                %     % Define the objective function
                %     model_fun = @(x, TE) x(1) * exp(-TE / x(2));
                % 
                %     try
                %         params = lsqcurvefit(model_fun, x0,  TEs(:), signal(:), lb, ub, [], [], [], [], [], options);
                % 
                %         S0_map(i, j)   = params(1);
                %         V1_map(i, j)   = 0;
                %         V2_map(i, j)   = 0;
                %         V3_map(i, j)   = 1;
                %         T2_1_map(i, j) = 20;
                %         T2_2_map(i, j) = 80;
                %         T2_3_map(i, j) = params(2);
                % 
                %     catch
                %         continue
                %     end

                % elseif (seg_wm + seg_gm + seg_mye) > 0.95
                % 
                %     v1 = seg_mye / (seg_wm + seg_mye+ seg_gm);
                %     v2 = (seg_wm + seg_gm) / (seg_wm + seg_mye+ seg_gm);
                % 
                %     % Inizializzazione: [S0, v1, v2, T2_1, T2_2]
                %     x0 = double([max(signal), v1, v2, 20, 80]);
                % 
                %     % Vincoli: v1 e v2 devono essere in [0,1], T2 positivi, S0 positivo
                %     lb = double([0, 0, 0, 10, 50]);
                %     ub = double([Inf, 1, 1, 50, 110]);
                % 
                %     % Vincolo lineare: v1 + v2 = 1
                %     Aeq = [0, 1, 1, 0, 0];
                %     beq = 1;
                % 
                %     % Modello
                %     model_fun = @(x, TE) x(1) * ( ...
                %         x(2)*exp(-TE / x(4)) + ...
                %         x(3)*exp(-TE / x(5)) );
                % 
                %     try
                % 
                %         params = lsqcurvefit(model_fun, x0, TEs(:), signal(:), lb, ub, [], [], Aeq, beq, [], options);
                % 
                %         S0_map(i, j)   = params(1);
                %         V1_map(i, j)   = params(2);
                %         V2_map(i, j)   = params(3);
                %         V3_map(i, j)   = 0;
                %         T2_1_map(i, j) = params(4);
                %         T2_2_map(i, j) = params(5);
                %         T2_3_map(i, j) = 2000;
                % 
                %     catch
                %         continue
                %     end

                % else

                    % Priors (inizializzazioni)
                    % v1 = seg_mye + 0.25 * seg_wm;      % prior frazione mielina
                    v1 = seg_mye;
                    % v1 = 0.25 * seg_wm;
                    % v2 = 0.75 * seg_wm + seg_gm;  % prior per comparto medio (~80ms)
                    v2 = seg_wm + seg_gm;
                    v3 = seg_csf;         % prior per comparto lento (~2000ms)

                    % Normalizzazione
                    total = v1 + v2 + v3;
                    if total > 0
                        v1 = v1 / total;
                        v2 = v2 / total;
                        v3 = v3 / total;
                    else
                        v1 = 0.2; v2 = 0.6; v3 = 0.2;
                    end

                    % Inizializzazione
                    x0 = [max(signal), v1, v2, v3, 20, 80, 2000];
                    % x0 = [max(signal), 0.3, 0.5, 0.2, 20, 80, 2000];
                    x0 = double(x0);

                    % Vincoli su frazioni (clamp a [0,1] per sicurezza)
                    % lb = double([0,   max(0,v1-delta), max(0,v2-delta), max(0,v3-delta),   5,  65, 1800]);
                    % ub = double([Inf, min(1,v1+delta), min(1,v2+delta), min(1,v3+delta),  60, 100, 3000]);
                    % ub = double([Inf, 1, 1, min(1,v3+delta), 100 , 200, 2500]);
                    % lb = double([0,  0, 0, 0,  0,  0, 0]);
                    % ub = double([Inf, 1 , 1 , 1,  50, 110, 3000]);
                    lb = [0,     0,   0,   0,   0,   65,  2000];
                    ub = [Inf,   1,   1,   1, 40, 200, inf];


                    %
                    Aeq = [0 1 1 1 0 0 0];
                    beq = 1;

                    model_fun = @(x, TE) x(1) * ( ...
                        x(2)*exp(-TE/x(5)) + ...
                        x(3)*exp(-TE/x(6)) + ...
                        x(4)*exp(-TE/x(7)) );

                    try
                        params = lsqcurvefit(model_fun, x0, TEs(:), signal(:), lb, ub, [], [], Aeq, beq, [], options);

                        S0_map(i, j)   = params(1);
                        V1_map(i, j)   = params(2);
                        V2_map(i, j)   = params(3);
                        V3_map(i, j)   = params(4);
                        T2_1_map(i, j) = params(5);
                        T2_2_map(i, j) = params(6);
                        T2_3_map(i, j) = params(7);

                    catch
                        continue
                    end
                % end
            end
        end
    end
end

% Clean-up NaNs e Infs
maps = {S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map};
for m = 1:numel(maps)
    maps{m}(isnan(maps{m}) | isinf(maps{m})) = 0;
end
[S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map] = deal(maps{:});

% Calcolo dei residui
for t = 1:num_echoes
    predicted = S0_map .* ( ...
        V1_map .* exp(-TEs(t) ./ T2_1_map) + ...
        V2_map .* exp(-TEs(t) ./ T2_2_map) + ...
        V3_map .* exp(-TEs(t) ./ T2_3_map));
    actual = images(:,:,28,t);
    residuals(:,:, t) = (actual - predicted).^2;
end

invalid_mask = (S0_map <= 0) | isnan(S0_map);
for t = 1:num_echoes
    temp = residuals(:,:,t);
    temp(invalid_mask) = 0;
    residuals(:,:,t) = temp;
end

for i = 1:rows
    for j = 1:cols
        RSS(i, j) = sum(residuals(i, j, :));
    end
end

mean_residual = mean(RSS(:), 'omitnan');

end
