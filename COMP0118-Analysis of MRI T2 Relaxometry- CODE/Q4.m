%% Q4
% Classical Bootstrap procdure to analyze the uncertainty of the parameters
% and estimation of AIC corrected for different models.
% Estimate the noise std.

%% Load the data

% read in the images
r=load_nii('case01-qt2_reg.nii');
r.img(r.img(:)<0)=0;
images = r.img;

TEs=load('case01-TEs.txt');
% disp(TEs);
brain_mask = load_nii('case01-mask.nii');
img_masked = images.*brain_mask.img;

% Select a slice
slice_num = round(size(images,3)/2); % A central slice

% Load the segmentation file
seg = load_nii('case01-seg.nii');

% Load values

load('T2_NNLS_bicomponent_results.mat','RSS');
RSS_NNLS = RSS;

load('T2_NLLS_bicomponent_results.mat', 'RSS');
RSS_NLLS = RSS;

load('T2_NLLS_one_comp','RSS');
RSS_NLLS_one_comp = RSS;

% Number of parameters in each model
n_parameters_NNLS = 2;
n_parameters_NLLS = 4;
n_parameters_NLLS_one_comp = 2;

%% Parametric Bootstrap

i = 45;
j = 45;
k = 28;

% Select a voxel
signal = double(squeeze(images(i, j, k, :)));

if all(signal > 0)
    if(brain_mask.img(i, j, k) > 0)

        % Starting point: [S0, v1, T2_1, T2_2]
        p0 = [max(signal), 0.5, 20, 70];
        lb = [0, 0, 10, 110];
        ub = [inf, 1, 110, 2500];

        % Non-linear model: S(TE) = S0 * [v1 * exp(-TE/T2_1) + (1 - v1) * exp(-TE/T2_2)]
        objective_fun = @(p, TE) p(1) * (p(2) * exp(-TE / p(3)) + (1 - p(2)) * exp(-TE / p(4)));

        % Fit with non-linear least squares
        options = optimset('Display', 'off', 'MaxFunEvals', 1000);


        params = lsqcurvefit(objective_fun, p0, TEs(:), signal(:), lb, ub, options);

    end
end


% Choose the number of iterations
T = 1000;

% Extract the estimated signal with the parameters obtained
signal_est = params(1) * (params(2) * exp(-TEs/params(3)) + ...
    (1 - params(2)) * exp(-TEs / params(4)));

% Compute the standard deviation
sigma = sqrt((1 / (length(TEs) - length(params))) * sum(signal - (signal_est)').^2);

% Initialize the matrix where we are going to store the samples
total_samples = zeros(length(params), T);

% Apply the parametric bootstrap algorithm
parfor t = 1 : T

    % Create samples from noise distribution N(0, sigma)
    E_samples = normrnd(0, sigma, 1, length(TEs));

    % Synthesize bootstrap data set
    A_hat = ((signal_est) + E_samples)';
    
      model_fun = @(x, TE) x(1) * ( ...
        x(2)*exp(-TE/x(3)) + ...
        (1 - x(2))*exp(-TE/x(4)));


    % Now run the fitting
    % try
         params = lsqcurvefit(model_fun, p0, TEs(:), A_hat(:), lb, ub, [], [], [], [], [], options);
    % catch
    %     continue;
    % end

    % Store the values
    total_samples(:, t) = params;
end

% Initialize the matrix where we are going to store the 95% range
ninety_five_range = zeros(4, 2);

% Sort in ascending order the samples
sorted_samples = sort(total_samples, 2);

for i= 1: 4
    sorted_samples_no_zeros = sorted_samples(i, sorted_samples(i, :) ~= 0);
    ninety_five_range(i, 1) = sorted_samples(i, round(length(sorted_samples_no_zeros) * 0.025));
    ninety_five_range(i, 2) = sorted_samples(i, round(length(sorted_samples_no_zeros) * (1 - 0.025)));

end


%% AIC maps

% Initialize the colour maps
colour_map_AICc = zeros(96, 96, 3);

% For each voxel select the model with the lowest AIC
for i = 1 : 96
    for j = 1 : 96
        if(brain_mask.img(i, j, slice_num) > 0)
            if(min(images(i, j, slice_num,:) > 0))

                AICc_NNLS(i, j) = Compute_AIC_Corrected(RSS_NNLS(i, j, slice_num), n_parameters_NNLS, length(TEs));
                AICc_NLLS(i, j) = Compute_AIC_Corrected(RSS_NLLS(i, j, slice_num), n_parameters_NLLS, length(TEs));

                if AICc_NNLS(i, j) < AICc_NLLS(i, j)
                    colour_map_AICc(i, j,:) = [0 0 1];
                elseif AICc_NNLS(i, j) > AICc_NLLS(i, j)
                    colour_map_AICc(i, j, :) = [1 0 0];
                end
            end
        end
    end
end

figure;
imagesc(rot90(flipud(colour_map_AICc))); 
title('AICc map');
axis image;



% %% 
% 
% map_AIC_total = zeros(96, 96, 55, 1);
% 
% for i = 1 : 96
%     for j = 1 : 96
%         for k = 1 : 55
%             if(brain_mask.img(i, j, k) > 0)
%                 if(min(images(i, j, k,:) > 0))
% 
%                     AIC_NNLS(i, j, k) = Compute_AIC_Corrected(RSS_NNLS(i, j, k), n_parameters_NNLS, length(TEs));
%                     AIC_NLLS(i, j, k) = Compute_AIC_Corrected(RSS_NLLS(i, j, k), n_parameters_NLLS, length(TEs));
% 
%                     if AIC_NNLS(i, j, k) < AIC_NLLS(i, j, k)
%                       map_AIC_total(i, j, k, :) = 1;
%                     elseif AIC_NNLS(i, j, k) > AIC_NLLS(i, j, k)
%                       map_AIC_total(i, j, k, :) = 2;
%                     end
%                 end
%             end
%         end
%     end
% end
% 

%% Estimate the noise std

% Initialize array for noise standard deviation at each TE
num_TEs = size(images, 4);
noise_std_per_TE = zeros(num_TEs, 1);

% Maximum intensity threshold to consider for noise estimation
intensity_threshold = 100;

% Loop over each TE image
for t = 1:num_TEs
    img = images(:,:,:,t);
    
    % Extract only background voxels with values below the threshold
    background_voxels = double(img(background_mask > 0));
    valid_voxels = background_voxels(background_voxels < intensity_threshold);
    
    % Compute the standard deviation of noise
    noise_std_per_TE(t) = std(valid_voxels);
end

noise_std = mean(noise_std_per_TE(:))
