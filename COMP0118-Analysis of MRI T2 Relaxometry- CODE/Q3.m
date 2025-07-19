%% Q3
% Extend the models for fitting two compartments and use the segmentations
% provided to find the averages of the parameters for WM, GM, CSF

%% Load the data
r=load_nii('case01-qt2_reg.nii');
r.img(r.img(:)<0)=0;
images = r.img;
brain_mask = load_nii('case01-mask.nii');
TEs=load('case01-TEs.txt');

% Load the segmentation file
seg = load_nii('case01-seg.nii');

% Choose a central slice to visualize
slice_num = round(size(images, 3) / 2);

%% Two Compartments, model: NLLS with estimated T2

% Estimate T2 maps using a bi-component NLLS model
[T2_1_map, T2_2_map, S0_map, V1_map, residuals, mean_residual, RSS] = estimateT2_multipoint_NLLS_bicomponent(images, TEs, brain_mask.img);
mean_residual

save('T2_NLLS_bicomponent_results.mat', ...
     'T2_1_map', 'T2_2_map', 'S0_map', 'V1_map', 'residuals', 'mean_residual', 'RSS');


% Display T2_1 map (short T2 component)
figure;
imagesc(rot90(flipud(T2_1_map(:,:,slice_num)))); 
% title('Estimated Short T2 (T2_1) with NLLS');
colorbar;
axis image;
clim([0 100]);  % Optional: set color limits

% Display T2_2 map (long T2 component)
figure;
imagesc(rot90(flipud(T2_2_map(:,:,slice_num)))); 
% title('Estimated Long T2 (T2_2) with NLLS');
colorbar;
axis image;
clim([0 2500]); 


% Display V1 map (volume fraction of compartment 1)
figure;
imagesc(rot90(flipud(V1_map(:,:,slice_num)))); 
title('Compartment 1 Volume Fraction (V1) with NLLS');
colorbar;
axis image;
clim([0 1]);  

% Display S0 map (initial signal intensity)
figure;
imagesc(rot90(flipud(S0_map(:,:,slice_num)))); 
title('Estimated S0 with NLLS');
colorbar;
axis image;


fprintf('\n===== Fit with NLLS =====\n');
Analyze_all_tissues(T2_1_map, T2_2_map, S0_map, V1_map, seg);

%% Two Compartments, model: NLLS with fixed T2
[T2_1_map, T2_2_map, S0_map, V1_map, residuals, mean_residual, RSS] = estimateT2_fixedT2s(images, TEs, brain_mask.img);
mean_residual

save('T2_fixed_bicomponent_results.mat', ...
     'T2_1_map', 'T2_2_map', 'S0_map', 'V1_map', 'residuals', 'mean_residual', 'RSS');

% Display V1 map (volume fraction of compartment 1)
figure;
imagesc(rot90(flipud(V1_map(:,:,slice_num)))); 
title('Compartment 1 Volume Fraction (V1) with fixed T2');
colorbar;
axis image;
clim([0 1]);  

% Display S0 map (initial signal intensity)
figure;
imagesc(rot90(flipud(S0_map(:,:,slice_num)))); 
title('Estimated S0 with fixed T2');
colorbar;
axis image;

fprintf('\n===== Fit with NLLS(fixedT2) =====\n');
Analyze_all_tissues(T2_1_map, T2_2_map, S0_map, V1_map, seg);

%% Two Compartments, model: NNLS with fixed T2
[S0_map, V1_map, V2_map, residuals, mean_residual, RSS] = estimateT2_NNLS_2comp(images, TEs, brain_mask.img);
mean_residual

save('T2_NNLS_bicomponent_results.mat', ...
      'S0_map', 'V1_map', 'residuals', 'mean_residual', 'RSS');

% Display V1 map (volume fraction of compartment 1 - T2 = 60ms)
figure;
imagesc(rot90(flipud(V1_map(:,:,slice_num)))); 
title('Compartment 1 Volume Fraction (T2 = 60 ms)');
colorbar;
axis image;
clim([0 1]);  

% Display V2 map (volume fraction of compartment 2 - T2 = 2000ms)
figure;
imagesc(rot90(flipud(V2_map(:,:,slice_num)))); 
title('Compartment 2 Volume Fraction (T2 = 2000 ms)');
colorbar;
axis image;
clim([0 1]);

% Display S0 map (initial signal intensity)
figure;
imagesc(rot90(flipud(S0_map(:,:,slice_num)))); 
title('Estimated S0 with fixed T2 (2-component NNLS)');
colorbar;
axis image;

[rows, cols, slices, num_echoes] = size(images);
T2_1_map = 60 .* ones(rows, cols, slices, 'double');
T2_2_map = 200 .* ones(rows, cols, slices, 'double');

fprintf('\n===== Fit with NNLS(fixedT2) =====\n');
Analyze_all_tissues(T2_1_map, T2_2_map, S0_map, V1_map, seg);