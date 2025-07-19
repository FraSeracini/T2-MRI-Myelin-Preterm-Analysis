%% Q2
% We develop and compare different one-compartment models for the estimate
% of the S0 and T2 parameters

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

%% Two-point algorithm

% Select the two images
id_img1 = 1;
id_img2 = 19;

% Estimate the parameters
tic;
[T2_2point, mean_residual] = estimateT2_twopoints(img_masked(:,:,:,id_img1), img_masked(:,:,:,id_img2), TEs(id_img1), TEs(id_img2));
time_2point = toc;

% Visualize the T2_map
figure;
imagesc(rot90(flipud(T2_2point(:,:,slice_num))));
title('Estimated T2 (ms)');
colorbar;
axis image;
clim([10 2500]);

%% Linear Least Square

% Estimate the parameters
tic;
[T2_linear, S0_linear] = estimateT2_multipoint_linear(img_masked, TEs);
time_linear = toc;

% Calculate the error
[residuals,mean_residual, RSS] = evaluate_model(T2_linear, S0_linear, img_masked, TEs, time_linear, seg);
time_linear
mean_residual

% Visualize the results
figure;
subplot(1,2,1);
imagesc(rot90(flipud(T2_linear(:,:,slice_num))));
title('T2 by Linear LS'); 
colorbar; 
clim([0 2500]);

subplot(1,2,2);
imagesc(rot90(flipud(S0_linear(:,:,slice_num))));
colorbar;
title('S0 Map by Linear LS');

%% Weighted least square

% Estimate the parameters
tic;
[T2_weighted, S0_weighted] = estimateT2_multipoint_weighted(img_masked, TEs);
time_weighted = toc;
[residuals, mean_residual, RSS] = evaluate_model(T2_weighted, S0_weighted, img_masked, TEs, time_weighted, seg)

mean_residual
time_weighted

% visualize the results
figure;
subplot(1,2,1);
imagesc(rot90(flipud(T2_weighted(:,:,slice_num))));title('T2 by Weighted LS'); 
colorbar; 
clim([0 2500]);

subplot(1,2,2);
imagesc(rot90(flipud(S0_weighted(:,:,slice_num))));colorbar;
title('S0 Map by Weighted LS');

%% Non-negative least-squares

% Estimate the parameters
tic;
[T2_NNLS, S0_NNLS] = estimateT2_multipoint_NNLS(img_masked, TEs);
time_NNLS = toc;
[residuals, mean_residual, RSS] = evaluate_model(T2_NNLS, S0_NNLS, img_masked, TEs, time_NNLS, seg);

mean_residual
time_NNLS

% Visualize the results
figure;
subplot(1,2,1);
imagesc(rot90(flipud(T2_NNLS(:,:,slice_num))));title('T2 by NNLS'); 
colorbar; 
clim([0 2500]);

subplot(1,2,2);
imagesc(rot90(flipud(S0_NNLS(:,:,slice_num))));colorbar;
title('S0 Map by NNLS');

%% Non-linear least squares

% Estimate the parameters
tic;
[T2_NLLS, S0_NLLS] = estimateT2_multipoint_NLLS(img_masked, TEs);
time_NLLS = toc;
[residuals, mean_residual_NLLS, RSS] = evaluate_model(T2_NLLS, S0_NLLS, img_masked, TEs, time_NLLS, seg);

mean_residual_NLLS
time_NLLS

save('T2_NLLS_one_comp.mat', ...
 'T2_NLLS','S0_NLLS', 'residuals', 'mean_residual_NLLS', 'RSS');

figure;
subplot(1,2,1);
imagesc(rot90(flipud(T2_NLLS(:,:,slice_num))));title('T2 by NLLS'); 
colorbar; 
clim([0 2500]);

subplot(1,2,2);
imagesc(rot90(flipud(S0_NLLS(:,:,slice_num))));colorbar;
title('S0 Map by NLLS');

%% 

% 5. non-linear least squares with fmincon
% tic;
% [T2_fmincon, S0_fmincon] = estimateT2_multipoint_fmincon(img_masked, TEs);
% time_fmincon = toc;
% [residuals, mean_residual_fmincon, RSS] = evaluate_model(T2_fmincon, S0_fmincon, img_masked, TEs, time_fmincon, seg);
% 
% save('T2_fmincon_one_comp.mat', ...
%      'T2_NLLS','S0_NLLS', 'residuals', 'mean_residual_fmincon', 'RSS');
% 
% figure;
% subplot(1,2,1);
% imagesc(rot90(fliplr(T2_NLLS(:,:,slice_num))));title('T2 by fmincon'); 
% colorbar; 
% clim([0 2500]);
% 
% subplot(1,2,2);
% imagesc(rot90(fliplr(S0_NLLS(:,:,slice_num))));colorbar;
% title('S0 Map by fmincon');

%% Visualize all the results together

% Print the fit time
fprintf('Least Square method: %.2f seconds\n', time_linear);
fprintf('Weighted least square method: %.2f seconds\n', time_weighted);
fprintf('non-negative least-squares method: %.2f seconds\n', time_NNLS);
fprintf('non-linear least squares method: %.2f seconds\n', time_NLLS);
% fprintf('non-linear least squares method with fmincon: %.2f seconds\n', time_fmincon);

% Print and compare the errors of the models
fprintf('Least Square method error: %.2f\n', error_linear);
fprintf('Weighted least square method error: %.2f\n',error_weighted);
fprintf('non-negative least-squares error: %.2f\n', error_NNLS);
fprintf('non-linear least squares error: %.2f\n', error_residual_NLLS);
% fprintf('non-linear least squares score with fmincon: %.2f\n', mean_residual_fmincon);