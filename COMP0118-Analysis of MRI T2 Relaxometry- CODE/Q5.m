%% Q5
% Extend previous models to multiple comaprtments and compare their
% performance. 
% Use AICc to choose the best model.

%% read in the images
r=load_nii('case01-qt2_reg.nii');
r.img(r.img(:)<0)=0;
images = r.img;
brain_mask = load_nii('case01-mask.nii');
TEs=load('case01-TEs.txt');

% Load the segmentation file
seg = load_nii('case01-seg.nii');


%% 3 compartments NLLS

% Fit across the brain
[S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map, mean_residual, RSS_NLLS_3comp] = ...
    threeCompartments_estimateT2s(images, TEs, brain_mask.img);

% Save results
save('T2_NLLS_3comp_results.mat', ...
     'S0_map', 'V1_map', 'V2_map', 'V3_map', 'T2_1_map', 'T2_2_map', 'T2_3_map', 'mean_residual', 'RSS_NLLS_3comp');

disp(['Mean residual: ', num2str(mean_residual)]);

% Visualize results
slice_num = round(size(images, 3) / 2); 
slice_num = 28;

% Volume Fractions
figure;
imagesc(rot90(flipud(V1_map(:,:,slice_num))));
title('Volume Fraction (V1)');
colorbar;
axis image;
clim([0 1]);

figure;
imagesc(rot90(flipud(V2_map(:,:,slice_num))));
title('Volume Fraction (V2)');
colorbar;
axis image;
clim([0 1]);

figure;
imagesc(rot90(flipud(V3_map(:,:,slice_num))));
title('Volume Fraction (V3)');
colorbar;
axis image;
clim([0 1]);

% S0
figure;
imagesc(rot90(flipud(S0_map(:,:,slice_num))));
title('Estimated S0');
colorbar;
axis image;

% T2 maps
figure;
imagesc(rot90(flipud(T2_1_map(:,:,slice_num))));
title('Estimated T2_1 (ms)');
colorbar;
axis image;
clim([0 150]);  

figure;
imagesc(rot90(flipud(T2_2_map(:,:,slice_num))));
title('Estimated T2_2 (ms)');
colorbar;
axis image;
clim([0 300]);  

figure;
imagesc(rot90(flipud(T2_3_map(:,:,slice_num))));
title('Estimated T2_3 (ms)');
colorbar;
axis image;
clim([0 3000]); 

%% 3 compartments NNLS
% Fit across the brain
[S0_map, V1_map, V2_map, V3_map, mean_residual, RSS_NNLS_3comp] = ...
    estimateT2_NNLS_3comp(images, TEs, brain_mask.img);

% Save results
save('T2_NNLS_3comp_results.mat', ...
     'S0_map', 'V1_map', 'V2_map', 'V3_map', 'mean_residual', 'RSS_NNLS_3comp');


% Visualize results
slice_num = round(size(images, 3) / 2);  % slice centrale

% Volume Fractions
figure;
imagesc(rot90(flipud(V1_map(:,:,slice_num))));
title('Volume Fraction (V1) - T2 = 20 ms');
colorbar;
axis image;
clim([0 1]);

figure;
imagesc(rot90(flipud(V2_map(:,:,slice_num))));
title('Volume Fraction (V2) - T2 = 80 ms');
colorbar;
axis image;
clim([0 1]);

figure;
imagesc(rot90(flipud(V3_map(:,:,slice_num))));
title('Volume Fraction (V3) - T2 = 2000 ms');
colorbar;
axis image;
clim([0 1]);

% S0
figure;
imagesc(rot90(flipud(S0_map(:,:,slice_num))));
title('Estimated S0 (NNLS)');
colorbar;
axis image;

%% 4 compartments NLLS

% Fit across the brain
[S0_map, V1_map, V2_map, V3_map, V4_map, ...
 T2_1_map, T2_2_map, T2_3_map, T2_4_map, mean_residual, RSS_NLLS_4comp] = ...
    fourCompartments_estimateT2s(images, TEs, brain_mask.img);

% Save results
save('T2_NLLS_4comp_results.mat', ...
     'S0_map', 'V1_map', 'V2_map', 'V3_map', 'V4_map', 'T2_1_map', 'T2_2_map', 'T2_3_map', 'T2_4_map', 'mean_residual', 'RSS_NLLS_4comp');


% Visualize results
slice_num = round(size(images, 3) / 2); 

% Volume Fractions
figure;
imagesc(rot90(flipud(V1_map(:,:,slice_num))));
title('Volume Fraction (V1)');
colorbar;
axis image;
clim([0 1]);

figure;
imagesc(rot90(flipud(V2_map(:,:,slice_num))));
title('Volume Fraction (V2)');
colorbar;
axis image;
clim([0 1]);

figure;
imagesc(rot90(flipud(V3_map(:,:,slice_num))));
title('Volume Fraction (V3)');
colorbar;
axis image;
clim([0 1]);

figure;
imagesc(rot90(flipud(V4_map(:,:,slice_num))));
title('Volume Fraction (V4)');
colorbar;
axis image;
clim([0 1]);

% S0 
figure;
imagesc(rot90(flipud(S0_map(:,:,slice_num))));
title('Estimated S0');
colorbar;
axis image;

% T2 maps
figure;
imagesc(rot90(flipud(T2_1_map(:,:,slice_num))));
title('Estimated T2_1 (ms)');
colorbar;
axis image;
clim([0 100]);

figure;
imagesc(rot90(flipud(T2_2_map(:,:,slice_num))));
title('Estimated T2_2 (ms)');
colorbar;
axis image;
clim([0 200]);

figure;
imagesc(rot90(flipud(T2_3_map(:,:,slice_num))));
title('Estimated T2_3 (ms)');
colorbar;
axis image;
clim([0 1500]);

figure;
imagesc(rot90(flipud(T2_4_map(:,:,slice_num))));
title('Estimated T2_4 (ms)');
colorbar;
axis image;
clim([0 3000]);  

%% 10 compartments NNLS
[S0_map, V_maps, mean_residual, RSS_NNLS_10comp] = ...
    estimateT2_NNLS_10comp(images, TEs, brain_mask.img);

% save results
save('T2_NNLS_10comp_results.mat', ...
     'S0_map', 'V_maps', 'mean_residual', 'RSS_NNLS_10comp');

% Visualize results
slice_num = round(size(images, 3) / 2); 

% T2 labels
T2_labels = [10, 20, 30, 50, 80, 120, 300, 500, 1000, 2000];

% Volume fractions
for idx = 1:10
    figure;
    imagesc(rot90(flipud(V_maps(:,:,slice_num, idx))));
    title(['Volume Fraction (V' num2str(idx) ') - T2 = ' num2str(T2_labels(idx)) ' ms']);
    colorbar;
    axis image;
    clim([0 1]);
end

% S0
figure;
imagesc(rot90(fliplr(S0_map(:,:,slice_num))));
title('Estimated S0 (NNLS 10-comp)');
colorbar;
axis image;
