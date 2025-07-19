%% Q6
% Use segmentations to create priors to speed up the fitting process

%% read in the images
r=load_nii('case01-qt2_reg.nii');
r.img(r.img(:)<0)=0;
images = r.img;
brain_mask = load_nii('case01-mask.nii');
TEs=load('case01-TEs.txt');

% Load the segmentation file
segmentation = load_nii('case01-seg.nii');
seg = segmentation.img;

%% NLLS with priors

tic
[S0_map, V1_map, V2_map, V3_map, T2_1_map, T2_2_map, T2_3_map, mean_residual, RSS_NLLS_3comp_priors] = ...
    priors_threeCompartments(images, TEs, brain_mask.img, seg);
time_NLLS_3comp_priors = toc
mean_residual

save('T2_NLLS_3comp_priors_results.mat', ...
     'S0_map', 'V1_map', 'V2_map', 'V3_map', 'T2_1_map', 'T2_2_map', 'T2_3_map', 'mean_residual', 'RSS_NLLS_3comp_priors', 'time_NLLS_3comp_priors');


