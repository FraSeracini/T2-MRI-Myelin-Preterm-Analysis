%% Q1
% We identify and visualize how the images change with varying echo-times
% and we plot some time-intensity curves for some single voxels, ROIs and
% the full brain to investigate if the data follow a monotonic and
% mon-exponential decay

%% Load the data
% read in the images
r=load_nii('case01-qt2_reg.nii');
r.img(r.img(:)<0)=0;
images = r.img;

TEs=load('case01-TEs.txt');
% disp(TEs);

% Select a slice
slice_num = round(size(images,3)/2); % A central slice


%% Visualize how the imaging data change with different echo-times
VaryingEchoTime(images, TEs, slice_num)

%% Visualize the Time-Intensity curve for a single voxel

% Choose a voxel
x = 45;
y = 22;
z = 45;

% Extract the signal for all the TEs
signal_voxel = squeeze(images(x, y, z, :)); 

% Plot the time intensity curve for a single voxel
figure;
plot(TEs, signal_voxel, 'x-', 'LineWidth', 0.5);
xlabel('Echo Time (ms)');
ylabel('Signal Intensity');
title('Time-Intensity Curve for one voxel');
grid on;

%% Visualize the Time-Intensity curve for a ROI

% Choose a ROI
x_range = 22:50; 
y_range = 35:55; 
z_slice = 22; 

% Compute mean intensity for each TE in the ROI
roi_signal = squeeze(mean(mean(images(x_range, y_range, z_slice, :), 1), 2));

% Plot the time intensity curve for the selected ROI
figure;
plot(TEs, roi_signal, 'x-', 'LineWidth', 0.5);
xlabel('Echo Time (ms)');
ylabel('Mean Signal Intensity');
title('Time-Intensity Curve for selected ROI');
grid on;

%% Visualize the Time-Intensity curve for the White Matter

% Load the segmentation file
seg = load_nii('case01-seg.nii');

% Select the white matter mask
wm_mask = seg.img(:, :, :, 4);

% Initialize
roi_signal = zeros(size(TEs)); 

for t = 1:length(TEs)
    current_image = squeeze(images(:, :, z_slice, t)); 
    wm_voxels = wm_mask(:, :, z_slice) > 0.7; % Select only the voxels that have a probability > 0.7
    roi_signal(t) = mean(current_image(wm_voxels)); 
end

% Plot the time intensity curve for the white matter
figure;
plot(TEs, roi_signal, 'x-', 'LineWidth', 0.5);
xlabel('Echo Time (ms)');
ylabel('Mean Signal Intensity');
title('Time-Intensity Curve for WM');
grid on;

%% %% Visualize the mean Time-Intensity curve across all the brain

% Mean of the signals for the brain
mean_signal = squeeze(mean(mean(mean(images,1),2),3)); 

% Plot the time intensity curve for the mean brain signal
figure;
plot(TEs, mean_signal, 'x-', 'LineWidth', 0.5);
xlabel('Echo Time (ms)');
ylabel('Mean Signal Intensity');
title('Average Time-Intensity Curve across the brain');
grid on;

