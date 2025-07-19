function VaryingEchoTime(img, TE, slice_num)
% VaryingEchoTime Visualizes the effect of varying echo times in a T2-weighted MRI scan.
%
%   VaryingEchoTime(img, TE, slice_num) displays a series of grayscale images
%   from a selected slice (slice_num) of a 4D MRI dataset (img) acquired at 
%   different echo times (TE). 
%
%   The function performs two visualizations:
%     1. A grid of static subplots showing how signal intensity decays with 
%        increasing TE for the selected slice, using radiological orientation.
%
%     2. A video animation ('T2_animation.mp4') that cycles through the same
%        slice across all TE values, illustrating the temporal evolution of
%        signal decay.
%
%   Inputs:
%     - img: 4D matrix of MRI images [x, y, z, t]
%     - TE:  Vector containing echo time values (in milliseconds)
%     - slice_num: Index of the slice to visualize

%% One subplot for each echo-time
figure;
colormap gray;

for t = 1:length(TE)
    subplot(4, ceil(length(TE)/4), t); 
    imagesc(rot90(flipud(img(:,:,slice_num,t)), 1));
    axis off;
    title(['TE = ', num2str(TE(t)), ' ms']);
end


%% Video animation

% Create video objecct
v = VideoWriter('T2_animation.mp4', 'MPEG-4');
v.FrameRate = 2; 
open(v);

figure;
colormap gray;

% Loop over all echo times to generate and show each frame
for t = 1:length(TE)
   
    % Flip left-right and rotate 90° counterclockwise for radiological view
    imagesc(rot90(fliplr(img(:,:,slice_num,t))));
    
    axis off;
    title(['TE = ', num2str(TE(t)), ' ms']);
    drawnow;

    frame = getframe(gcf);   % Capture the current frame
    writeVideo(v, frame);    % Write the frame to the video
end

close(v);  % Finalize and close the video file