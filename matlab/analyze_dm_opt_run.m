load('../calibration_forward.mat'); 

% load('DMoptimization_run4.mat'); 

frames_sum = sum(save_frames, 2);
frames_sum = sum(frames_sum, 3);
frames_sum = squeeze(frames_sum);

figure; 
plot(frames_sum); 

[~, indx] = sort(frames_sum(14000:end), 'descend'); 
indx = indx(1:100) + 14000 - 1; 
avg_wfs_x = mean(save_wfs_dx(indx, :), 1)'; 
avg_wfs_y = mean(save_wfs_dy(indx, :), 1)'; 
okay = (abs(avg_wfs_x-mx) + abs(avg_wfs_y-my)) < 3;
avg_wfs_x = avg_wfs_x(okay);
avg_wfs_y = avg_wfs_y(okay);

figure;
scatter3(mx(okay), my(okay), avg_wfs_x - mx); 
title('average dx for top solutions'); 

figure
scatter3(mx(okay), my(okay), avg_wfs_y - my); 
title('average dy for top solutions'); 

figure; 
avg_dmcommand = mean(save_dmcommand(indx, :), 1); 
colors = jet(200); 
colorindx = round((avg_dmcommand + 0.15) / 0.3 * 200); 
colorindx = min(colorindx, 200); 
colorindx = max(colorindx, 1); 
scatter(dmx, dmy, 250, colors(colorindx, :), 'filled')
title('average DM control signals for top solutions'); 

Best_DMcommand = avg_dmcommand; 
save('../Best_DMcommand4.mat', 'Best_DMcommand', ...
	'avg_wfs_x', 'avg_wfs_y', 'okay'); 
% save the absolute centroid positions
% convert to relative later, depends on forward-calibration.

movie_frames = single(zeros(256, 256, 1, 15000/4)); 
maxx = max(save_frames, [], 'all'); 
maxx = maxx / 4.0; 
minn = min(save_frames, [], 'all');  
for k = 1:4:15000
	frame = save_frames(k, :, :); 
	frame = frame + save_frames(k+1, :, :);
	frame = frame + save_frames(k+2, :, :);
	frame = frame + save_frames(k+3, :, :);
	frame = (frame-minn)/(maxx-minn)/4.0 * 254 + 1; 
	frame = max(frame, 1); 
	frame = min(frame, 255); 
	movie_frames(:,:, 1, (k-1)/4+1) = frame; 
end
cmap = repmat(linspace(0,1,255)', 1, 3);
mov = immovie(movie_frames,cmap);