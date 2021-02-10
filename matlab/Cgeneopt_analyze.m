load('../data/calibration_forward.mat'); 
load('../data/calibration_flat.mat'); 
load('../data/DMoptimization_run5.mat'); 

% this takes the output of a genetic optimization run, 
% runs some diagnostics, 
% and outputs the 'best' DM control signal (not used usually)
% and the corresponding WFS reading 
% the latter is assumed to be 'flat'; 
% modifications are made from it

figure; 
frames_sum = sum(save_sumstd(:, 2:10), 2); 
plot(frames_sum); 

nlenslets = sum(mask); 
[~, indx] = sort(frames_sum(9000:end), 'descend'); 
indx = indx(1:100) + 9000 - 1; 
avg_wfs_x = mean(double(save_wfs_dx(indx, :)), 1)'; 
avg_wfs_y = mean(double(save_wfs_dy(indx, :)), 1)'; 
avg_wfs_x = avg_wfs_x(cmask(1:nlenslets)); 
avg_wfs_y = avg_wfs_y(cmask(1:nlenslets)); 

figure;
scatter3(mx, my, avg_wfs_x - mx); 
title('average dx for top solutions'); 

figure
scatter3(mx, my, avg_wfs_y - my); 
title('average dy for top solutions'); 

figure; 
avg_dmcommand = mean(save_dmcommand(indx, :), 1); 
colors = jet(200); 
colorindx = round((avg_dmcommand + 0.15) / 0.3 * 200); 
colorindx = min(colorindx, 200); 
colorindx = max(colorindx, 1); 
[dmx, dmy] = dm_actuator_to_xy();
scatter(dmx, dmy, 250, colors(colorindx, :), 'filled')
title('average DM control signals for top solutions'); 

Best_DMcommand = avg_dmcommand; 
genecalib = [avg_wfs_x avg_wfs_y];
save('../data/calibration_geneopt.mat','Best_DMcommand','genecalib'); 
% save the absolute centroid positions
% convert to relative later, depends on forward-calibration.

movie_frames = single(zeros(256, 256, 1, 10000/4)); 
maxx = max(save_frames, [], 'all'); 
maxx = maxx / 4.0; 
minn = min(save_frames, [], 'all');  
for k = 1:4:10000
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
implay(mov);