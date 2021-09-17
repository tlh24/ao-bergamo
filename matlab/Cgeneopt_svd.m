% run Cforward_fit first, to populate/clean the variables. 
% namely: 'A' (wavefront dx dy matrix) and 'v' (control signals)
% if 0 
% 	[U, S, V] = svd(A', 'econ');
% 
% 	plot(diag(S))
% 
% 	% A' = U*(S*V'); 
% 
% 	for i = 1:10
% 		vx = V(1:nc, i); 
% 		vy = V(nc+1:2*nc, i); 
% 		figure;
% 		subplot(1,2,1);
% 		scatter3(mx, my, vx);
% 		title(['dx for svd component ' num2str(i)])
% 		subplot(1,2,2);
% 		scatter3(mx, my, vy);
% 		title(['dy for svd component ' num2str(i)])
% 	end
% 	% unsurprisingly, the first order is tip / tilt.. 
% 	% both! (SVD is agnostic up to a rotation..)
% 
% 	% but these orthogonal axes only matter with respect to 
% 	% image metrics -- look at that now.
% end

load('../data/calibration_forward.mat'); % forward, mx, my
% load('../data/calibration_flat.mat'); % calib
fname = '1080nm_orangePSbeads_long1';
load(['../rundata/DMoptimization_' fname '.mat']); 
nc = numel(mx); 



% 
% figure; 
% S2 = S(1:10, :);
% V2 = S2 * V'; 
% traj = wf * V2';
% trajsmooth = sgolayfilt(double(traj), 3, 41, [], 1); 
% plot(trajsmooth(:, 1:10)); 
% legend
% title('trajectory of SVD components 1-10 of wf sensor under GA opt');

% note: this works because V is an orthonormal matrix
% and so S * V' is orthogonal
% hence their matrix inverse is ~= the transpose!
% this seems to suggest that we really only need ~4 components
% to optimize the image quality, at least in these synthetic samples. 
% Need to take repeated measurements to make sure this is true!
% (This is assuming that the genetic optimizer is unbiased)

% now, need to output a 'best wavefront' matrix, 
frames_sum = sum(save_sumstd(:, 7:12),2); 
plot(frames_sum)
[x, indx] = sort(frames_sum, 'desc'); 
hold on
plot(indx(1:400), x(1:400), 'ro'); 
best_wfs_x = median(save_wfs_dx(indx(1:400), :), 1);
best_wfs_y = median(save_wfs_dy(indx(1:400), :), 1);

% save these in the file.. 
Best_DMcommand = double(median(save_dmcommand(indx(1:400), :), 1));
genecalib = double([best_wfs_x' best_wfs_y']); 

% as well as these first 10 SVD components for online correction.
wf = [(save_wfs_dx-best_wfs_x) (save_wfs_dy-best_wfs_y) ...
	ones(size(save_wfs_dx,1), 1)]; 
wf = wf(1:end-1, :);
frames_sum = sum(save_sumstd(:, 7:12), 2); 
coef = wf \ frames_sum(1:end-1); 
pred = wf * coef; 
figure; 
plot(frames_sum(1:end-1));
hold on
plot(pred)

% % try SVD implicitly weighted based on GA optim. 
[U, S, V] = svd(wf, 'econ'); 
for i = 1:20
	vx = V(1:nc, i); 
	vy = V(nc+1:2*nc, i); 
	figure;
	subplot(1,2,1);
	red = max(min((vx + 0.2)/0.4, 1.0), 0.0); 
	blue = max(min((0.2 - vx)/0.4, 1.0), 0.0); 
	scatter(mx, my, 750*ones(nc,1), [red zeros(nc,1) blue], 'filled');
	axis([200 1700 200 1700]); 
	axis square; 
	set(gca,'Color','k')
	title(['dx for GA svd component ' num2str(i)])
	subplot(1,2,2);
	red = max(min((vy + 0.2)/0.4, 1.0), 0.0); 
	blue = max(min((0.2 - vy)/0.4, 1.0), 0.0); 
	scatter(mx, my, 750*ones(nc,1), [red zeros(nc,1) blue], 'filled');
	axis([200 1700 200 1700]); 
	axis square; 
	set(gca,'Color','k')
	title(['dy for GA svd component ' num2str(i)])
	set(gcf, 'Position',  [100, 100, 1500, 750])
end

% do SVD implicitly weighted based on GA optim & 
% difference from the 'best' corrected wavefront. 
% so, when we do online adjustment, we tweak the wavefront by 
% these vectors; setting everything to zero = 'best'. 
[U, S, V] = svd(double(wf(:,1:end-1)), 'econ'); 

svd_wfs_x = V(1:nc, 1:10); 
svd_wfs_y = V(nc+1:nc*2, 1:10); 

% need to get the permissible ranges for scaling these wavefronts. 
K = U*S; 
meanK = mean(K, 1); 
stdK = std(K, 1);
figure; 
subplot(1,2,1);
plot(meanK(1:20));
subplot(1,2,2);
plot(stdK(1:20));

% allow scaling by +-1.5 std, except for the first component, which = 2nd
svd_scl = stdK(1:10) * 8.0; 
svd_scl(1) = svd_scl(2); 
% --------
save(['../data/calibration_' fname '_geneopt.mat'], ...
'Best_DMcommand','genecalib','svd_wfs_x','svd_wfs_y','svd_scl'); 

scatter3(best_wfs_x, best_wfs_y, coef(1:nc)/1e5, 'o'); 