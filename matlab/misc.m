figure
hold on
% nsamp = 24
for i=1:nsamp
	red = (i-1)/(nsamp-1); 
	blue = 1.0-red; 
	plot(save_sumstd(1:k, i), 'color', [red, 0.0, blue], ...
		'LineWidth', 2);
end
% figure;
% plot(std(save_dmcommand(1:k, :)'))
figure; 
plot(diff(save_sumstd(1:k, :)'))

figure; 
save_time2 = reshape(save_time(1:k, :)', k*nsamp, 1);
hist(diff(save_time2(1:end-20)), 100)
title('histogram of inter-frame interval. (should be ~2.5ms)'); 
% ----------

load ../data/calibration_forward.mat
load ../rundata/centroids.mat

cmask = cmask(1:1100); 
x = x(cmask, :); 
y = y(cmask, :); 
x = x - mx; 
y = y - my; 
cmaskr = cmaskr(cmask); 
x = x(cmaskr, :); 
y = y(cmaskr, :); 
wf = double([x' y']); 

[U,S,V] = svd(wf, 'econ'); 
nc = size(x, 1); 
mx = mx(cmaskr); 
my = my(cmaskr); 

for i = 1:20
	vx = V(1:nc, i); 
	vy = V(nc+1:2*nc, i); 
	figure;
	subplot(1,2,1);
	red = max(min((vx + 0.15)/0.3, 1.0), 0.0); 
	blue = max(min((0.15 - vx)/0.3, 1.0), 0.0); 
	scatter(mx, my, 750*ones(nc,1), [red zeros(nc,1) blue], 'filled');
	axis([200 1700 200 1700]); 
	axis square; 
	set(gca,'Color','k')
	title(['dx for GA svd component ' num2str(i)])
	subplot(1,2,2);
	red = max(min((vy + 0.15)/0.3, 1.0), 0.0); 
	blue = max(min((0.15 - vy)/0.3, 1.0), 0.0); 
	scatter(mx, my, 750*ones(nc,1), [red zeros(nc,1) blue], 'filled');
	axis([200 1700 200 1700]); 
	axis square; 
	set(gca,'Color','k')
	title(['dy for GA svd component ' num2str(i)])
	set(gcf, 'Position',  [100, 100, 1500, 750])
end

% -------
% for looking at SVD components from random perturbation..
% run after Cforward_fit. 
for i = 20:-1:1
	vx = V(1:nc, i); 
	vy = V(nc+1:2*nc, i); 
	wavefront_plot_circles(mx, my, vx, vy, ...
		['Random perturb svd component ' num2str(i)]); 
	set(gcf, 'Position',  [100, 100, 1500, 750])
end
% useful!  seems like a good basis for optimization. 
% need to project this into actuator space, though: 
% we don't have time for closed-loop control in the online case. 
% the columns of V correspond to singular vectors, so transpose.
Vdm = V' * Cforward; 
[dmx,dmy] = dm_actuator_to_xy(); 
for i = 20:-1:1
	va = Vdm(i, :)'; 
	figure; 
	red = clamp((200*va+0.5), 0, 1); 
	blue = clamp((0.5 - 200*va), 0, 1); 
	scatter(dmx, dmy, 2000*ones(97,1), ...
		[red zeros(97,1) blue], 'filled'); 
	axis([-6 6 -6 6]); 
	axis square; 
	set(gca,'Color','k')
	title(['DM actuator for wf SVD component ' num2str(i)]); 
end

% I think, because the SVD bases are orthogonal, it *should*
% be OK to do line optimization in this dimension (e.g. the wavefront
% dim)... so, we do iterative optimization along each dimension, in a
% series, with the un-modified GA loop, just with variance on only on dim
% at a time. 
% Skip the defocus dim, of course. 
% tip-tilt is automaticall removed, and does not occurr in the bases. 
