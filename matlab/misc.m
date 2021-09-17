figure
hold on
for i=1:12
	red = (i-1)/11.0; 
	blue = 1.0-red; 
	plot(save_sumstd(1:k, i), 'color', [red, 0.0, blue], ...
		'LineWidth', 2);
end
% figure;
% plot(std(save_dmcommand(1:k, :)'))
figure; 
plot(diff(save_sumstd(1:k, :)'))

figure; 
save_time2 = reshape(save_time(1:k, :)', k*12, 1);
hist(diff(save_time2(1:end-20)), 100)
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
