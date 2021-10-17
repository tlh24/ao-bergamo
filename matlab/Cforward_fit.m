load ../rundata/centroids_full.mat
% x = double(x); 
% y = double(y); 
% v = double(v); 
% n = min(384000, size(x,2)); % because my laptop has limited memory. 
% x = x(:, 1:n); 
% y = y(:, 1:n); 
% v = v(:, 1:n);
dx = diff(x, 1, 2); 
dy = diff(y, 1, 2); 

sx = sum(abs(dx) + abs(dy), 2);
figure; 
plot(sx); 
title('summed dx/dt and dy/dt'); 
xlabel('centroid no'); 

mx = mean(x, 2); 
my = mean(y, 2); 
dx = x - mx; 
dy = y - my;
st = std(dx, [], 2) + std(dy, [], 2); 

% ignore all centroids that don't move
% .. or move too much. 
cmaskr = (sx > median(sx(sx > 0))/8) .* (st < 10); 
cmaskr = cmaskr>0;
hold on
plot(cmaskr * 1e7, 'r')
load('../data/calibration_forward.mat', 'cmask'); 
cmask = cmask(1:1100); 
plot(cmask * 1e7, 'k')
legend('summed dx/dt and dy/dt', 'new cmask', 'old cmask'); 
disp(['number of active centroids: ' num2str(sum(cmaskr))]);
answer = input('Update cmask? (1 / 0): '); 
if answer == 1
	cmask = cmaskr; 
else
	% load in cmask from file -- to keep the other files 
	% (geneopts) valid.
	% need to zero all noisy data anyway, so it doesn't slip into cforward. 
	cmaskrr = repmat(cmaskr, 1, size(x, 2)); 
	x = x .* cmaskrr; 
	y = y .* cmaskrr; 
	clear cmaskrr; 
end
x = x(cmask, :); 
y = y(cmask, :); 
mx = mean(x, 2); 
my = mean(y, 2); 
dx = x - mx; 
dy = y - my;

nc = size(x, 1); 
nt = size(x, 2); 
A = [dx; dy; ones(1,nt)]; 
% ignore bad frames too
amx = max(abs(A(1:end-1, :)), [], 1); 
if 0 
	tmask = (amx < 20) .* (rand(1,nt) > 0.5); 
	tmask = tmask > 0; 
else 
	tmask = (amx < 20); 
end
A = A(:, tmask); 
v = v(:, tmask); 

amx = amx(tmask); 
% figure; 
% plot(amx)
% title('max absolute dx and dy for all valid centroids');
% xlabel('time'); 
amx = std(A(1:end-1, :), [], 2); 
figure; 
plot(amx)
title('std of dx and dy for all time');
xlabel('valid centroid (repeats once)'); 

C = A'\v';  
pred = A' * C;  

err = v' - pred; 

if 1
	% SVD yields the same coefficient matrix. 
	[U,S,V] = svd(A', 0); 
	C2 = V * (S^-1) * U' * v'; 

	pred2 = A' * C2; 
	err2 = v' - pred2; 
end

figure; 
subplot(1,3,1); 
imagesc(v', [-0.12 0.12]); 
title('DM command signal')
colorbar; 
subplot(1,3,2); 
imagesc(err, [-0.012 0.012]); 
title('DM command pred error')
colorbar;
subplot(1,3,3); 
imagesc(pred, [-0.12 0.12]); 
title('DM pred')
colorbar;

figure; 
plot(std(err, 1), 'r','Linewidth', 2)
hold on
plot(std(pred, [], 1), 'b','Linewidth', 4)
plot(std(v', [], 1), 'g','Linewidth', 2); 
title('red = std. of pred error, blue = std. of prediction, black=data'); 
xlabel('DM actuator number'); 

figure; 
imagesc(atan(zscore(C, 0, 'all'))); 
title('weight matrix, zscored & atan transformed');

% uh, it looks like some of the actuators are not observed 
% by the wavefront sensor. 
Cp = C(1:nc, :) + C(nc+1:nc*2, :); 
colors = jet(nc); 
if 0
    figure;
    for i = 1:97
        [q, indx] = sort(Cp(:, i)); 
        subplot(1,2,1)
        scatter(mx(indx), my(indx), abs(q)*1200, colors, 'filled'); 
        subplot(1,2,2)
        plot(q); 
        pause; 
    end
end

% alright, confused, let's also look at the mean weight in C. 
Cmn = mean(Cp, 2); 
figure; 
subplot(1,2,1)
[q, indx] = sort(Cmn); 
scatter(mx(indx), -my(indx), ceil(abs(q*1e11)+0.001), colors, 'filled'); 
title('mean weight of C (fwd trans mtx) * 1e11 per centroid'); 

% Check by looking at the std of data on each centroid
subplot(1,2,2)
Cstd = std(dx, [], 2) + std(dy, [], 2); 
[q, indx] = sort(Cstd); 
scatter(mx(indx), -my(indx), abs(q * 500 + 0.1), colors, 'filled'); 
title('std(dx) + std(dy) for all time per centroid'); 
Cforward = C; 
cmask2 = zeros(3000, 1); 
cmask2(1:1100) = cmask; 
cmask = cmask2>0; 

VS = V*S; % multiply by the SVD values.
VS = VS / 500; % scale VS by +- 1.0 to drive the synthetic wavefront. 
VS = VS(:, 1:97); % we don't care about higher mechanical modes (no degrees of freedom in the system)
VS(nc*2+1, :) = 0.0; 

answer = input('Save calibration_forward.mat? (1 / 0): '); 
if answer == 1
	save('../data/calibration_forward.mat', 'Cforward', 'cmask', 'mx', 'my', 'cmaskr','VS');
end

% also need to save for pytorch model fitting: 
% note hard-coded flat calibration!! (Needs to be a hdf5 file, might as
% well splice in here)
load('../data/calibration_960nm_PSbeads_long4_geneopt.mat', 'genecalib', 'Best_DMcommand')
flatwf = ones(nc*2+1, 1); % last value is bias / 1
flatwf(1:nc) = genecalib(:,1) - mx; % let's hope the microscope optics haven't moved between these...
flatwf(nc+1:nc*2) = genecalib(:,2) - my; 
flatbad = abs(flatwf) > 10;
flatwf(flatbad) = 0.0; 
VS = VS'; % matlab saves in Fortran order... numpy reads in C order.
save('../rundata/centroids_full_cleaned.mat','-v7.3',...
	'A','v','cmask','VS','flatwf')


answer = input('plot the SVD modes? (1 / 0)');
if answer == 1
    for j = 40:-1:1
        wavefront_plot_circles(mx, my, V(1:nc, j), V(nc+1:2*nc, j), ['svd wavefront ' num2str(j)])
    end
end
