load ../centroids.mat
x = double(x); 
y = double(y); 
v = double(v); 
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

% ignore all data that doesn't move.. 
% .. or moves too much. 
mask = (sx > median(sx)/4) .* (st < 2); 
mask = mask>0;
x = x(mask, :); 
y = y(mask, :); 
mx = mean(x, 2); 
my = mean(y, 2); 
dx = x - mx; 
dy = y - my;

nc = size(x, 1); 
nt = size(x, 2); 
A = [dx; dy; ones(1,nt)]; 
% ignore bad frames too
amx = max(abs(A(1:end-1, :)), [], 1); 
tmask = amx < 10; 
A = A(:, tmask); 
v = v(:, tmask); 
amx = amx(tmask); 
figure; 
plot(amx)
title('max absolute dx and dy for all valid centroids');
xlabel('time'); 
amx = std(A(1:end-1, :), [], 2); 
figure; 
plot(amx)
title('std of dx and dy for all time');
xlabel('valid centroid (repeats once)'); 

C = A'\v';  
pred = A' * C; 

err = v' - pred; 

% SVD yields the same coefficient matrix. 
% [U,S,V] = svd(A', 0); 
% C2 = V * (S^-1) * U' * v'; 
% 
% pred2 = A' * C2; 
% err2 = v' - pred2; 

figure; 
subplot(1,2,1); 
imagesc(v'); 
title('DM command signal')
colorbar; 
subplot(1,2,2); 
imagesc(err); 
title('DM command pred error')
colorbar;

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
scatter(mx(indx), -my(indx), abs(q*1e11), colors, 'filled'); 
title('mean weight of C (fwd trans mtx) * 1e11 per centroid'); 

% cx = mean(mx); 
% cy = mean(my); 
% rad = sqrt( (mx - cx).^2 + (my - cy).^2 ); 
% plot(sort(rad))
% radius of 925 should do it. 

% still confused.  
% Try looking at the std of data on each centroid? 
subplot(1,2,2)
Cstd = std(dx, [], 2) + std(dy, [], 2); 
[q, indx] = sort(Cstd); 
scatter(mx(indx), -my(indx), abs(q * 500), colors, 'filled'); 
title('std(dx) + std(dy) for all time per centroid'); 
Cforward = C; 
save('../calibration_forward.mat', 'Cforward', 'mask');

load('../calibration_flat.mat');
