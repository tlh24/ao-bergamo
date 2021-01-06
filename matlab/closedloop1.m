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
% ignore all data that doesn't move.. 
mask = sx > mean(sx)/2; 
x = x(mask, :); 
y = y(mask, :); 

mx = mean(x, 2); 
my = mean(y, 2); 
dx = x - mx; 
dy = y - my; 

nt = size(x, 2); 
nc = size(x, 1); 
A = [dx; dy; ones(1,nt)]; 
% ignore blips here too
amx = max(abs(A(1:end-1, :)), [], 1); 
figure; 
plot(amx)
title('max absolute dx and dy for all centroids');
xlabel('time'); 
amx = std(A(1:end-1, :), [], 2); 
figure; 
plot(amx)
title('std of dx and dy for all time');
xlabel('centroid'); 

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
figure;
colors = jet(nc); 
Cp = C(1:nc, :) + C(nc+1:nc*2, :); 
for i = 1:97
    [q, indx] = sort(Cp(:, i)); 
    subplot(1,2,1)
    scatter(mx(indx), my(indx), abs(q)*1200, colors, 'filled'); 
    subplot(1,2,2)
    plot(q); 
    pause; 
end

% alright, confused, let's also look at the mean weight in C. 
Cmn = mean(Cp, 2); 
figure; 
subplot(1,2,1)
[q, indx] = sort(Cmn); 
scatter(mx(indx), my(indx), abs(q * 5000), colors, 'filled'); 

cx = mean(mx); 
cy = mean(my); 
rad = sqrt( (mx - cx).^2 + (my - cy).^2 ); 
plot(sort(rad))
% radius of 925 should do it. 

% still confused.  
% Try looking at the std of data on each centroid? 
subplot(1,2,2)
Cstd = std(dx, [], 2) + std(dy, [], 2); 
[q, indx] = sort(Cstd); 
scatter(mx(indx), my(indx), abs(q * 500), colors, 'filled'); 

% I would say that if the mean weight is off, 
% we should exclude the centroid from the control algo. 
goodmask = (Cmn < (mean(Cmn) + std(Cmn))) .* (Cmn > (mean(Cmn) - std(Cmn)));

% aand ... do it again. 
x = x(goodmask>0, :); 
y = y(goodmask>0, :); 

mx = mean(x, 2); 
my = mean(y, 2); 
dx = x - mx; 
dy = y - my; 

nt = size(x, 2); 
nc = size(x, 1); 
A = [dx; dy; ones(1,nt)]; 
C = A'\v';  
pred = A' * C; 

err = v' - pred; 

figure; 
subplot(1,2,1); 
imagesc(v'); 
title('DM command signal')
colorbar; 
subplot(1,2,2); 
imagesc(err); 
title('DM command pred error')
colorbar;

colors = jet(nc); 
Cp = C(1:nc, :) + C(nc+1:nc*2, :); 
Cmn = mean(Cp, 2); 
figure; 
[q, indx] = sort(Cmn); 
scatter(mx(indx), my(indx), abs(q * 5000), colors, 'filled'); 