load ../centroids.mat
x = double(x); 
y = double(y); 
v = double(v); 
dx = diff(x, 1, 2); 
dy = diff(y, 1, 2); 

sx = sum(abs(dx) + abs(dy), 2);
plot(sx); 
% ignore all data that doesn't move.. 
mask = sx > 500; 

x = x(mask, :); 
y = y(mask, :); 

mx = mean(x, 2); 
my = mean(y, 2); 
dx = x - mx; 
dy = y - my; 

A = [dx; dy]; 
C = A'\v'; 
pred = A' * C; 

err = v' - pred; 

subplot(1,2,1); 
imagesc(v'); 
colorbar; 
subplot(1,2,2); 
imagesc(err); 
colorbar;

        