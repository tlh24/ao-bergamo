load ../centroids.mat
x = double(x); 
y = double(y); 
v = double(v); 
dx = diff(x, 1, 2); 
dy = diff(y, 1, 2); 

sx = sum(abs(dx) + abs(dy), 2);
plot(sx); 
% ignore all data that doesn't move.. 
mask = sx > 1000; 

x = x(mask, :); 
y = y(mask, :); 

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

subplot(1,2,1); 
imagesc(v'); 
colorbar; 
subplot(1,2,2); 
imagesc(err); 
colorbar;

% uh, it looks like some of the actuators are not observed 
% by the wavefront sensor. 
figure;
colors = jet(97); 
Cp = abs(C(1:nc, :)) + abs(C(nc+1:nc*2, :)); 
for i = 1:97
    [q, indx] = sort(abs(Cp(:, i))); 
    indx = indx(end-1:end-30); 
    scatter(mx(indx), my(indx), 30, colors(i, :)); 
    hold on; 
end