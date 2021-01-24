load('../centroids.mat')
x = double(x); 
y = double(y); 
xm = mean(x, 2); 
ym = mean(y, 2); 
v = std(x, 0, 2) + std(y, 0, 2);
mask = (xm > 0) .* (v < 1) .* (ym > 0); 
mask = mask > 0; 
xm = xm(mask); 
ym = ym(mask);

nlenslets = size(xm, 1) 
xcalibctr = mean(xm);
ycalibctr = mean(ym);

figure
scatter(xm, ym, 'b'); 
hold on
plot(xcalibctr, ycalibctr, 'ko'); 

figure
dist = (xm.'-xm).^2 + (ym.'-ym).^2;
imagesc(dist > 10); 
title('distance matrix between centroids; should be no off-axis'); 

figure
plot(std(x(mask, :), 0, 2) + std(y(mask, :), 0, 2));
title('std of centroid positions'); 

calib = [xm, ym]; 
save('../calibration_flat.mat', 'calib'); 