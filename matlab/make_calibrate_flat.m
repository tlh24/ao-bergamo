load('../data/calibrated_flat_14mm300um_Feb4.mat')
x = double(x); 
y = double(y); 
xm = mean(x, 2); 
ym = mean(y, 2); 
v = std(x, 0, 2) + std(y, 0, 2);
mask = (xm > 0) .* (v < 1) .* (ym > 0); 
mask = mask > 0; 
figure
scatter(xm, ym, 'r'); 
xm = xm(mask); 
ym = ym(mask);
hold on
scatter(xm, ym, 'b'); 

nlenslets = size(xm, 1) 
xcalibctr = mean(xm);
ycalibctr = mean(ym);
plot(xcalibctr, ycalibctr, 'ko'); 

figure
dist = (xm.'-xm).^2 + (ym.'-ym).^2;
imagesc(dist > 10); 
title('distance matrix between centroids; should be no off-axis'); 

figure
plot(std(x(mask, :), 0, 2) + std(y(mask, :), 0, 2));
title('std of centroid positions'); 

calib = [xm, ym]; 
save('../data/calibration_flat.mat', 'calib'); 
clear mask % because it may confuse downstream scripts -- 
% shwfs reads in the centroids coalesced. 