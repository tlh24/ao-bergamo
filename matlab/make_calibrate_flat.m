load centroids_flat_50hz.mat
x = double(x); 
y = double(y); 
xm = mean(x, 2); 
ym = mean(y, 2); 
mask = xm > 0; 
xm = xm(mask); 
ym = ym(mask);
mask = ym > 0; 
xm = xm(mask); 
ym = ym(mask);

nlenslets = size(xm, 1); 
xcalibctr = mean(xm);
ycalibctr = mean(ym);

figure
scatter(xm, ym, 'b'); 
hold on
plot(xcalibctr, ycalibctr, 'ko'); 

aperture = 900; 

dist = (xm.'-xm).^2 + (ym.'-ym).^2;
imagesc(dist > 1); 
figure
scatter(xm, ym); 
calib = [xm, ym]; 
save('calibration_flat.mat', 'calib'); 