function [] = wavefront_plot_circles(wfmx, wfmy, wfdx, wfdy, titl)
% plot the wavefront dx and dy relative to the centroids 
% wfmx wfmy, scaled from red to blue.  
figure;
nc = numel(wfmx); 

% remove tip, tilt, and piston.
wfdx = wfdx - mean(wfdx); 
wfdy = wfdy - mean(wfdy); 

% scale both plots by the same min/max so they are comparable.
minx = min(min(wfdx), min(wfdy)); 
maxx = max(max(wfdx), max(wfdy)); 
wfdxo = (wfdx - minx)./(maxx - minx); 
wfdyo = (wfdy - minx)./(maxx - minx); 

subplot(1,2,1);
red = wfdxo; 
blue = 1.0 - wfdxo; 
scatter(wfmx, wfmy, 750*ones(nc,1), [red zeros(nc,1) blue],...
	'filled');
axis([200 1700 200 1700]); 
axis square; 
set(gca,'Color','k')
title(['dx for ' titl])

subplot(1,2,2);
red = wfdyo; 
blue = 1.0 - wfdyo; 
scatter(wfmx, wfmy, 750*ones(nc,1), [red zeros(nc,1) blue],...
	'filled');
axis([200 1700 200 1700]); 
axis square; 
set(gca,'Color','k')
title(['dy for ' titl]); 

set(gcf, 'Position',  [100, 100, 1500, 750])

% microlens array has a focal length of 14.6mm
% pixel pitch is 5.5um
% relay has a mag of 2.5x to fit the sensor pupil
text(0, 0, ['min (blue) ' num2str(minx*5.5/14600 / 2.5) ...
	' max (red) ' num2str(maxx*5.5/14600 / 2.5) ' radians'] , ...
	'FontSize', 24)

% really need the quantitative wavefront here! 
cx = mean(wfmx); 
cy = mean(wfmy); 
wfmxo = wfmx - cx; 
wfmyo = wfmy - cy; 
lo = min(min(wfmxo), min(wfmyo)) - 30; 
hi = max(max(wfmxo), max(wfmyo)) + 30; 
wfmxo = ((wfmxo - lo) / (hi - lo) - 0.5)*2; 
wfmyo = ((wfmyo - lo) / (hi - lo) - 0.5)*2; 
[Z,dZx,dZy] = zernikes_and_derivatives_cartesian_OSA(...
	wfmxo, wfmyo, 13,'zero');

% now need to fit dZx and dZy -> r
dZx = squeeze(dZx); 
dZy = squeeze(dZy); 
Z = squeeze(Z); 
A = [dZx; dZy];
B = [wfdx*5.5/14600/2.5; wfdy*5.5/14600/2.5]; % convert to radians.
coef = A\B;
pred = A*coef;
predx = pred(1:533);
predy = pred(534:end);
figure; 
plot(B); 
hold on; 
plot(pred, 'r'); 
% figure; 
% subplot(1,2,1)
% scatter3(wfmx, wfmy, predx); 
% subplot(1,2,2)
% scatter3(wfmx, wfmy, predy); 

% figure; 
% plot(coef); 
% title('coef');
% kill tip / tilt and piston
% coef(1:3) = 0; 

% now just need to get the units right for the wavefront. 
% wfdx and wfdy are in radians (angle)
% so, orders 2 and 3 in Z go +- 2.0 over the unit pupil
% dZx and dZy accordingly are slope 2. 
% to get everything in real units, need to scale by the size of the real
% pupil. 
wfscale = (hi-lo) * 5.5e-6 / 2.0; % ~ 3.8mm <-> 1 normalized
wavefront = Z*coef .* wfscale; 

figure; 
T = delaunay(wfmx*5.5, wfmy*5.5); 
trisurf(T, wfmx*5.5, wfmy*5.5, wavefront, 'EdgeColor', 'none')
shading interp
colorbar
xlabel('centroid X position, um');
ylabel('centroid Y position, um');
view(2)
set(gcf, 'Position',  [100, 100, 700, 650])
title(['waefront for ' titl])
