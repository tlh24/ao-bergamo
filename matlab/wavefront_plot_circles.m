function [] = wavefront_plot_circles(wfmx, wfmy, wfdx, wfdy, titl)
% plot the wavefront dx and dy relative to the centroids 
% wfmx wfmy, scaled from red to blue.  
figure;
nc = numel(wfmx); 

% scale both plots by the same min/max so they are comparable.
minx = min(min(wfdx), min(wfdy)); 
maxx = max(max(wfdx), max(wfdy)); 
wfdx = (wfdx - minx)./(maxx - minx); 
wfdy = (wfdy - minx)./(maxx - minx); 

subplot(1,2,1);
red = wfdx; 
blue = 1.0 - wfdx; 
scatter(wfmx, wfmy, 750*ones(nc,1), [red zeros(nc,1) blue],...
	'filled');
axis([200 1700 200 1700]); 
axis square; 
set(gca,'Color','k')
title(['dx for ' titl])

subplot(1,2,2);
red = wfdy; 
blue = 1.0 - wfdy; 
scatter(wfmx, wfmy, 750*ones(nc,1), [red zeros(nc,1) blue],...
	'filled');
axis([200 1700 200 1700]); 
axis square; 
set(gca,'Color','k')
title(['dy for ' titl]); 

set(gcf, 'Position',  [100, 100, 1500, 750])