function plot_shwfs_wavefront()
% realy need to plot this relative to a caibrated flat. 
load('../data/calibration_forward.mat', ...
	'cmask', 'Cforward', 'mx', 'my'); 

% get the data. 
mmf = memmapfile('../shared_centroids.dat',...
	'Format','single','Offset',0,'Repeat',6000);

dat = mmf.Data; 
datx = dat(1:2:end); 
daty = dat(2:2:end); 
dx = datx(cmask); 
dy = daty(cmask); 

% need some extra cleanup. 
nc = sum(cmask); 
good = sum(abs(Cforward(1:nc,:)), 2) > 0.00001;

mx = mx(good); 
my = my(good); 
dx = dx(good); 
dy = dy(good); 

load('../data/calibration_960nm_PSbeads_long4_geneopt.mat', ...
	'genecalib'); 
calx = genecalib(good, 1); 
caly = genecalib(good, 2); 

wavefront_plot_circles(mx, my, dx-calx, dy-caly, 'SHWFS'); 