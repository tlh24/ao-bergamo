function cmaes_si_reset(src,evt,varargin)

global cm; 

cm.udps = dsp.UDPSender('RemoteIPPort', 13131, ...
    'RemoteIPAddress','10.123.1.152'); 
 
cmaes_init; 
% run cmaes once to get an initial population.
cm.costs = zeros(cm.lambda, 1); 
cm.samples = cmaes_run(cm.costs);

cm.mask = zeros(300, 512); 
for r = 1:300
	for c = 1:512
		if abs(r-150) < 60 && abs(c-256) < 80
			cm.mask(r,c) = 1; 
		end
	end
end

cm.k = 0; 
cm.averageFrames = 4; 