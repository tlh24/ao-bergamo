function cmaes_scanimage(src,evt,varargin)

global cm; 

lastStripe = src.hSI.hDisplay.stripeDataBuffer{src.hSI.hDisplay.stripeDataBufferPointer}; % get the pointer to the last acquired stripeData
channels = lastStripe.roiData{1}.channels; % get the channel numbers stored in stripeData

for idx = 1:length(channels)  
	frame{idx} = lastStripe.roiData{1}.imageData{idx}{1}; % extract all channels
	if idx == 1  
		frm = single(frame{idx}); 
		vd = sum(sum(frm));   
		% vd = sum(sum(frm .* cm.mask')); 
	end   
end

k = cm.k; 
md = mod(k, cm.averageFrames+1); 
popn = floor(k / (cm.averageFrames+1)); 
popn = mod(popn, cm.lambda+1); 
if md == 0
	% skip the first frame, it's a transient. 
	cm.vd = 0; 
end
if md > 0 && md < cm.averageFrames
	cm.vdarr(md) = vd; 
end
if md == cm.averageFrames
	cm.vdarr(md) = vd; 
	vds = mean(cm.vdarr(1:cm.averageFrames)); 
	if popn > 0 && popn <= cm.lambda
		cm.costs(popn) = -1.0*vds; 
		disp(['popn ' num2str(popn) ' cost ' num2str(vds)]); 
	end
	if popn == cm.lambda
		cm.samples = cmaes_run(cm.costs);
		disp('generating new samples'); 
	end
	if popn < cm.lambda
		% update the DM.
		cm.udps(single([2.7182818; cm.samples(popn+1, :)'] ));
	end
end
cm.k = k + 1; 

% note: 
% if everything is behaving well (disconneted from internet, e.g), then
% running SI at 73hz / 300 x 512 resolution seems to mostly reliably give
% the right intensity reading after two clocks.  
% this means that in the test above, you'll still get on average 4 clocks
% of reliable data out of the 5 in the 50% duty cycle.  
