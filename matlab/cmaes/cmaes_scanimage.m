function cmaes_scanimage(src,evt,varargin)

global genop; 

lastStripe = src.hSI.hDisplay.stripeDataBuffer{src.hSI.hDisplay.stripeDataBufferPointer}; % get the pointer to the last acquired stripeData
channels = lastStripe.roiData{1}.channels; % get the channel numbers stored in stripeData

for idx = 1:length(channels)  
    frame{idx} = lastStripe.roiData{1}.imageData{idx}{1}; % extract all channels
    if idx == 1  
        frm = single(frame{idx}); 
% 		mn = mean(mean(frm)); 
% 		st = std(std(frm)); 
% 		genop.avstd = 0.99 * genop.avstd + 0.01 * st; 
% 		msk = frm > mn + 0.5*genop.avstd; 
% 		vd = sum(sum(frm));   
		vd = sum(sum(frm .* genop.mask')); 
    end   
end

k = genop.k; 
N = genop.N; 
md = mod(k, genop.averageFrames+1); 
if md == 0
    % skip the first frame, it's a transient. 
	genop.vd = 0; 
    genop.k = k + 1; 
end
if md > 0 && md < genop.averageFrames
	genop.vdarr(md) = vd; 
    genop.k = k + 1; 
end
if md == genop.averageFrames
genop.vdarr(md) = vd; 
vds = median(genop.vdarr); 

genop.DMcommandHist = [genop.DMcommandHist genop.DMcommand]; 
genop.DMcommandVd = [genop.DMcommandVd vds]; 
genop.DMcommandK = [genop.DMcommandK k]; 
if 0
    if mean(k-genop.DMcommandK) > 1200
        genop.bleach_correct = genop.bleach_correct * 0.99461; 
		bleach_correct = genop.bleach_correct; 
		bleach_correct
    end
    if mean(k-genop.DMcommandK) < 100 && k > 1000
        genop.bleach_correct = genop.bleach_correct * 1.00137;
		bleach_correct = genop.bleach_correct; 
		bleach_correct
    end
end
agedecay = 1-((k - genop.DMcommandK) / genop.bleach_correct); 
% this, roughly, should mirror the photobleaching rate
% for the given excitation wavelength.
% seems with 7% power at 950nm, it halves over 10k frames. 
% or goes from ~ 14 to 10 over 5k frames, slope of ~= 
[~, indx] = sort(genop.DMcommandVd .* agedecay, 'descend'); 
genop.DMcommandVd = genop.DMcommandVd(indx(1:genop.npop)); 
genop.DMcommandHist = genop.DMcommandHist(:, indx(1:genop.npop)); 
genop.DMcommandK = genop.DMcommandK(:, indx(1:genop.npop)); 

% now we can pick & select a new DM command.
temperature = genop.temperatures(min(k+1, N)); 
pick = floor(rand(1) * (genop.npop-1)) + 1; 
father = genop.DMcommandHist(:,pick); 
pick2 = pick; 
while pick2 == pick
    pick2 = floor(rand(1) * (genop.npop-1)) + 1; 
end
mother = genop.DMcommandHist(:,pick2); 
recomb = randn(97, 1); 
% recomb = zeros(97, 1); % removing recombination strongly affects convergence. 
recomb = (rand(1)-0.5) * 2*pi; 
kid = father .* (genop.dmangle > recomb) + mother .* (genop.dmangle <= recomb); 
% kid = father .* (recomb > 0) + mother .* (recomb < 0); 
noise = (randn(97,1)*temperature) .* (rand(97, 1) > 0.85);
tax = 0; 
% noise = zeros(97, 1); 
% tax = mod(floor(k/200), 27)+3; 
% tax = floor(rand(1) * 82)+4;  
% noise(tax) = randn(1)*temperature; 

genop.DMcommand = reshape(kid+noise, 97, 1); 

% tcmd = zeros(97, 1); 
% tax = floor(k/300)+1
% tcmd(tax) = 4*sin(k/10);
% genop.udps(single([3.1415926; zeros(97,1)])); 
% genop.udps(single([3.1415926; DMcommand_save])); 
genop.udps(single([2.7182818; genop.DMcommand + genop.hidden]));

s = sprintf('%1.5e %1.5e %1.5e %1.5e %1.5e', ...
	genop.DMcommandVd(1), mean(genop.DMcommandVd), temperature, vd, tax); % display the best one. 
disp(s); 
genop.k = k + 1; 
end % if mod

if k > 0
	genop.save_dmcommand(k, :) = genop.DMcommand;
	genop.save_vd(k) = vd;
	genop.save_time(k) = toc; 
end

if mod(k,700) == 699
	% nuke it b/c of bleaching..
	genop.DMcommandVd = zeros(1, genop.npop);
end
% note: 
% if everything is behaving well (disconneted from internet, e.g), then
% running SI at 73hz / 300 x 512 resolution seems to mostly reliably give
% the right intensity reading after two clocks.  
% this means that in the test above, you'll still get on average 4 clocks
% of reliable data out of the 5 in the 50% duty cycle.  
