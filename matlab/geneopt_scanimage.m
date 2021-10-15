function geneopt_scanimage(src,evt,varargin)

global genop; 

lastStripe = src.hSI.hDisplay.stripeDataBuffer{src.hSI.hDisplay.stripeDataBufferPointer}; % get the pointer to the last acquired stripeData
channels = lastStripe.roiData{1}.channels; % get the channel numbers stored in stripeData

for idx = 1:length(channels)  
    frame{idx} = lastStripe.roiData{1}.imageData{idx}{1}; % extract all channels
    if idx == 1  
        frm = frame{idx}; 
        vd = sum(sum(double(frm) ));   
    end   
end

k = genop.k; 
N = genop.N; 

genop.DMcommandHist = [genop.DMcommandHist genop.DMcommand]; 
genop.DMcommandVd = [genop.DMcommandVd vd]; 
genop.DMcommandK = [genop.DMcommandK k]; 
if 1
    if mean(k-genop.DMcommandK) > 450
        genop.bleach_correct = genop.bleach_correct * 0.99461; 
    end
    if mean(k-genop.DMcommandK) < 100 && k > 1000
        genop.bleach_correct = genop.bleach_correct * 1.00137; 
    end
end
agedecay = 1-((k - genop.DMcommandK) / genop.bleach_correct); 
% this, roughly, should mirror the photobleaching rate
% for the given excitation wavelength.
% seems with 7% power at 950nm, it halves over 10k frames. 
% or goes from ~ 14 to 10 over 5k frames, slope of ~= 
[~, indx] = sort(genop.DMcommandVd .* agedecay, 'descend'); 
genop.DMcommandVd = genop.DMcommandVd(indx(1:100)); 
genop.DMcommandHist = genop.DMcommandHist(:, indx(1:100)); 
genop.DMcommandK = genop.DMcommandK(:, indx(1:100)); 

% now we can pick & select a new DM command.
temperature = genop.temperatures(mod(k-1,N)+1); 
pick = floor(rand(1) * 99) + 1; 
father = genop.DMcommandHist(:,pick); 
pick2 = pick; 
while pick2 == pick
    pick2 = floor(rand(1) * 99) + 1; 
end
mother = genop.DMcommandHist(:,pick2); 
recomb = randn(97, 1); 
kid = father .* (recomb > 0) + mother .* (recomb < 0); 
noise = (randn(97,1)*temperature) .* (rand(97, 1) > 0.83); 
genop.DMcommand = reshape(kid+noise, 97, 1); 

genop.udps(single([3.1415926; genop.DMcommand])); 

genop.save_dmcommand(k, :) = genop.DMcommand;
genop.save_vd(k) = vd;
genop.save_time(k) = toc; 

disp([genop.DMcommandVd(1) mean(genop.DMcommandVd) temperature*1e7 vd]); % display the best one. 

genop.k = k + 1; 
