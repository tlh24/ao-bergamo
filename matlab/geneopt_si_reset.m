function geneopt_si_reset(src,evt,varargin)

global genop; 

genop.udps = dsp.UDPSender('RemoteIPPort', 13131, ...
    'RemoteIPAddress','10.123.1.152'); 

N = 4e4; 
genop.N = N;   
genop.bleach_correct = 30000; 
genop.npop = 25; 

genop.save_dmcommand = single(zeros(N, 97));
genop.save_vd = single(zeros(N, 1)); 
genop.save_time = single(zeros(N, 1)); 

genop.DMcommand = zeros(97, 1);
genop.DMcommandHist = repmat(genop.DMcommand, 1, genop.npop); 
genop.DMcommandVd = zeros(1, genop.npop);
genop.DMcommandK = zeros(1, genop.npop);
if 0
	genop.DMcommand = DMcommand_save; 
	genop.DMcommandHist = DMcommand_Hist; 
end
genop.starttemp = 0.37; % start from zero DM command: 0.005
genop.endtemp = 0.07; 
genop.temperatures = linspace(genop.starttemp, genop.endtemp, N); 
genop.k = 0; 
genop.vs = 0;
genop.averageFrames = 9; 
genop.vdarr = zeros(genop.averageFrames, 1);
genop.avstd = 0.0; 
genop.hidden = [randn(20, 1); zeros(77, 1)] * 0.65; 
genop.hidden = zeros(97, 1); 

genop.mask = zeros(300, 512); 
for r = 1:300
	for c = 1:512
		if abs(r-150) < 50 && abs(c-256) < 60
			genop.mask(r,c) = 1; 
		end
	end
end
tic; 