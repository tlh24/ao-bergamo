function geneopt_si_reset(src,evt,varargin)

global genop; 

genop.udps = dsp.UDPSender('RemoteIPPort', 13131, ...
    'RemoteIPAddress','10.123.1.152'); 

N = 10e3; 
genop.N = N;   
genop.bleach_correct = 7000; 

genop.save_dmcommand = single(zeros(N, 97));
genop.save_vd = single(zeros(N, 1)); 
genop.save_time = single(zeros(N, 1)); 

genop.DMcommand = zeros(97, 1);
genop.DMcommandHist = repmat(genop.DMcommand, 1, 100); 
genop.DMcommandVd = zeros(1, 100);
genop.DMcommandK = zeros(1, 100);

genop.starttemp = 0.007; % start from zero DM command: 0.005
genop.endtemp = 0.002; 
genop.temperatures = linspace(genop.starttemp, genop.endtemp, N); 
genop.k = 1; 
tic; 