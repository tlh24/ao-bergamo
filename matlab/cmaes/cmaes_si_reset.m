function cmaes_si_reset(src,evt,varargin)

gobal cm; 

cm.udps = dsp.UDPSender('RemoteIPPort', 13131, ...
    'RemoteIPAddress','10.123.1.152'); 
 
cmaes_init; 

cm.k = 0; 
cm.averageFrames = 4; 