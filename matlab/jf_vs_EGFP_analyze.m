egfp = zeros(20, 4);
jf = egfp; 

load jf669_vs_EGFP_940_0.mat
egfp(:,1) = dat(:,1); 
load jf669_vs_EGFP_940_6.mat
egfp(:,2) = dat(:,1); 
load jf669_vs_EGFP_940_10.mat
egfp(:,3) = dat(:,1); 
load jf669_vs_EGFP_940_20.mat
egfp(:,4) = dat(:,1); 

load jf669_vs_EGFP_1130_0.mat
jf(:,1) = dat(:,1); 
load jf669_vs_EGFP_1130_6.mat
jf(:,2) = dat(:,1); 
load jf669_vs_EGFP_1130_10.mat
jf(:,3) = dat(:,1); 
load jf669_vs_EGFP_1130_20.mat
jf(:,4) = dat(:,1); 

x = [0 0.5 1.0 2.0]; 
data = [mean(egfp, 1) ; mean(jf, 1)]'; 
se = [std(egfp, 1) ; std(jf, 1)]';
bar(data)
hold on
er = errorbar(data, 2*se);
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 