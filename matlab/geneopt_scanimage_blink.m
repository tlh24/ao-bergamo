function geneopt_scanimage_blink(src,evt,varargin)

global genop; 

k = genop.k; 
N = genop.N; 


genop.DMcommand = genop.Best_DMcommand; 
% if mod(floor(k/45), 2) == 0
% 	genop.DMcommand = genop.DMcommand_save; 
% end

genop.udps(single([2.7182818;; genop.DMcommand]));

genop.k = k+1; 
end % else
