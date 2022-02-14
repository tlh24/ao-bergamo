% convert actuator no to XY coordinates.
% I coded this up in C ... guess have to do it again. 
function [dmx,dmy] = dm_actuator_to_xy()

x = ones(5,1)*5; 
x = [x; ones(7,1) * 4]; 
x = [x; ones(9,1) * 3]; 
x = [x; ones(11,1) * 2]; 
x = [x; ones(11,1) * 1]; 
x = [x; ones(11,1) * 0]; 
x = [x; ones(11,1) * -1]; 
x = [x; ones(11,1) * -2]; 
x = [x; ones(9,1) * -3]; 
x = [x; ones(7,1) * -4]; 
x = [x; ones(5,1) * -5]; 

y = [-2:2];
y = [y  [-3:3]]; 
y = [y  [-4:4]]; 
y = [y  [-5:5]]; 
y = [y  [-5:5]]; 
y = [y  [-5:5]]; 
y = [y  [-5:5]]; 
y = [y  [-5:5]]; 
y = [y  [-4:4]];
y = [y  [-3:3]]; 
y = [y  [-2:2]]; 

dmx = x; 
dmy = y'; 