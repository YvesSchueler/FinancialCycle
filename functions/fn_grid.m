function [grid,step ] = fn_grid(q_low,q_high,resolution)
%produces the grid along which the spectral density is estimated
%step: defines x-axis for graphs

%Coherent financial cycles for G-7 countries: Why credit can be an asset
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2018/03/27
omega_step = 2048;

if q_low==0 && q_high ~=0

grid = (2*pi/(q_high)):2*pi/omega_step:pi-pi/(omega_step);
    
elseif q_low ~=0 && q_high ==0

grid = 0:2*pi/omega_step:(2*pi)/q_low;
    
elseif q_low==0 && q_high ==0

grid = 0:2*pi/omega_step:pi-pi/(omega_step);
    
else
    
grid = (2*pi/(q_high)):2*pi/omega_step:(2*pi/q_low);

end


step = (max(grid) - min(grid))/resolution;                    % define steps for frequency display in graphs


end

