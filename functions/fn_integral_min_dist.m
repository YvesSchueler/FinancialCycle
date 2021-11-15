function [min_distance]=fn_integral_min_dist(coh,peak,grid,pct)
% computes the perctentage of the area below the spectral density with
% minimal distance.

%input
%coh: cohesion mesure (vector)
%peak: max value around which to select the area
%grid: grid of the periodogram
%pct: percentage of area

%output
%min_distance: place in grid (ii and ii+jj), years (y_high, y_low),
%              distance in grid-points, area share

%Coherent financial cycles for G-7 countries: Why credit can be an asset
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2018/03/27

total_area = trapz(coh);

output_summary = zeros(size(grid,1)-1,6);
for ii = 1:size(grid,2)-1
    jj = 1;
    share = trapz(coh(ii:ii+jj));
    while share/total_area <= pct && ii+jj ~= size(grid,2) 
        jj = jj+1;
        share = trapz(coh(ii:ii+jj));
    end
output_summary(ii,:) = [ii jj+ii (grid(ii)).^(-1)*pi/2 (grid(jj+ii)).^(-1)*pi/2 jj share/total_area];
    
end

temp = output_summary(find(output_summary(:,end) > pct),:); %filter via the pct rule (algorithm also yields areas smaller than pct at the end points

temp = temp(find(temp(:,3) >= (peak)),:); %make sure peak is included
temp = temp(find(temp(:,4) <= (peak)),:);

min_distance = temp(find(temp(:,end-1) == min(temp(:,end-1))),:); %select minimum distance in grid points

%chances that distance is the same... 

min_distance = min_distance(end,:); %take the one with highest frequency



end

