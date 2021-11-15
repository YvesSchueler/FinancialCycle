function [f_hat_save,f_hat_save_separate] = pw_cohesion(dY,grid,factor,idx)
%Purpose: Computes power cohesion
%
%Coherent financial cycles for G-7 countries: Why credit can be an asset
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2018/03/27
dy = dY(:,logical(idx));
M = size(dy,2);



f_hat_save = zeros(length(grid),1);
f_hat_save_separate = [];

kk = 0;
for ii=1:M-1
    for jj=ii+1:M
    kk = kk + 1;
    y=dy(:,ii);
    x=dy(:,jj);
    inv=(x~=x)|(y~=y); % clean missing or invalid data points
    x(inv)=[];
    y(inv)=[];        
    T = length(y);
    gammahat = xcov(y,x)/T;  
    gammahat = gammahat./(std(y)*std(x));
    t_cov = length(gammahat);
    Mwind = round(sqrt(length(y))*factor);                     
    w = parzenwin(2*Mwind+1);
    ti=floor(t_cov/2);
    tt=(ti+Mwind):-1:(ti-Mwind);
    cov_trunc = gammahat(tt,:).*(w);
    f_hat=(2/(2*pi))*abs((cov_trunc'*exp(-1i*(-Mwind:Mwind)'*grid))'); 
    

    f_hat_save = f_hat_save + f_hat;
    f_hat_save_separate = [f_hat_save_separate f_hat.*(std(y)*std(x))];
    end
end

f_hat_save = f_hat_save/kk;



end

