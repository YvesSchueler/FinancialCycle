function [dY_w_fcast,length_fcast] = fn_forecast_ar(temp)
%Purpose: Forecast to extend sample before filtering
%
%Coherent financial cycles for G-7 countries: Why credit can be an asset
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2018/03/27

[T]=length(temp);
%forecast the financial cycle
length_fcast = min(floor(2*T/3),40);
FC_fcast = zeros(length_fcast,1);
[w_hat, A_hat, ~, ~, ~, ~]=arfit(temp,1,8);
p_opt = length(A_hat);
FC_fcast(1,1) = w_hat + A_hat*flipud(temp(end-p_opt+1:end,1));
if p_opt > 1
for ii = 2:p_opt
FC_fcast(ii,1) = w_hat + A_hat(1,1:ii-1)*flipud(FC_fcast(1:ii-1,1)) + A_hat(1,ii:end)*flipud(temp(end-p_opt+ii:end,1));
end
end
for jj=1+p_opt:length_fcast
    FC_fcast(jj,1) = w_hat + A_hat*flipud(FC_fcast(jj-p_opt:jj-1,1));
end

dY_w_fcast = [temp; FC_fcast];
end