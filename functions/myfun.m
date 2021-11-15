%Coherent financial cycles for G-7 countries: Why credit can be an asset
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2018/03/27
function f = myfun(x,FCyc_s,FCyc_us)
f = sum(((FCyc_s*x(1))+x(2) - FCyc_us).^2);
end