function [alphaM, betaM, alphaH, betaH, alphaN, betaN] = gatingVariable(V)
% gatingVariable calculates gating variables, m, h and n.
% 
% [alphaM, betaM, alphaH, betaH, alphaN, betaN] = gatingVariable(V)
% 
% Parameters
% ----------
% V: array [L, 1]
%   Membrane potential [mV]
%
% Returns
% -------
% alphaM, betaM: array [length(V), 1]
%   Rate constant of sodium channel
% alphaH, betaH: array [length(V), 1]
%   Rate constant of sodium channel
% alphaN, betaN: array [length(V), 1]
%   Rate constant of potassium channel
%
    alphaM = 0.1*(25 - V) ./ (exp((25-V)/10) - 1);
    betaM  = 4*exp(-V/18);
    alphaH = 0.07*exp(-V/20);
    betaH  = 1 ./ (exp((30-V)/10) + 1);
    alphaN = 0.01.*(10 - V) ./ (exp((10-V)/10) - 1);
    betaN  = 0.125.*exp(-V/80);
end