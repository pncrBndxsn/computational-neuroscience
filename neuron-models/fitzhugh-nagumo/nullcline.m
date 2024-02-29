function [VNullcline, wNullcline] = nullcline(V, I, a, b, c)
% nullcline calculates nullclines of FitzHugh-Nagumo model.
% 
% [VNullcline, nNullcline] = nullcline(V, I, a, b, c)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% I: double
%   External stimulus [pA]
% a: double
% b: double
% c: double
% 
% Returns
% -------
% VNullcline: array
%   V-nullcline
% wNullcline: array
%   w-nullcline
%
    VNullcline = V.*(a - V).*(V - 1) + I;
    wNullcline = (b/c).*V;
end