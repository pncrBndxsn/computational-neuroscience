function [dVdt, dwdt] = vectorField(V, w, I, a, b, c)
% vectorField calculates vector field of FitzHugh-Nagumo model.
% 
% [dVdt, dwdt] = vectorField(V, I, a, b, c)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% w: array
%   Recovery variable
% I: double
%   External stimulus [pA]
% a: double
% b: double
% c: double
% 
% Returns
% -------
% dVdt: array
%   Time derivative of V
% dwdt: array
%   Time derivative of w
%
    dVdt = V.*(a - V).*(V - 1) - w + I;
    dwdt = b*V - c*w;
end