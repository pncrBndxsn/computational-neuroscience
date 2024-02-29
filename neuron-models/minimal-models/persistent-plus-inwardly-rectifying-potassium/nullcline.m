function [VNullcline, nNullcline] = nullcline(V, I, gKir, EK, gK, Vh, kh, Vn, kn)
% nullcline calculates nullclines.
% 
% [VNullcline, hNullcline] = nullcline(V, I, gKir, EK, gK, Vh, kh, Vn, kn)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% I: double
%   External stimulus [pA]
% C: double
%   Membrane capacitance [μF]
% gKir: double
%   Inwardly rectifying potassium conductance [nS]
% EK: double
%   Sodium equilibrium potential [mV]
% gK: double
%   Sodium conductance [nS]
% Vh, Vn: double
% kh, kn: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = h or n
%
% Returns
% -------
% VNullcline: array
%   V-nullcline
% n_null: array
%   n-nullcline
%
    hInf = 1 ./ (1 + exp((Vh-V)./kh));
    nInf = 1 ./ (1 + exp((Vn-V)./kn));

    VNullcline = I ./ (gK.*(V-EK)) - gKir.*hInf ./ gK;
    nNullcline = nInf;
end