function [VNullcline, hNullcline] = nullcline(V, I, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh)
% nullcline calculates nullclines of h-current plus inwardly rectifying potassium model.
% 
% [VNullcline, hNullcline] = nullcline(V, I, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% I: double
%   External stimulus [pA]
% C: double
%   Membrane capacitance [μF]
% gL: double
%   Leakage conductance [nS]
% EL: double
%   Resting potential [mV]
% gKir: double
%   Inwardly rectifying potassium conductance [nS]
% EK: double
%   Sodium equilibrium potential [mV]
% gh: double
%   Conductance of h-current [nS]
% Eh: double
%   Equilibrium potential of h-current [mV]
% VhKir, Vh: double
% khKir, kh: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = hKir or h
% 
% Returns
% -------
% VNullcline: array
%   V-nullcline
% hNullcline: array
%   h-nullcline
%
    hKirInf = 1 ./ (1 + exp((VhKir-V)./khKir));
    hInf    = 1 ./ (1 + exp((Vh-V)./kh));

    VNullcline = (I - gL.*(V-EL) - gKir.*hKirInf.*(V-EK)) ./ (gh.*(V-Eh));
    hNullcline = hInf;
end