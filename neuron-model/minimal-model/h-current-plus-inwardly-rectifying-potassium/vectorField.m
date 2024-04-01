function [dVdt, dhdt] = vectorField(V, h, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, Cbase, Camp, Vmax, sig)
% vectorField calculates vector field of h-current plus inwardly rectifying potassium model.
% 
% [dVdt, dhdt] = vectorField(V, h, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, Cbase, Camp, Vmax, sig)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% h: array
%   Inactivation variable
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
%   Conductance of Ih [nS]
% Eh: double
%   Equilibrium potential of Ih [mV]
% VhKir, Vh: double
% khKir, kh: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = hKir or h
% Cbase, Camp, Vmax, sig: double
%   Parameters for voltage-sensitive time constant [ms]
%   tauH = Cbase + Camp.*exp(-((Vmax-V)./sig).^2)
% 
% Returns
% -------
% dVdt: array
%   Time derivative of V
% dhdt: array
%   Time derivative of h
%
    hKirInf = 1 ./ (1 + exp((VhKir-V)./khKir));
    hInf = 1 ./ (1 + exp((Vh-V)./kh));
    tauH = Cbase + Camp.*exp(-((Vmax-V)./sig).^2);

    dVdt = (I - gL.*(V-EL) - gKir.*hKirInf.*(V-EK) - gh.*h.*(V-Eh)) ./ C;
    dhdt = (hInf-h) ./ tauH;
end