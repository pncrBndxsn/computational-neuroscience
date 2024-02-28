function [dVdt, dhdt] = vectorField(V, h, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, Cbase, Camp, Vmax, sig)
% vectorField calculates vector field of persistent sodium plus h-current model.
% 
% [dVdt, dhdt] = vectorField(V, h, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, Cbase, Camp, Vmax, sig)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% h: array
%   Na^+ inactivation variable
% I: double
%   External stimulus [pA]
% C: double
%   Membrane capacitance [μF]
% gL: double
%   Leakage conductance [nS]
% EL: double
%   Resting potential [mV]
% gNa: double
%   Sodium conductance [nS]
% ENa: double
%   Sodium equilibrium potential [mV]
% gh: double
%   Conductance of Ih [nS]
% Eh: double
%   Equilibrium potential of Ih [mV]
% Vm, Vh: double
% km, kh: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
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
    mInf = 1 ./ (1 + exp((Vm-V)./km));
    hInf = 1 ./ (1 + exp((Vh-V)./kh));
    tauH = Cbase + Camp.*exp(-((Vmax-V)./sig).^2);

    dVdt = (I - gL.*(V-EL) - gNa.*mInf.*(V-ENa) - gh.*h.*(V-Eh)) ./ C;
    dhdt = (hInf-h) ./ tauH;
end