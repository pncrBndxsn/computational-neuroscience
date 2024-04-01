function [dVdt, dhdt] = vectorField(V, h, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tauH)
% vectorField calculates vector field of transient sodium model.
% 
% [dVdt, dhdt] = vectorField(V, h, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tauH)
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
% Vm, Vh: double
% km, kh: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
% tauH: double
%   Time constant of hInf [ms]
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

    dVdt = (I - gL.*(V-EL) - gNa.*(mInf.^3).*h.*(V-ENa)) ./ C;
    dhdt = (hInf-h) ./ tauH;
end