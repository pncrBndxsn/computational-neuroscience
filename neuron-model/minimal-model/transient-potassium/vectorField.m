function [dVdt, dmdt] = vectorField(V, m, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tauM)
% vectorField calculates vector field of transient potassium (A-current) model.
% 
% [dVdt, dmdt] = vectorField(V, m, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tauM)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% m: array
%   K^+ activation variable
% I: double
%   External stimulus [pA]
% C: double
%   Membrane capacitance [μF]
% gL: double
%   Leakage conductance [nS]
% EL: double
%   Resting potential [mV]
% gA: double
%   Potassium conductance [nS]
% EK: double
%   Potassium equilibrium potential [mV]
% Vm, Vh: double
% km, kh: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
% tauM: double
%   Time constant of mInf [ms]
%
% Returns
% -------
% dVdt: array
%   Time derivative of V
% dmdt: array
%   Time derivative of m
%
    mInf = 1 ./ (1 + exp((Vm-V)./km));
    hInf = 1 ./ (1 + exp((Vh-V)./kh));

    dVdt = (I - gL.*(V-EL) - gA.*m.*hInf.*(V-EK)) ./ C;
    dmdt = (mInf-m) ./ tauM;
end