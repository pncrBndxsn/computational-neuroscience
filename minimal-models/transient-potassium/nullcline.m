function [VNullcline, mNullcline] = nullcline(V, I, gL, EL, gA, EK, Vm, km, Vh, kh)
% nullcline calculates nullclines of transient potassium (A-current) model.
% 
% [VNullcline, mNullcline] = nullcline(V, I, gL, EL, gA, EK, Vm, km, Vh, kh)
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
% gA: double
%   Potassium conductance [nS]
% EK: double
%   Potassium equilibrium potential [mV]
% Vm, Vh: double
% km, kh: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
%
% Returns
% -------
% VNullcline: array
%   V-nullcline
% mNullcline: array
%   m-nullcline
%
    mInf = 1 ./ (1 + exp((Vm-V)./km));
    hInf = 1 ./ (1 + exp((Vh-V)./kh));

    VNullcline = (I - gL.*(V-EL)) ./ (gA.*hInf.*(V-EK));
    mNullcline = mInf;
end