function [VNullcline, hNullcline] = nullcline(V, I, gL, EL, gNa, ENa, Vm, km, Vh, kh)
% nullcline calculates nullclines of transient sodium model.
% 
% [VNullcline, hNullcline] = nullcline(V, I, gL, EL, gNa, ENa, Vm, km, Vh, kh)
% 
% Parameters
% ----------
% V: array
%   membrane potential [mV]
% I: double
%   external stimulus [pA]
% gL: double
%   leakage conductance [nS]
% EL: double
%   resting potential [mV]
% gNa: double
%   sodium conductance [nS]
% ENa: double
%   sodium equilibrium potential [mV]
% Vm, Vh: double
% km, kh: double
%   parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
%
% Returns
% -------
% VNullcline: array
%   V-nullcline
% hNullcline: array
%   h-nullcline
%
    mInf = 1 ./ (1 + exp((Vm-V)./km));
    hInf = 1 ./ (1 + exp((Vh-V)./kh));

    VNullcline = (I - gL.*(V-EL)) ./ (gNa.*(mInf.^3).*(V-ENa));
    hNullcline = hInf;
end