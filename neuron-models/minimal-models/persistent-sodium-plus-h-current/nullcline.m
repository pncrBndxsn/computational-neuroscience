function [VNullcline, hNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh)
% nullcline calculates nullclines of persistent sodium plus h-current model.
% 
% [VNullcline, hNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh)
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
% gNa: double
%   Sodium conductance [nS]
% ENa: double
%   Sodium equilibrium potential [mV]
% gh: double
%   Conductance of h-current [nS]
% Eh: double
%   Equilibrium potential of h-current [mV]
% Vm, Vh: double
% km, kh: double
%   Parameters for steady-state activation (or inactivation) curves
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

    VNullcline = (I - gL.*(V-EL) - gNa.*mInf.*(V-ENa)) ./ (gh.*(V-Eh));
    hNullcline = hInf;
end