function [VNullcline, nNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn)
% nullcline calculates nullclines of persistent sodium plus potassium model.
% 
% [VNullcline, nNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn)
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
% gK: double
%   Potassium conductance [nS]
% EK: double
%   Potassium equilibrium potential [mV]
% Vm, Vn: double
% km, kn: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or n
%
% Returns
% -------
% VNullcline: array
%   V-nullcline
% wNullcline: array
%   w-nullcline
%
    mInf = 1 ./ (1 + exp((Vm-V)./km));
    nInf = 1 ./ (1 + exp((Vn-V)./kn));

    VNullcline = (I - gL.*(V-EL) - gNa.*mInf.*(V-ENa)) ./ (gK.*(V-EK));
    nNullcline = nInf;
end