function [dVdt, dndt] = vectorField(V, n, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tauN)
% vectorField calculates vector field of persistent sodium plus potassium model.
% 
% [dVdt, dndt] = vectorField(V, n, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tauN)
% 
% Parameters
% ----------
% V: array
%   Membrane potential [mV]
% n: array
%   K^+ activation variable
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
% tauN: double
%   Time constant of nInf [ms]
%
% Returns
% -------
% dVdt: array
%   Time derivative of V
% dndt: array
%   Time derivative of n
%
    mInf = 1 ./ (1 + exp((Vm-V)./km));
    nInf = 1 ./ (1 + exp((Vn-V)./kn));

    dVdt = (I - gL.*(V-EL) - gNa.*mInf.*(V-ENa) - gK.*n.*(V-EK)) ./ C;
    dndt = (nInf-n) ./ tauN;
end