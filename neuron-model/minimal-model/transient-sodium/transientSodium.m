function dXdt = transientSodium(X, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tauH)
% transientSodium creates a function handle of transient sodium model.
% 
% dXdt = transientSodium(X, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tauH)
% 
% Parameters
% ----------
% X: array
%   X = [V, h]
%   V: Membrane potential [mV]
%   h: Na^+ inactivation variable
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
% dXdt: array
%   Transient sodium model
%
    mInf = 1 ./ (1 + exp((Vm-X(1))./km));
    hInf = 1 ./ (1 + exp((Vh-X(1))./kh));

    dXdt = zeros(2,1);
    dXdt(1) = (I - gL*(X(1)-EL) - gNa*(mInf^3)*X(2)*(X(1)-ENa)) / C;
    dXdt(2) = (hInf-X(2)) / tauH;
end