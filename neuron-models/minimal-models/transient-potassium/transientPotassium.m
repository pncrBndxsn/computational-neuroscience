function dXdt = transientPotassium(X, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tauM)
% transientPotassium creates a function handle of transient potassium (A-current) model.
% 
% dXdt = transientPotassium(X, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tauM)
% 
% Parameters
% ----------
% X: array
%   X = [dVdt, dhdt]
%   dVdt: Time derivative of V
%   dmdt: Time derivative of m
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
% dXdt: array
%   transient potassium (A-current) model
%
    mInf = 1 ./ (1 + exp((Vm-X(1))./km));
    hInf = 1 ./ (1 + exp((Vh-X(1))./kh));

    dXdt = zeros(2,1);
    dXdt(1) = (I - gL*(X(1)-EL) - gA*X(2)*hInf*(X(1)-EK)) / C;
    dXdt(2) = (mInf-X(2)) / tauM;
end