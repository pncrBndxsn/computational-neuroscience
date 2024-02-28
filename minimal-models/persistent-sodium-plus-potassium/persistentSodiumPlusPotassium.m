function dXdt = persistentSodiumPlusPotassium(X, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tauN)
% persistentSodiumPlusPotassium creates a function handle of persistent sodium plus potassium model.
% 
% dXdt = persistentSodiumPlusPotassium(X, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tauN)
% 
% Parameters
% ----------
% X: array
%   X = [V, n]
%   V: Membrane potential [mV]
%   n: K^+ activation variable
% I: double
%   External stimulus [pA]
% C: double
%   Eembrane capacitance [μF]
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
% dXdt: array
%   Persistent sodium plus potassium model
%
    mInf = 1 ./ (1 + exp((Vm-X(1))./km));
    nInf = 1 ./ (1 + exp((Vn-X(1))./kn));

    dXdt = zeros(2,1);
    dXdt(1) = (I - gL*(X(1)-EL) - gNa*mInf*(X(1)-ENa) - gK*X(2)*(X(1)-EK)) / C;
    dXdt(2) = (nInf - X(2)) / tauN;
end