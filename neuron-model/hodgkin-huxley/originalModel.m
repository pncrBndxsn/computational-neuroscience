function dXdt = originalModel(X, I, C, gL, EL, gNa, ENa, gK, EK)
% originalModel create a function handle of original Hodgkin-Huxley model.
% 
% dXdt = originalModel(X, I, C, gL, EL, gNa, ENa, gK, EK)
% 
% Parameters
% ----------
% X: array [1, 4]
%   X = [V, m, h, n]
%   V: double
%     Membrane potential [mV]
%   m, h, n: double
%     Gating variable
% I: double
%   External stimulus [pA]
% C: double
%   Membrane capacitance [μF/cm^2]
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
%
% Returns
% -------
% dXdt: array [4, 1]
%   Original Hodgkin-Huxley model
%
    [alphaM, betaM, alphaH, betaH, alphaN, betaN] = gatingVariable(X(1));

    dXdt = zeros(4,1);
    dXdt(1) = (I - gL*(X(1) - EL) - gNa*(X(2)^3)*X(3)*(X(1) - ENa) - gK*(X(4)^4)*(X(1) - EK)) / C;
    dXdt(2) = alphaM*(1 - X(2)) - betaM*X(2);
    dXdt(3) = alphaH*(1 - X(3)) - betaH*X(3);
    dXdt(4) = alphaN*(1 - X(4)) - betaN*X(4);
end