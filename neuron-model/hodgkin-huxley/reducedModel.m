function dXdt = reducedModel(X, I, C, gL, EL, gNa, ENa, gK, EK, p)
% reducedModel create a function handle of reduced Hodgkin-Huxley model.
% 
% dXdt = reducedModel(X, I, C, gL, EL, gNa, ENa, gK, EK, p)
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
% p: array [1, 2]
%   Least-squares fit polynomial coefficients
%
% Returns
% -------
% dXdt: array [4, 1]
%   Reduced Hodgkin-Huxley model
%
    [alphaM, betaM, alphaH, betaH, alphaN, betaN] = gatingVariable(X(1));
    tauM = 1./(alphaM + betaM);  mInf = alphaM .* tauM;
    tauH = 1./(alphaH + betaH);  hInf = alphaH .* tauH;
    tauN = 1./(alphaN + betaN);  nInf = alphaN .* tauN;

    dXdt = zeros(4,1);
    dXdt(1) = (I - gL*(X(1) - EL) - gNa*(mInf^3)*(p(2) + p(1)*X(4))*(X(1) - ENa) - gK*(X(4)^4)*(X(1) - EK)) / C;
    dXdt(2) = (mInf - X(2))/tauM;
    dXdt(3) = (hInf - X(3))/tauH;
    dXdt(4) = (nInf - X(4))/tauN;
end