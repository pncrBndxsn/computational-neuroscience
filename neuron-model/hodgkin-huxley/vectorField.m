function [dVdt, dndt] = vectorField(V, n, I, C, gL, EL, gNa, ENa, gK, EK, p)
% vectorField calculates vector field of reduced Hodgkin-Huxley model.
% 
% [dVdt, dndt] = vectorField(V, n, I, C, gL, EL, gNa, ENa, gK, EK, p)
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
% p: array [1, 2]
%   Least-squares fit polynomial coefficients
%
% Returns
% -------
% dVdt: array
%   Time derivative of V
% dndt: array
%   Time derivative of n
%
    [alphaM, betaM, ~, ~, alphaN, betaN] = gatingVariable(V);
    tauM = 1./(alphaM + betaM);  mInf = alphaM .* tauM;
    tauN = 1./(alphaN + betaN);  nInf = alphaN .* tauN;

    dVdt = (I - gL.*(V - EL) - gNa.*(mInf.^3).*(p(2) + p(1).*n).*(V - ENa) - gK*(n.^4).*(V - EK)) ./ C;
    dndt = (nInf - n) ./ tauN;
end