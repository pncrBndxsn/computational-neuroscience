function [VNullcline, nNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, p)
% nullcline calculates nullclines of reduced Hodgkin-Huxley model.
% 
% [VNullcline, nNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, p)
% 
% Parameters
% ----------
% V: double
%   Membrane potential [mV]
% I: double
%   External stimulus [pA]
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
% VNullcline: double
%   V-nullcline
% nNullcline: double
%   n-nullcline
%
    syms n

    [alphaM, betaM, ~, ~, alphaN, betaN] = gatingVariable(V);
    tauM = 1/(alphaM + betaM);  mInf = alphaM * tauM;
    tauN = 1/(alphaN + betaN);  nInf = alphaN * tauN;

    solutionHH = vpasolve(I - gL*(V-EL) - gNa*(mInf^3)*(p(2)+p(1)*n)*(V-ENa) - gK*(n^4)*(V-EK) == 0, n, [-Inf Inf]);
    VNullcline = double(max(solutionHH));
    nNullcline = nInf;
end