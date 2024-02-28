function [dVdt, dndt] = vectorField(V, n, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tauN)
% vectorField calculates vector field.
% 
% [dVdt, dndt] = vectorField(V, n, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tauN)
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
% gKir: double
%   Inwardly rectifying potassium conductance [nS]
% EK: double
%   Sodium equilibrium potential [mV]
% gK: double
%   Sodium conductance [nS]
% Vh, Vn: double
% kh, kn: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = h or n
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
    hInf = 1 ./ (1 + exp((Vh-V)./kh));
    nInf = 1 ./ (1 + exp((Vn-V)./kn));

    dVdt = (I - gKir.*hInf.*(V-EK) - gK.*n.*(V-EK)) ./ C;
    dndt = (nInf-n) ./ tauN;
end