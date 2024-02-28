function dXdt = persistentPlusInwardlyRectifyingPotassium(X, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tauN)
% persistentPlusInwardlyRectifyingPotassium creates a function handle of persistent plus inwardly rectifying potassium model.
% 
% dXdt = persistentPlusInwardlyRectifyingPotassium(X, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tauN)
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
% dXdt: array
%   Persistent plus inwardly rectifying potassium model
%
    hInf = 1 ./ (1 + exp((Vh-X(1))./kh));
    nInf = 1 ./ (1 + exp((Vn-X(1))./kn));

    dXdt = zeros(2,1);
    dXdt(1) = (I - gKir*hInf*(X(1)-EK) - gK*X(2)*(X(1)-EK)) / C;
    dXdt(2) = (nInf-X(2)) / tauN;
end