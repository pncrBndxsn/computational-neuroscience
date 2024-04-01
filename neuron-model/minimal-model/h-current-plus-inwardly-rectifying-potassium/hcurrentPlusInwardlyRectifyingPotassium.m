function dXdt = hcurrentPlusInwardlyRectifyingPotassium(X, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, Cbase, Camp, Vmax, sig)
% hcurrentPlusInwardlyRectifyingPotassium creates a function handle of h-current plus inwardly rectifying potassium model.
%
% dXdt = hcurrentPlusInwardlyRectifyingPotassium(X, I, C, gL, EL, gKir, EK, gh, Eh, VhKir, khKir, Vh, kh, Cbase, Camp, Vmax, sig)
%
% Parameters
% ----------
% X: array
%   X = [dVdt, dhdt]
%   dVdt: Time derivative of V
%   dhdt: Time derivative of h
% I: double
%   External stimulus [pA]
% C: double
%   Membrane capacitance [μF]
% gL: double
%   Leakage conductance [nS]
% EL: double
%   Resting potential [mV]
% gKir: double
%   Inwardly rectifying potassium conductance [nS]
% EK: double
%   Sodium equilibrium potential [mV]
% gh: double
%   Conductance of h-current [nS]
% Eh: double
%   Equilibrium potential of h-current [mV]
% VhKir, Vh: double
% khKir, kh: double
%   Parameters of steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = hKir or h
% Cbase, Camp, Vmax, sig: double
%   Parameters for voltage-sensitive time constant
%   tauH = Cbase + Camp.*exp(-((Vmax-V)./sig).^2)
%
% Returns
% -------
% dXdt: array
%   h-current plus inwardly rectifying potassium model
%
    hKirInf = 1 ./ (1 + exp((VhKir-X(1))./khKir));
    hInf = 1 ./ (1 + exp((Vh-X(1))./kh));
    tauH = Cbase + Camp.*exp(-((Vmax-X(1))./sig).^2);

    dXdt = zeros(2,1);
    dXdt(1) = (I - gL*(X(1)-EL) - gKir*hKirInf*(X(1)-EK) - gh*X(2)*(X(1)-Eh)) / C;
    dXdt(2) = (hInf-X(2)) / tauH;
end