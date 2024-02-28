function dXdt = persistentSodiumPlusHcurrent(X, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, Cbase, Camp, Vmax, sig)
% persistentSodiumPlusHcurrent create a function handle of persistent sodium plus h-current model.
% 
% dXdt = persistentSodiumPlusHcurrent(X, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, Cbase, Camp, Vmax, sig)
%
% Parameters
% ----------
% X: array
%   X = [dVdt, dhdt]
%   dVdt: time derivative of V
%   dhdt: time derivative of h
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
% gh: double
%   Conductance of h-current [nS]
% Eh: double
%   Equilibrium potential of h-current [mV]
% Vm, Vh: double
% km, kh: double
%   Parameters for steady-state activation (or inactivation) curves
%   pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
% Cbase, Camp, Vmax, sig: double
%   Parameters for voltage-sensitive time constant [ms]
%   tauH = Cbase + Camp.*exp(-((Vmax-V)./sig).^2)
%
% Returns
% -------
% dXdt: array
%   Persistent sodium plus h-current model
%
    mInf = 1 ./ (1 + exp((Vm-X(1))./km));
    hInf = 1 ./ (1 + exp((Vh-X(1))./kh));
    tauH = Cbase + Camp.*exp(-((Vmax-X(1))./sig).^2);

    dXdt = zeros(2,1);
    dXdt(1) = (I - gL*(X(1)-EL) - gNa*mInf*(X(1)-ENa) - gh*X(2)*(X(1)-Eh)) / C;
    dXdt(2) = (hInf-X(2)) / tauH;
end