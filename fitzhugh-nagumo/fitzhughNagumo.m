function dXdt = fitzhughNagumo(X, I, a, b, c)
% fitzhughNagumo create a function handle of FitzHugh-Nagumo model.
% 
% dXdt = fitzhughNagumo(X, I, a, b, c)
% 
% Parameters
% ----------
% X: array [1, 2]
%   X = [V, w]
%   V: double
%     Membrane potential [mV]
%   w: double
%     Recovery variable
% I: double
%   External stimulus [pA]
%
% Returns
% -------
% dXdt: array [2, 1]
%   FitzHugh-Nagumo model
%
    dXdt = zeros(2,1);
    dXdt(1) = X(1).*(a - X(1)).*(X(1) - 1) - X(2) + I;
    dXdt(2) = b*X(1) - c*X(2);
end