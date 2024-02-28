%% Persistent sodium plus potassium model

close all;
clear;
clc;

%% Parameters for Persistent sodium plus potassium model
I = 5.0;                               % External stimulus [pA]
C = 1.0;                               % Membrane capacitance [μF]
gL =   8.0;  gNa = 20.0;  gK =  10.0;  % Membrane conductance [nS]
EL = -80.0;  ENa = 60.0;  EK = -90.0;  % Resting or equilibrium potential [mV]

% Parameters for steady-state activation (or inactivation) curves
% pInf = 1./(1 + (exp(Vp-V)./kp)), p = m or n
Vm = -20.0;  Vn = -25.0;
km =  15.0;  kn =   5.0;

tauN = 1.0;  % Time constant of nInf

%% Solve Persistent sodium plus potassium model.
tmin = 0.0;  tmax = 100.0;
interval = [tmin tmax];
X0 = [-20.0, 0.2];
dXdt = @(t, x) persistentSodiumPlusPotassium(x, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tauN);
[t1, X1] = ode45(dXdt, interval, X0);

%% Caluculate nullclines of persistent sodium plus potassium model.
xmin = -85.0;  xmax = 20.0;
ymin =  -0.1;  ymax =  0.7;
V = linspace(xmin, xmax, 1000);
[VNullcline, nNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn);

%% Caluculate vector field of persistent sodium plus potassium model.
x = linspace(xmin, xmax, 30);
y = linspace(ymin, ymax, 30);
[X, Y] = meshgrid(x, y);
[dVdt, dndt] = vectorField(X, Y, I, C, gL, EL, gNa, ENa, gK, EK, Vm, km, Vn, kn, tauN);

%% Plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(X, Y, dVdt, dndt, 4, Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, VNullcline, 'k-', LineWidth=2.0);
plot(V, nNullcline, 'k--', LineWidth=2.0);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2.0);
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$ \rm K^+ \ activation, $ \it n', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, FontSize=16);
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2.0);
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, FontSize=16);
grid on;
