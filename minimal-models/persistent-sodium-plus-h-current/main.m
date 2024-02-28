%% Persistent sodium plus h-current model

close all;
clear;
clc;

%% Parameters for persistent sodium plus h-current model
I = -1.0;                              % External stimulus [pA]
C = 1.0;                               % Membrane capacitance [μF]
gL =   1.3;  gNa =  0.9;  gh =   3.0;  % Membrane conductance [nS]
EL = -80.0;  ENa = 20.0;  Eh = -43.0;  % Resting or equilibrium potential [mV]

% Parameters for steady-state activation curves
% pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
Vm = -54.0;  Vh = -75.0;
km =   9.0;  kh =  -5.5;

% Parameters for voltage-sensitive time constant [ms]
% tauH = Cbase + Camp.*exp(-((Vmax-V)./sig).^2)
Cbase =  100.0;
Camp  = 1000.0;
Vmax  =  -75.0;
sig   =   15.0;

%% Solve persistent sodium plus h-current model.
tmin = 0.0;  tmax = 1000.0;
interval = [tmin tmax];
X0 = [-60.0, 0.04];
dXdt = @(t, x) persistentSodiumPlusHcurrent(x, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, Cbase, Camp, Vmax, sig);
[t1, X1] = ode45(dXdt, interval, X0);

%% Caluculate nullcline of persistent sodium plus h-current model.
xmin = -100.0;  xmax = -44.0;
ymin =   -0.1;  ymax =   1.0;
V = linspace(xmin, xmax, 1000)';
[VNullcline, hNullcline] = nullcline(V, I, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh);

%% Caluculate vector field of persistent sodium plus h-current model.
x = linspace(xmin, xmax, 30);
y = linspace(ymin, ymax, 30);
[X, Y] = meshgrid(x, y);
[dVdt, dhdt] = vectorField(X, Y, I, C, gL, EL, gNa, ENa, gh, Eh, Vm, km, Vh, kh, Cbase, Camp, Vmax, sig);

%% Plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(X, Y, dVdt, dhdt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, VNullcline, 'k-', LineWidth=2);
plot(V, hNullcline, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2);
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$ \rm inactivation, $ \it h', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16, YDir='reverse');
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2);
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;
