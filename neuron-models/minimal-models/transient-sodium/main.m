%% Transient sodium model

close all;
clear;
clc;

%% Parameters for transient sodium model
I = 0.5;                  % External stimulus [pA]
C = 1.0;                  % Membrane capacitance [μF]
gL =   1.0;  gNa = 10.0;  % Membrane conductance [nS]
EL = -70.0;  ENa = 60.0;  % Resting or equilibrium potential [mV]

% Parameters for steady-state activation (or inactivation) curves
% pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
Vm = -40.0;  Vh = -42.0;
km =  15.0;  kh =  -7.0;

tauH = 5.0;  % Time constant of hInf [ms]

%% Solve transient sodium model.
tmin = 0.0;  tmax = 100.0;
interval = [tmin tmax];
X0 = [-50.0, 0.8];
dXdt = @(t, x) transientSodium(x, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tauH);
[t1, X1] = ode45(dXdt, interval, X0);

%% Caluculate nullclines of transient sodium model.
xmin = -70.0;  xmax = 50.0;
ymin =   0.0;  ymax =  1.0;
V = linspace(xmin, xmax, 1000)';
[VNullcline, hNullcline] = nullcline(V, I, gL, EL, gNa, ENa, Vm, km, Vh, kh);

%% Caluculate vector field of transient sodium model.
x = linspace(xmin, xmax, 30);
y = linspace(ymin, ymax, 30);
[X, Y] = meshgrid(x, y);
[dVdt, dhdt] = vectorField(X, Y, I, C, gL, EL, gNa, ENa, Vm, km, Vh, kh, tauH);

%% Plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(X, Y, dVdt, dhdt, 4, Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, VNullcline, 'k', LineWidth=2);
plot(V, hNullcline, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('Membrane voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$ \rm Na^+ \ inactivation, $ \it h', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16, YDir='reverse');
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), LineWidth=2)
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('Membrane voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;
