%% transient potassium (A-current) model

close all;
clear;
clc;

%% Parameters for transient potassium (A-current) model
I = 10.6;                 % External stimulus [pA]
C = 1.0;                  % Membrane capacitance [μF]
gL =   0.2;  gA =   5.0;  % Membrane conductance [nS]
EL = -60.0;  EK = -80.0;  % Resting or equilibrium potential [mV]

% Parameters for steady-state activation curves
% pInf = 1./ (1 + (exp(Vp-V)./kp)), p = m or h
Vm = -45.0;  Vh = -66.0;
km =  10.0;  kh = -10.0;

tauM = 20.0;  % Time constant of mInf [ms]

%% Solve transient potassium (A-current) model.
tmin = 0.0;  tmax = 500.0;
interval = [tmin tmax];
X0 = [-25.0, 0.9];
f = @(t, x) transientPotassium(x, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tauM);
[t1, X1] = ode45(f, interval, X0);

%% Caluculate nullclines of transient potassium (A-current) model.
xmin = -80.0;  xmax = 0.0;
ymin =   0.0;  ymax = 1.0;
V = linspace(xmin,xmax,1000)';
[VNullcline, mNullcline] = nullcline(V, I, gL, EL, gA, EK, Vm, km, Vh, kh);

%% Caluculate vector field of transient potassium (A-current) model.
x = linspace(xmin, xmax, 30);
y = linspace(ymin, ymax, 30);
[X, Y] = meshgrid(x, y);
[dVdt, dmdt] = vectorField(X, Y, I, C, gL, EL, gA, EK, Vm, km, Vh, kh, tauM);

%% Plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(X, Y, dVdt, dmdt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, VNullcline, 'k-', LineWidth=2.0);
plot(V, mNullcline, 'k--', LineWidth=2.0);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2.0)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$ \rm K^+ \ activation, $ \it m', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2.0)
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;
