%% persistent plus inwardly rectifying potassium model

close all;
clear;
clc;

%% Parameters for persistent plus inwardly rectifying potassium model
I = 68.0;                 % External stimulus [pA]
C = 1.0;                  % Membrane capacitance [μF]
gKir =  20.0;  gK = 2.0;  % Membrane conductance [nS]
EK   = -80.0;             % Potassium equilibrium potential [mV]

% Parameters for steady-state activation curves
% pInf = 1 ./ (1 + (exp(Vp-V)./kp)), p = h or n
Vh = -80.0;  Vn = -40.0;
kh = -12.0;  kn =   5.0;

tauN = 5.0;  % time constant of nInf [ms]

%% Solve persistent plus inwardly rectifying potassium model.
tmin = 0.0;  tmax = 100.0;
interval = [tmin tmax];
X0 = [-60.0, 0.0];
dXdt = @(t, x) persistentPlusInwardlyRectifyingPotassium(x, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tauN);
[t1, X1] = ode45(dXdt, interval, X0);

%% Caluculate nullclines of persistent plus inwardly rectifying potassium model.
xmin = -80.0;  xmax = 30.0;
ymin =  -0.1;  ymax =  1.0;
V = linspace(xmin, xmax, 1000);
[VNullcline, nNullcline] = nullcline(V, I, gKir, EK, gK, Vh, kh, Vn, kn);

%% Caluculate vector field of persistent plus inwardly rectifying potassium model.
x = linspace(xmin, xmax, 30);
y = linspace(ymin, ymax, 30);
[X, Y] = meshgrid(x, y);
[dVdt, dndt] = vectorField(X, Y, I, C, gKir, EK, gK, Vh, kh, Vn, kn, tauN);

%% Plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(X, Y, dVdt, dndt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, VNullcline, 'k-', LineWidth=2);
plot(V, nNullcline, 'k--', LineWidth=2);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ylabel('$ \rm K^+ \ activation, $ \it n', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2)
xlim([tmin tmax]);
ylim([xmin xmax]);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('Membrane Voltage, $ V $ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;
