%% FitzHugh-Nagumo model

close all;
clear;
clc;

%% Set parameters for FitzHugh-Nagumo model.
I =  0.0;  % External stimulus [pA]
a = -0.1;
b =  0.01;
c =  0.02;

%% Solve FitzHugh-Nagumo model.
tmin = 0.0;  tmax = 300.0;
interval = [tmin tmax];

f = @(t, X) fitzhughNagumo(X, I, a, b, c);
initializeX = [0.0, 0.2];
[t1, X1] = ode45(f, interval, initializeX);

%% Calculate nullclines of FitzHugh-Nagumo model.
xmin = -0.5;  xmax = 1.2;
ymin = -0.1;  ymax = 0.3;
V = linspace(xmin, xmax, 100);
[VNullcline, wNullcline] = nullcline(V, I, a, b, c);

%% Calculate vector field of FitzHugh-Nagumo model.
x = linspace(xmin, xmax, 30);
y = linspace(ymin, ymax, 30);
[X, Y] = meshgrid(x, y);
[dVdt, dwdt] = vectorField(X, Y, I, a, b, c);

%% Plot
figure(1); hold on;
subplot(2,1,1); hold on;
quiver(X, Y, dVdt, dwdt, 4,  Marker='.', Alignment='head', ShowArrowHead='off');
plot(V, VNullcline, 'k-', LineWidth=2.0);
plot(V, wNullcline, 'k--', LineWidth=2.0);
plot(X1(:,1), X1(:,2), 'r-', LineWidth=2.0)
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('Membrane Voltage, $V$ [mV]', Interpreter='latex');
ylabel('$ \rm Recovery, $ \it w', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax,'FontSize',16);
grid on;

subplot(2,1,2); hold on;
plot(t1, X1(:,1), '-', LineWidth=2.0);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('Membrane Voltage, $V$ [mV]', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax,'FontSize',16);
grid on;