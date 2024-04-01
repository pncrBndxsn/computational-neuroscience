%% Original and reduced Hodgkin-Huxley model

close all;
clear;
clc;

%% Set parameters of Hodgkin-Huxley model.
I1 = 10.0;  % External stimulus of original Hodgkin-Huxley [pA]
I2 = 5.0;  % External stimulus of reduced Hodgkin-Huxley [pA]
C = 1.0;  % Membrane capacitance [μF]
gL =  0.3;  gNa = 120.0;  gK =  36.0;  % Membrane conductance [nS]
EL = 10.6;  ENa = 120.0;  EK = -12.0;  % Resting, equilibrium potential [mV]

%% Integral interval
tmin = 0.0;  tmax = 50.0;
interval = [tmin tmax];

%% Solve original Hodgkin-Huxley model.
f = @(t, X) originalModel(X, I1, C, gL, EL, gNa, ENa, gK, EK);
initializeX1 = [0.0, 0.1, 0.6, 0.3];

[t1, X1] = ode45(f, interval, initializeX1);

%% Linear approximation by least squares method
p = polyfit(X1(:,4), X1(:,3), 1);

%% Solve reduced Hodgkin-Huxley model.
g = @(t, X) reducedModel(X, I2, C, gL, EL, gNa, ENa, gK, EK, p);
initializeX2 = [0.0, 0.1, 0.6, 0.3];

[t2, X2] = ode45(g, interval, initializeX2);

%% Calculate nullclines of reduced Hodgkin-Huxley model.
xmin = -11.0;  xmax = 119;
ymin =   0.0;  ymax = 1.0;
V = linspace(xmin, xmax, 200);

VNullcline = zeros(length(V), 1);
nNullcline = zeros(length(V), 1);
for i = 1:length(V)
    [VNullcline(i), nNullcline(i)] = nullcline(V(i), I2, gL, EL, gNa, ENa, gK, EK, p);
end

%% Calculate vector field of reduced Hodgkin-Huxley model.
x = linspace(xmin, xmax, 30);
y = linspace(ymin, ymax, 30);
[X, Y] = meshgrid(x, y);
[dVdt, dndt] = vectorField(X, Y, C, I2, gL, EL, gNa, ENa, gK, EK, p);

%% Plot
figure(1); hold on;
subplot(3,1,1); hold on;
plot(t1, X1(:,1), '-', LineWidth=2.0);
ylabel('Membrane Voltage, $V$ [mV]', Interpreter='latex');
title('Original Hodgkin-Huxley', Interpreter='latex')
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;

subplot(3,1,2); hold on;
plot(t1, X1(:,2), 'r-', LineWidth=2.0);
plot(t1, X1(:,3), 'g-', LineWidth=2.0);
plot(t1, X1(:,4), 'b-', LineWidth=2.0);
ylabel('Gating Variables', Interpreter='latex');
legend('$m$', '$h$', '$n$', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;

subplot(3,1,3); hold on;
plot(t1, X1(:,3)+X1(:,4), '-', LineWidth=2.0);
ylim([0.0 2.0]);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('$h + n$', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;


figure(2); hold on;
subplot(1,1,1); hold on;
plot(X1(:,4), X1(:,3), '.', LineWidth=2.0);
plot(X1(:,4), p(2)+p(1)*X1(:,4), '-', LineWidth=2.0);
xlabel('$ \rm K^+ \ activation, $ \it n', Interpreter='latex');
ylabel('$ \rm Na^+ \ inactivation, $ \it h', Interpreter='latex');
text(0.45, 0.5, ['$h=-$',num2str(abs(p(1))),'$n+$',num2str(abs(p(2)))], FontSize=16, Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;
axis equal


figure(3); hold on;
subplot(1,1,1); hold on;
quiver(X, Y, dVdt, dndt, 4, Marker='.', Alignment='head', ShowArrowHead='off');
plot(X2(:,1), X2(:,4), '-', LineWidth=2.0);
plot(V, VNullcline, 'k-', LineWidth=2.0);
plot(V, nNullcline, 'k--', LineWidth=2.0);
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('Membrane Voltage, $V$ [mV]', Interpreter='latex');
ylabel('$ \rm K^+ \ activation, $ \it n', Interpreter='latex');
legend('', 'Reduced', Interpreter='latex');
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;


figure(4); hold on;
subplot(1,1,1); hold on;
plot(t1, X1(:,1), '-', LineWidth=2.0);
plot(t2, X2(:,1), '--', LineWidth=2.0);
xlabel('$t$ [ms]', Interpreter='latex');
ylabel('Membrane Voltage, $V$ [mV]', Interpreter='latex');
legend('Original','Reduced', Interpreter='latex')
ax = gca;
ax.TickLabelInterpreter='latex';
set(ax, FontSize=16);
grid on;