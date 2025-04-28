% filepath: c:\Users\nc24951\Downloads\numerical_simulations_M2.m
clc;
clear;

% Define helper functions
function val = h1_star(s2)
    val = s2; % Modify as needed
end

function val = h2_star(s1)
    val = s1; % Modify as needed
end

function val = g_star(zeta)
    val = zeta; % Modify as needed
end

% Define the whole system
function dydt = system_whole(t, y, params)
    p = y(1);
    d = y(2);
    zeta = y(3);
    s1 = y(4);
    s2 = y(5);
    
    gamma_bar = params(1);
    alpha = params(2);
    k_bar = params(3);
    h1_bar = params(4);
    h2_bar = params(5);
    c1_bar = params(6);
    d1_bar = params(7);
    c2_bar = params(8);
    d2_bar = params(9);
    d_bar = params(10);
    g1_bar = params(11);
    g2_bar = params(12);
    
    dp_dt = p * (1 - (p + d)) - gamma_bar * p;
    dd_dt = k_bar * d * (1 - (p + d)) - gamma_bar * d + alpha * g_star(zeta);
    dzeta_dt = (g1_bar - g2_bar) + h1_bar * h1_star(s1) - h2_bar * h2_star(s2) - d_bar * zeta;
    ds1_dt = c1_bar * p - d1_bar * s1;
    ds2_dt = c2_bar * d - d2_bar * s2;
    
    dydt = [dp_dt; dd_dt; dzeta_dt; ds1_dt; ds2_dt];
end

% Define the QSS system
function dydt = system_qss(t, y, params)
    p = y(1);
    d = y(2);
    zeta = y(3);
    
    gamma_bar = params(1);
    alpha = params(2);
    k_bar = params(3);
    h1_bar = params(4);
    h2_bar = params(5);
    c1_bar = params(6);
    d1_bar = params(7);
    c2_bar = params(8);
    d2_bar = params(9);
    d_bar = params(10);
    g1_bar = params(11);
    g2_bar = params(12);
    dp_dt = p * (1 - (p + d)) - gamma_bar * p;
    dd_dt = k_bar * d * (1 - (p + d)) - gamma_bar * d + alpha * g_star(zeta);
    dzeta_dt = (g1_bar - g2_bar) + h1_bar * c1_bar / d1_bar * p - h2_bar * c2_bar / d2_bar * d - d_bar * zeta;
    
    dydt = [dp_dt; dd_dt; dzeta_dt];
end

% Initial conditions
y0_whole = [1; 1; 0; 0; 0]; % Adjust as necessary
y0_qss = y0_whole(1:3);

% Parameters
gamma_bar = 0;
alpha = 769/26;
k_bar = 30;
g1_bar = 0;
h1_bar = 539.9792;
d_bar = 0;
g2_bar = 0;
h2_bar = 747.6635;
c1_bar = 3.1152;
d1_bar = 9.9688;
c2_bar = 1.8691;
d2_bar = 86.0851;
r = (h1_bar * c1_bar * d2_bar) / (h2_bar * c2_bar * d1_bar);

params_whole = [gamma_bar, alpha, k_bar, h1_bar, h2_bar, c1_bar, d1_bar, c2_bar, d2_bar, d_bar, g1_bar, g2_bar];
params_qss = params_whole;

% Time span
t_span = [0, 40];
t_eval = linspace(t_span(1), t_span(2), 1000);

% Solve the whole system
[t_whole, y_whole] = ode45(@(t, y) system_whole(t, y, params_whole), t_eval, y0_whole);

% Extract solution for the whole system
p_sol_whole = y_whole(:, 1);
d_sol_whole = y_whole(:, 2);
zeta_sol_whole = y_whole(:, 3);
s1_sol = y_whole(:, 4);
s2_sol = y_whole(:, 5);

% Solve the QSS system
[t_qss, y_qss] = ode45(@(t, y) system_qss(t, y, params_qss), t_eval, y0_qss);

% Extract solution for the QSS system
p_sol_qss = y_qss(:, 1);
d_sol_qss = y_qss(:, 2);
zeta_sol_qss = y_qss(:, 3);

% Colorblind-friendly colors
blue       = [0, 114, 178]/255;
vermillion = [213, 94, 0]/255;
green      = [0, 158, 115]/255;
orange     = [230, 159, 0]/255;
purple     = [204, 121, 167]/255;
skyblue    = [86, 180, 233]/255;

% Plot results
figure;

fontSize = 20;
fontSize2 = 42;
lineWidth = 3;

% --- P population ---
subplot(2, 2, 1);
plot(t_whole, p_sol_whole, 'Color', blue, 'LineWidth', lineWidth, 'DisplayName', 'Whole system');
hold on;
%plot(t_qss, p_sol_qss, 'Color', skyblue, 'LineWidth', lineWidth, 'DisplayName', 'QSS');
yline(1 / (1 + r), '--', 'Color', vermillion, 'LineWidth', lineWidth, 'DisplayName', 'Steady state');
legend('FontSize', 14);
xlabel('Time', 'FontSize', fontSize);
ylabel('Values', 'FontSize', fontSize);
title('P population', 'FontSize', fontSize2);


% --- D population ---
subplot(2, 2, 2);
plot(t_whole, d_sol_whole, 'Color', blue, 'LineWidth', lineWidth, 'DisplayName', 'Whole system');
hold on;
%plot(t_qss, d_sol_qss, 'Color', skyblue, 'LineWidth', lineWidth, 'DisplayName', 'QSS');
yline(r / (1 + r), '--', 'Color', vermillion, 'LineWidth', lineWidth, 'DisplayName', 'Steady state');
legend('FontSize', 14);
xlabel('Time', 'FontSize', fontSize);
ylabel('Values', 'FontSize', fontSize);
title('D population', 'FontSize', fontSize2);


% --- Antithetic Control ---
subplot(2, 2, 3);
plot(t_whole, zeta_sol_whole, 'Color', blue, 'LineWidth', lineWidth, 'DisplayName', 'Whole system');
hold on;
%plot(t_qss, zeta_sol_qss, 'Color', skyblue, 'LineWidth', lineWidth, 'DisplayName', 'QSS');
yline(0, '--', 'Color', vermillion, 'LineWidth', lineWidth, 'DisplayName', 'Steady state');
legend('FontSize', 14);
xlabel('Time', 'FontSize', fontSize);
ylabel('Values', 'FontSize', fontSize);
title('Antithetic Control', 'FontSize', fontSize2);


% --- Signalling Molecules ---
subplot(2, 2, 4);
plot(t_whole, s1_sol, 'Color', orange, 'LineWidth', lineWidth, 'DisplayName', 'S1');
hold on;
yline(c1_bar / (d1_bar * (1 + r)), '--', 'Color', orange, 'LineWidth', lineWidth, 'DisplayName', 'Steady state for S1');
plot(t_whole, s2_sol, 'Color', green, 'LineWidth', lineWidth, 'DisplayName', 'S2');
yline(c2_bar * r / (d2_bar * (1 + r)), '--', 'Color', green, 'LineWidth', lineWidth, 'DisplayName', 'Steady state for S2');
legend('FontSize', 14);
xlabel('Time', 'FontSize', fontSize);
ylabel('Values', 'FontSize', fontSize);
title('Signalling Molecules', 'FontSize', fontSize2);


% Plot ratio of populations
ratio_whole = (p_sol_whole) ./ (d_sol_whole + p_sol_whole);
ratio_qss = (p_sol_qss) ./ (d_sol_qss + p_sol_qss);

%figure;
%plot(t_whole, ratio_whole, 'b', 'DisplayName', 'Whole system');
%hold on;
%plot(t_qss, ratio_qss, 'c', 'DisplayName', 'QSS');
%yline(1/(1+r), '--r', 'At steady state');
%legend;
%xlabel('Time');
%ylabel('Ratio');
%title('Ratio of the populations of cells');