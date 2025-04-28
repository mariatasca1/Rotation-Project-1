% Define functions
function val = f_star(s3)
    val = s3; % Modify as needed
end

function val = h1_star(s2)
    val = s2; % Modify as needed
end

function val = h2_star(s1)
    val = s1; % Modify as needed
end

function val = g_star(zeta)
    val = zeta; % Modify as needed
end

% Define ODE systems
function dydt = system_whole(t, y, params)
    p = y(1);
    d = y(2);
    zeta = y(3);
    s1 = y(4);
    s2 = y(5);
    s3 = y(6);
    gamma_bar = params(1);
    k3_bar = params(2);
    k_bar = params(3);
    h1_bar = params(4);
    h2_bar = params(5);
    c1_bar = params(6);
    d1_bar = params(7);
    c2_bar = params(8);
    d2_bar = params(9);
    g_bar = params(10);
    d3_bar = params(11);
    g1_bar = params(12);
    g2_bar = params(13);
    d_bar = params(14);
   
    dp_dt = p * (1 - (p + d)) - gamma_bar * p + k3_bar * f_star(s3);
    dd_dt = k_bar * d * (1 - (p + d)) - gamma_bar * d;
    dzeta_dt = (g1_bar - g2_bar) + h1_bar * h1_star(s2) - h2_bar * h2_star(s1) - d_bar * zeta;
    ds1_dt = c1_bar * p - d1_bar * s1;
    ds2_dt = c2_bar * d - d2_bar * s2;
    ds3_dt = g_bar * g_star(zeta) - d3_bar * s3;
    
    dydt = [dp_dt; dd_dt; dzeta_dt; ds1_dt; ds2_dt; ds3_dt];
end

function dydt = system_qss(t, y, params)
    p = y(1);
    d = y(2);
    zeta = y(3);
    gamma_bar = params(1);
    k3_bar = params(2);
    k_bar = params(3);
    h1_bar = params(4);
    h2_bar = params(5);
    c1_bar = params(6);
    d1_bar = params(7);
    c2_bar = params(8);
    d2_bar = params(9);
    g_bar = params(10);
    d3_bar = params(11);
   
    dp_dt = p * (1 - (p + d)) - gamma_bar * p + k3_bar * g_bar * zeta / d3_bar;
    dd_dt = k_bar * d * (1 - (p + d)) - gamma_bar * d;
    dzeta_dt = h1_bar * c2_bar / d2_bar * d - h2_bar * c1_bar / d1_bar * p;
   
    dydt = [dp_dt; dd_dt; dzeta_dt];
end

% Initial conditions
y0_whole = [1; 1; 0; 0; 0;0];
y0_qss = y0_whole(1:3);

% Parameters
gamma_bar = 0;
k3_bar = 107/167;
k_bar = 38/51;
g1_bar = 0;
h1_bar = 270;
d_bar = 0;
g2_bar = 0;
h2_bar = 40;
c1_bar = 34;
d1_bar = 78;
c2_bar = 9;
d2_bar = 59;
g_bar = 12;
d3_bar = 44;

r = c1_bar * d2_bar * h2_bar / (c2_bar * d1_bar * h1_bar);
params_whole = [gamma_bar, k3_bar, k_bar, h1_bar, h2_bar, c1_bar, d1_bar, c2_bar, d2_bar, g_bar, d3_bar, g1_bar, g2_bar, d_bar];
params_qss = params_whole(1:11);

% Time span
tspan = [0 400];

% Solve the systems
[t_whole, y_whole] = ode45(@(t,y) system_whole(t, y, params_whole), tspan, y0_whole);
[t_qss, y_qss] = ode45(@(t,y) system_qss(t, y, params_qss), tspan, y0_qss);

% Extract solutions
p_sol_whole = y_whole(:,1);
d_sol_whole = y_whole(:,2);
zeta_sol_whole = y_whole(:,3);
s1_sol = y_whole(:,4);
s2_sol = y_whole(:,5);
s3_sol = y_whole(:,6);

p_sol_qss = y_qss(:,1);
d_sol_qss = y_qss(:,2);
zeta_sol_qss = y_qss(:,3);


% Colorblind-friendly colors
blue       = [0, 114, 178]/255;
vermillion = [213, 94, 0]/255;
green      = [0, 158, 115]/255;
orange     = [230, 159, 0]/255;
purple     = [204, 121, 167]/255;
skyblue    = [86, 180, 233]/255;

% Plot settings
figure;
tiledlayout(2,2);
fontSize = 20;
fontSize2 = 42;
lineWidth = 1.5;

% P population
nexttile;
plot(t_whole, p_sol_whole, 'Color', blue, 'LineWidth', lineWidth);
hold on;
plot(t_qss, p_sol_qss, 'Color', skyblue, 'LineWidth', lineWidth);
yline(1 / (1+r), '--', 'Color', vermillion, 'LineWidth', lineWidth);
legend('Whole system', 'QSS', 'Steady state', 'FontSize', 14);
xlabel('Time', 'FontSize', fontSize);
ylabel('Values', 'FontSize', fontSize);
title('P population', 'FontSize', 42);


% D population
nexttile;
plot(t_whole, d_sol_whole, 'Color', blue, 'LineWidth', lineWidth);
hold on;
plot(t_qss, d_sol_qss, 'Color', skyblue, 'LineWidth', lineWidth);
yline(r / (1+r), '--', 'Color', vermillion, 'LineWidth', lineWidth);
legend('Whole system', 'QSS', 'Steady state', 'FontSize', 14);
xlabel('Time', 'FontSize', fontSize);
ylabel('Values', 'FontSize', fontSize);
title('D population', 'FontSize', 42);


% Antithetic Control
nexttile;
plot(t_whole, zeta_sol_whole, 'Color', blue, 'LineWidth', lineWidth);
hold on;
plot(t_qss, zeta_sol_qss, 'Color', skyblue, 'LineWidth', lineWidth);
yline(0, '--', 'Color', vermillion, 'LineWidth', lineWidth);
legend('Whole system', 'QSS', 'Steady State', 'FontSize', 14);
xlabel('Time', 'FontSize', fontSize);
ylabel('Values', 'FontSize', fontSize);
title('Antithetic Control', 'FontSize', fontSize2);

% Signalling Molecules
nexttile;
plot(t_whole, s1_sol, 'Color', orange, 'LineWidth', lineWidth);
hold on;
plot(t_whole, s2_sol, 'Color', green, 'LineWidth', lineWidth);
plot(t_whole, s3_sol, 'Color', purple, 'LineWidth', lineWidth);
yline(c2_bar * r / (d2_bar * (1+r)), '--', 'Color', green, ...
    'LineWidth', lineWidth, 'DisplayName', 'Steady State S2');
yline(c1_bar / (d1_bar * (1+r)), '--', 'Color', orange, ...
    'LineWidth', lineWidth, 'DisplayName', 'Steady State S1');
legend('S1', 'S2', 'S3','Steady State', 'FontSize', 14);
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