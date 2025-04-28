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

function val = g_star(z1)
    val = z1; % Modify as needed
end



% Define ODE systems
function dydt = system_whole(t, y, params)
    p = max(y(1),0);
    d = max(y(2),0);
    z1 = max(y(3),0);
    z2 = max(y(4),0);
    s1 = max(y(5),0);
    s2 = max(y(6),0);
    s3 = max(y(7),0);
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
    k4 = params(12);
    d_bar = params(13);
    g1_bar = params(14);
    g2_bar = params(15);

    dp_dt = p * (1 - (p + d)) - gamma_bar * p + k3_bar * f_star(s3);
    dd_dt = k_bar * d * (1 - (p + d)) - gamma_bar * d;
    dz1_dt = g1_bar + h1_bar * s2 - d_bar * z1 - k4 * z1 * z2;
    dz2_dt = g2_bar + h2_bar * s1 - d_bar * z2 - k4 * z1 * z2;
    ds1_dt = c1_bar * p - d1_bar * s1;
    ds2_dt = c2_bar * d - d2_bar * s2;
    ds3_dt = g_bar * g_star(z1) - d3_bar * s3;
    
    dydt = [dp_dt; dd_dt; dz1_dt;dz2_dt; ds1_dt; ds2_dt; ds3_dt];
end

function dydt = system_qss(t, y, params)
    p = max(y(1),0);
    d = max(y(2),0);
    z1 = max(y(3),0);
    z2 = max(y(4),0);

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
    k4 = params(12);
    d_bar = params(13);

    dp_dt = p * (1 - (p + d)) - gamma_bar * p + k3_bar * g_bar * z1 / d3_bar;
    dd_dt = k_bar * d * (1 - (p + d)) - gamma_bar * d;
    dz1_dt = + h1_bar * c2_bar / d2_bar * d - d_bar * z1 - k4 * z1 * z2;
    dz2_dt = + h2_bar * c1_bar / d1_bar * p - d_bar * z2 - k4 * z1 * z2;
   
    dydt = [dp_dt; dd_dt; dz1_dt; dz2_dt];
end

% Initial conditions
y0_whole = [0.5; 0.2; 0; 0; 0;0;0];
y0_qss = y0_whole(1:4);

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
k4 = 0;
r = c1_bar * d2_bar * h2_bar / (c2_bar * d1_bar * h1_bar);
params_whole = [gamma_bar, k3_bar, k_bar, h1_bar, h2_bar, c1_bar, d1_bar, c2_bar, d2_bar, g_bar, k4,d_bar, d3_bar, g1_bar, g2_bar];
params_qss = params_whole(1:13);

% Time span
tspan = [0 400];

% Solve the systems
[t_whole, y_whole] = ode45(@(t,y) system_whole(t, y, params_whole), tspan, y0_whole);
[t_qss, y_qss] = ode45(@(t,y) system_qss(t, y, params_qss), tspan, y0_qss);

% Extract solutions
p_sol_whole = y_whole(:,1);
d_sol_whole = y_whole(:,2);
z1_sol_whole = y_whole(:,3);
z2_sol_whole = y_whole(:,4);
s1_sol = y_whole(:,5);
s2_sol = y_whole(:,6);
s3_sol = y_whole(:,7);

p_sol_qss = y_qss(:,1);
d_sol_qss = y_qss(:,2);
z1_sol_qss = y_qss(:,3);
z2_sol_qss = y_qss(:,4);


% Plot results
figure;
tiledlayout(2,2);

% P population
nexttile;
plot(t_whole, p_sol_whole, 'b');
%plot(t_qss, p_sol_qss, 'c');
hold on;
yline(1 / (1+r), '--r');
legend('Whole system', 'Steady state');
xlabel('Time'); ylabel('Values'); title('P population');

% D population
nexttile;
plot(t_whole, d_sol_whole, 'b');
%plot(t_qss, d_sol_qss, 'c');
hold on;
yline(r / (1+r), '--r');
legend('Whole system', 'Steady state');
xlabel('Time'); ylabel('Values'); title('D population');

% Antithetic Control
nexttile;
plot(t_whole, z1_sol_whole, 'b');
%plot(t_qss, z1_sol_qss, 'c');
hold on;
yline(0, '--r');
legend('Whole system', 'Steady state');
xlabel('Time'); ylabel('Values'); title('Z1');

nexttile;
plot(t_whole, z2_sol_whole, 'b');
%plot(t_qss, z2_sol_qss, 'c');
hold on;
yline(0, '--r');
legend('Whole system', 'Steady state');
xlabel('Time'); ylabel('Values'); title('Z2');


% Signalling molecules
nexttile;
plot(t_whole, s1_sol, 'y', t_whole, s2_sol, 'g', t_whole, s3_sol, 'b');
yline(c2_bar * r / (d2_bar * (1+r)), '--g', 'Steady State');
yline(c1_bar  / (d1_bar * (1+r)), '--y', 'Steady State');
hold on;
legend('S1', 'S2', 'S3');
xlabel('Time'); ylabel('Values'); title('Signalling Molecules');

% Plot ratio of populations
ratio_whole = (p_sol_whole) ./ (d_sol_whole + p_sol_whole);
ratio_qss = (p_sol_qss) ./ (d_sol_qss + p_sol_qss);

figure;
plot(t_whole, ratio_whole, 'b', 'DisplayName', 'Whole system');
hold on;
plot(t_qss, ratio_qss, 'c', 'DisplayName', 'QSS');
yline(1/(1+r), '--r', 'At steady state');
legend;
xlabel('Time');
ylabel('Ratio');
title('Ratio of the populations of cells');