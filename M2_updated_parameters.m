function dYdt = population_model(t, Y, params)
    % Unpack variables
    P = Y(1);
    D = Y(2);
    Z1 = Y(3);
    Z2 = Y(4);
    S1 = Y(5);
    S2 = Y(6);
    
    % Unpack parameters
    k_P = params.k_P;
    k_D = params.k_D;
    K = params.K;
    gamma = params.gamma;
    g1 = params.g1;
    g2 = params.g2;
    k1 = params.k1;
    k2 = params.k2;
    k4 = params.k4;
    d = params.d;
    c1 = params.c1;
    gamma_1 = params.gamma_1;
    c2 = params.c2;
    gamma_2 = params.gamma_2;
    dS1 = params.dS1;
    dS2 = params.dS2;

    % Time-dependent alpha
    alpha_t = 0.8 * exp(-0.033 * t);

    % Define functions
    f = @(Z1) alpha_t * Z1;
    h1 = @(S1) S1;
    h2 = @(S2) S2;

    % Equations
    dPdt = k_P * P * (1 - (P + D)/K) - gamma * P;
    dDdt = k_D * D * (1 - (P + D)/K) - gamma * D + f(Z1);
    dZ1dt = g1 + k1 * h1(S1) - d * Z1 - k4 * Z2 * Z1;
    dZ2dt = g2 + k2 * h2(S2) - d * Z2 - k4 * Z2 * Z1;
    dS1dt = c1 * P - (gamma_1 + dS1) * S1;
    dS2dt = c2 * D - (gamma_2 + dS2) * S2;

    % Output
    dYdt = [dPdt; dDdt; dZ1dt; dZ2dt; dS1dt; dS2dt];
end
% Initial conditions [P0, D0, Z10, Z20, S10, S20]
Y0 = [50000; 50000; 0; 0; 0; 0];

% Parameter values
params.k_P = 0.0144;
params.k_D = 0.0347;
params.K = 200000;
params.gamma = 0;
params.g1 = 0;
params.g2 = 0;
params.k1 = 5.2*10^(-6);
params.k2 = 0.072*10^(-6);
params.k4 = 13.468*10^(-6);
params.d = 49.4*10^(-9);
params.c1 = 0.75*10^(-11);
params.c2 = 450*10^(-11);
params.dS1 = 0.12*10^(-11);
params.dS2 = 0.0495*10^(-11);
params.gamma_1 = 0;
params.gamma_2 = 0.0433;

% Time span
tspan = [0 10000];

% Solve ODE
[t, Y] = ode45(@(t, Y) population_model(t, Y, params), tspan, Y0);

r = (params.k1 * params.c1 * (params.gamma_2+params.dS2)) / (params.k2 * params.c2 * (params.gamma_1+params.dS1));

figure;
sgtitle('Model 2');

subplot(2, 2, 1);
plot(t, Y(:,1), 'b','DisplayName', 'P');
yline(params.K / (1 + r), '--r', 'DisplayName', 'Steady state');
legend;
xlabel('Time');
ylabel('Values');
title('P population');

subplot(2, 2, 2);
plot(t, Y(:,2), 'b', 'DisplayName', 'D');
yline(r*params.K / (1 + r), '--r','DisplayName', 'Steady state');
legend;
xlabel('Time');
ylabel('Values');
title('D population');

subplot(2, 2, 3);
plot(t, Y(:,3), 'b', 'DisplayName', 'Z_1');
hold on;
plot(t, Y(:,4), 'c', 'DisplayName', 'Z_2');
legend;
xlabel('Time');
ylabel('Values');
title('Antithetic Control');

subplot(2, 2, 4);
plot(t, Y(:,5), 'b', 'DisplayName', 'S1');
hold on;
yline(params.c1 / ((params.dS1 + params.gamma_1)* (1 + r)), '--b', 'Steady state for S1');
plot(t, Y(:,6), 'g', 'DisplayName', 'S2');
yline(params.c2 * r / ((params.dS2 + params.gamma_2) * (1 + r)), '--g', 'Steady state for S2');
legend;
xlabel('Time');
ylabel('Values');
title('Signalling molecules');



