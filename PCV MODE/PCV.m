syms R C Tao_f K K_i w s;
assume(R > 0 & C > 0 & Tao_f > 0);

% System Transfer Function
G_f = 1 / (Tao_f*s + 1);
G = (1/R) / (s + 1/(R*C));

% Getting Characteristic polynomial
F_pi = K_i/s;
L = F_pi * G_f * G;
T_pi = simplify(L / (1 + L));
[~, Pcs] = numden(T_pi);
coefficent = coeffs(Pcs, s, "All");
disp("Characteristic polynomial for Close Loop of PI Controller");
fprintf("(%s)*s^2 + (%s)*s + (%s)\n\n", string(coefficent(1)), string(coefficent(2)), string(coefficent(3)));

F_p = K;
L_2 = F_p * G_f * G;
T_p = simplify(L_2 / (1 + L_2));
[~, Pcs_2] = numden(T_p);
coefficent_2 = coeffs(Pcs_2, s, "All");
disp("Characteristic polynomial for Close Loop of P Controller");
fprintf("(%s)*s^2 + (%s)*s + (%s)\n\n", string(coefficent_2(1)), string(coefficent_2(2)), string(coefficent_2(3)));
Delta = coefficent_2(2)^2 - 4 * coefficent_2(1) * coefficent_2(3);
sol = solve(Delta==0, K);
% Display the solutions for K and their conditions
disp("Solutions for K when Delta = 0:");
disp(sol);



% Script to plot PCV Mode simulation results
figure('Name', 'PCV Mode Waveforms', 'Color', 'w');

% Plot Filtered Pressure
subplot(3,1,1);
plot(out.P_vent.Time, out.P_vent.Data);
title('Filtered Airway Pressure');
ylabel('Pressure [cmH2O]');
grid on;

% Plot Tidal Volume
subplot(3,1,2);
plot(out.Tidal_Volume.Time, out.Tidal_Volume.Data);
title('Tidal Volume (V)');
ylabel('Volume [L]');
grid on;

% Plot Airflow
subplot(3,1,3);
plot(out.Airflow.Time, out.Airflow.Data);
title('Airflow (V'')');
ylabel('Flow [L/s]');
xlabel('Time [s]');
grid on;







