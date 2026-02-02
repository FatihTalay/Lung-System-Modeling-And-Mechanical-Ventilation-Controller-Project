% Parameter definitions
R = 4;         % Resistance
C = 0.0965;    % Compliance ( C * 1.5 for COPD)

%% 

tau_f = 0.1;


K_pi = (16 * R) / (25 * tau_f);
K_i  = (64 * R) / (25 * tau_f^2);
K_pd = (32 * R) / (5 * tau_f) - (1 / C);
K_d  = (23 * R) / 5 - (tau_f / C);


s = tf('s'); 

Gs = (1/R) / ( s + 1/(R*C));
F_PI = K_pi + K_i / s;
F_PD = K_pd + K_d * s;

Ts = minreal(F_PI * Gs / ( 1 + (F_PI + F_PD) * Gs ));