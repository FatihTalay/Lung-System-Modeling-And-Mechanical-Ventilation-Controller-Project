% Q1

% Model parameters defined below

% Resistance values
Rtube = 4.0;
Rtb = 0.212;
RbA = 0.082;

% Capacitor values
Ctr = 0.0016;
Cb = 0.013;
CA = 0.15;
Ccw = 0.2;

% KCL will be used to obtain the transfer function

syms v1 v2 v3 s Pmus vdot vout tf_sym

v1_eq = v1/Rtube == vdot;
out_eq = vdot/s == vout;
tf_eq = vout/Pmus == tf_sym;

KCL1 = vdot - s*Ctr*(v3-v1) - (v2-v1)/Rtb == 0;
KCL2 = (v2-v1)/Rtb - s*Cb*(v3-v2) - (v3-v2)/(RbA + 1/(s*CA)) == 0;
KCL3 = s*Ctr*(v3-v1) + s*Cb*(v3-v2) + (v3-v2)/(RbA + 1/(s*CA)) - s*Ccw*(Pmus-v3) == 0;

KCL_all = [KCL1 KCL2 KCL3];

sol_kcl = solve([KCL_all out_eq v1_eq tf_eq], [tf_sym v1 v2 v3 Pmus vdot]);

[tfnum, tfden] = numden(sol_kcl.tf_sym);

tfnum = sym2poly(tfnum);
tfden = sym2poly(tfden);

tfnum = tfnum/tfden(1);
tfden = tfden/tfden(1);

systf = tf(tfnum, tfden);

[z, p, k] = tf2zp(tfnum, tfden);
% pzmap(systf);

tol = sqrt(eps)*10000000;
systf_min = minreal(systf, tol);

R = 1/0.25;
C = 1/(R*2.59);

% Q1c
% subplot(2, 1, 1);
% plot(out.airflow.time, out.airflow.data)
% subplot(2, 1, 2);
% plot(out.tidal_vol.time, out.tidal_vol.data)

% Q1b
% plot(out.out_pmust.time, out.out_pmust.data)

% Q1d
% subplot(2, 1, 1);
% plot(out.tfout.time, out.tfout.data)
% subplot(2, 1, 2);
% plot(out.tidal_vol.time, out.tidal_vol.data)

% Q1f
% subplot(4, 1, 1);
% plot(out.tfminDout.time, out.tfminDout.data)
% subplot(4, 1, 2);
% plot(out.tfminout.time, out.tfminout.data)
% subplot(4, 1, 3);
% plot(out.airflow.time, out.airflow.data)
% subplot(4, 1, 4);
% plot(out.tidal_vol.time, out.tidal_vol.data)