%% B3: Apply B1 full-information regulator to the nonlinear inverted pendulum model
clear; close all; clc;

%% Parameters
M = 1;
F = 1;
L = 1;
g = 9.81;

alpha = 1;
omega = 0.1;

Amp = 0.5;
T   = 50;

tEnd = 200;
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

%% Linearized model used for design
A = [ 0  1  0  0;
      0 -1  0  0;
      0  0  0  1;
      0  1  g  0];

B = [0; 1; 0; -1];

P = [0; 1; 0; -1];

S = [0  0      0;
     0  0   omega;
     0 -omega  0];

Pw = [P, zeros(4,2)];

Cr = [1 0 0 0];
Qr = [0 -alpha 0];

%% FBI explicit solution Pi and Gamma
Pi = zeros(4,3);
Pi(1,2) = alpha;
Pi(2,3) = alpha*omega;
Pi(3,2) = -(alpha*omega^2)/(g + omega^2);
Pi(4,3) = -(alpha*omega^3)/(g + omega^2);

Gamma = [-1, -alpha*omega^2, alpha*omega];

% FBI residual check
res1 = norm(Pi*S - (A*Pi + B*Gamma + Pw), 'fro');
res2 = norm(Cr*Pi + Qr, 'fro');
fprintf('FBI residual1 = %.3e, residual2 = %.3e\n', res1, res2);

%% Controller gains K and L
desired_poles = [-2, -3, -4, -5];
K = -place(A, B, desired_poles);
Lgain = Gamma - K*Pi;

fprintf('K (used in mu = Kx + Lw) = [%s]\n', num2str(K, ' % .5f'));
fprintf('L = [%s]\n', num2str(Lgain, ' % .5f'));

%% Disturbance d1(t) square wave
d1_of_t = @(t) Amp * square_no_toolbox(t, T);

%% Nonlinear closed-loop simulation
f_nl = @(t,z) closed_loop_nonlinear(t,z,M,F,L,g,K,Lgain,omega,d1_of_t);

z0 = [0;0;0;0; 0;1];
[t,z] = ode15s(f_nl, [0 tEnd], z0, opts);

s     = z(:,1);
sdot  = z(:,2);
phi   = z(:,3);
phidot= z(:,4);
v1    = z(:,5);
v2    = z(:,6);

d1 = arrayfun(d1_of_t, t);
d2 = alpha*v1;
e  = s - d2;

%% Control input
u = zeros(size(t));
for kidx = 1:length(t)
    xlin = [s(kidx); sdot(kidx); phi(kidx); phidot(kidx)];
    w = [d1(kidx); v1(kidx); v2(kidx)];
    u(kidx) = K*xlin + Lgain*w;
end

%% Plots
figure('Color','w','Name','B3 Nonlinear: y(t), u(t)');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(t, s, 'LineWidth',1.5); hold on;
plot(t, d2, '--', 'LineWidth',1.2);
plot(t, e, ':', 'LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('s(t)');
legend('s(t)','d_2(t)=\alpha v_1(t)','e(t)=s-d_2','Location','best');
title('B3 (Nonlinear): tracking and regulation of s(t)');

nexttile;
plot(t, phi, 'LineWidth',1.5);
grid on; xlabel('t [s]'); ylabel('\phi(t) [rad]');
title('B3 (Nonlinear): pendulum angle \phi(t)');

nexttile;
plot(t, u, 'LineWidth',1.5); hold on;
plot(t, d1, '--', 'LineWidth',1.0);
grid on; xlabel('t [s]'); ylabel('\mu(t)');
legend('\mu(t)','d_1(t)','Location','best');
title('B3 (Nonlinear): control input \mu(t) and disturbance d_1(t)');

%% Regulation metrics
idx_ss = t > 150;
fprintf('Nonlinear: steady max |e| (t>150s) = %.3e\n', max(abs(e(idx_ss))));
fprintf('Nonlinear: steady RMS |e| (t>150s) = %.3e\n', rms(e(idx_ss)));

%% ================= helper functions =================
function dz = closed_loop_nonlinear(t,z,M,F,L,g,K,Lgain,omega,d1_of_t)
    s     = z(1);
    sdot  = z(2);
    phi   = z(3);
    phidot= z(4);
    v1    = z(5);
    v2    = z(6);

    d1 = d1_of_t(t);
    w  = [d1; v1; v2];

    xlin = [s; sdot; phi; phidot];
    u = K*xlin + Lgain*w;

    sddot = (u + d1 - F*sdot)/M;
    phiddot = (g/L)*sin(phi) - (1/L)*sddot*cos(phi);

    dv1 =  omega*v2;
    dv2 = -omega*v1;

    dz = [sdot; sddot; phidot; phiddot; dv1; dv2];
end

function y = square_no_toolbox(t, T)
    k = floor(t/(T/2));
    y = 1 - 2*mod(k,2);
end
