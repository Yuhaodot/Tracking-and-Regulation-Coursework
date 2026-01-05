%% B2: Linear simulation and plots
% Use the full-information controller from B1
% Simulate the linear plant with square-wave disturbance d1
% Plot s, phi, and control input u, and show that e = s - d2 converges to zero

clear; close all; clc;

%% 0. Parameters
g     = 9.81;
alpha = 1;
omega = 0.1;

Amp = 0.5;      % d1 square wave amplitude
T   = 50;       % d1 square wave period [s]

tEnd = 200;     % simulate horizon
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

%% 1. Linearized plant
A = [ 0  1  0  0;
    0 -1  0  0;
    0  0  0  1;
    0  1  g  0];

B = [0; 1; 0; -1];

% matched disturbance channel
P = [0; 1; 0; -1];

% output for plots: y = [s; phi]
C = [1 0 0 0;
    0 0 1 0];

%% 2. Exosystem for w = [d1; v1; v2]
S = [0  0      0;
    0  0   omega;
    0 -omega  0];

% Pw = [P, 0, 0]
Pw = [P, zeros(4,2)];

% regulated error e = s - alpha*v1
Cr = [1 0 0 0];
Qr = [0 -alpha 0];

%% 3. FBI explicit solution Pi and Gamma
Pi = zeros(4,3);
Pi(1,2) = alpha;
Pi(2,3) = alpha*omega;
Pi(3,2) = -(alpha*omega^2)/(g + omega^2);
Pi(4,3) = -(alpha*omega^3)/(g + omega^2);

Gamma = [-1, -alpha*omega^2, alpha*omega]; % 1x3

% FBI residual check
res1 = norm(Pi*S - (A*Pi + B*Gamma + Pw), 'fro');
res2 = norm(Cr*Pi + Qr, 'fro');
fprintf('FBI residual1 = %.3e, residual2 = %.3e\n', res1, res2);

%% 4. Controller gains K and L
desired_poles = [-2, -3, -4, -5];
K = -place(A, B, desired_poles); % pole placement using place
L = Gamma - K*Pi;

fprintf('K = [%s]\n', num2str(K, ' % .5f'));%print
fprintf('L = [%s]\n', num2str(L, ' % .5f'));

%% 5. Closed-loop simulation
d1_of_t = @(t) Amp * sign(sin(2*pi*t/T));   % square wave, no toolbox;

% Augmented state z = [x; v1; v2]
f_cl = @(t,z) closed_loop_linear(t,z,A,B,P,K,L,omega,d1_of_t);


% initial condition: x0 = 0, v1_0 = 0, v2_0 = 1 so v1 equals sin omega t
z0 = [zeros(4,1); 0; 1];

[t,z] = ode45(f_cl, [0 tEnd], z0, opts); % linear slover: LK, option:ode23, ode15s, ode15;
% ode45 by defauly

x  = z(:,1:4);
v1 = z(:,5);
v2 = z(:,6);

% Signals
d1 = d1_of_t(t);
d2 = alpha*v1;                % reference: sin
y  = (C*x.').';               % y=[s,phi]
s  = y(:,1);

phi= y(:,2);

% Control input u(t) = Kx + Lw
u = zeros(size(t)); %pre-allowcation
for kidx = 1:length(t)
    w = [d1(kidx); v1(kidx); v2(kidx)];
    u(kidx) = K*x(kidx,:).' + L*w;
end

% Regulation error e = s - d2
e = s-d2;

%% 6. Plot
figure('Color','w','Name','B2 Linear: y(t), u(t)');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(t, s, 'LineWidth',1.5); hold on;
plot(t, d2, '--', 'LineWidth',1.2);
plot(t, e, ':', 'LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('s(t)');
legend('s(t)','d_2(t)=\alpha v_1(t)','e(t)=s-d_2','Location','best');
title('B2 (Linear): tracking/regulation of s(t)');

nexttile;
plot(t, phi, 'LineWidth',1.5);
grid on; xlabel('t [s]'); ylabel('\phi(t) [rad]');
title('B2 (Linear): pendulum angle \phi(t)');

nexttile;
plot(t, u, 'LineWidth',1.5); hold on;
plot(t, d1, '--', 'LineWidth',1.0);
grid on; xlabel('t [s]'); ylabel('u(t)');
legend('u(t)','d_1(t)','Location','best');
title('B2 (Linear): control input u(t) and disturbance d_1(t)');

%% 7. Regulation metrics
idx_ss = t > 150;
fprintf('Steady-state |e(t)| max (t>150s) = %.3e\n', max(abs(e(idx_ss))));
fprintf(['Steady-state |e(t)| RMSE(root mean square error. ' ...
    'Please note that MSE=mean squre err, RMSE =root MSE ) (t>150s) = %.3e\n'], rms(e(idx_ss)));

%% ==================== helper function ====================
function dz = closed_loop_linear(t,z,A,B,P,K,L,omega,d1_of_t)
x  = z(1:4);
v1 = z(5);
v2 = z(6);

d1 = d1_of_t(t);
w  = [d1; v1; v2];

u  = K*x + L*w;

dx  = A*x + B*u + P*d1;
dv1 =  omega*v2;   % v1 is sin omega t, v2 is cos omega t
dv2 = -omega*v1;

dz = [dx; dv1; dv2];
end
