%% B1: Full-information regulator design
% State: x = [s; sdot; phi; phidot], input: u = mu
% Disturbance: d1(t) enters via P
% Reference: d2(t) = alpha*v1(t), error: e(t) = s(t) - d2(t)
% Exosystem: w = [d1; v1; v2], v1dot = omega*v2, v2dot = -omega*v1
% FBI equations: Pi*S = A*Pi + B*Gamma + Pw, and Cr*Pi + Qr = 0
% Control law: u(t) = K*x(t) + L*w(t), where L = Gamma - K*Pi
% K is chosen by pole placement to stabilize A + B*K

clear; clc;

%% 0) Parameters
g     = 9.81;
alpha = 1;
omega = 0.1;

% Disturbance settings
Amp = 0.5;
T   = 50;

%% 1) Linearized plant matrices
A = [ 0  1  0  0;
      0 -1  0  0;
      0  0  0  1;
      0  1  g  0];

B = [0; 1; 0; -1];

% Disturbance channel
P = [0; 1; 0; -1];

% Exosystem matrix for w = [d1; v1; v2]
S = [0  0      0;
     0  0   omega;
     0 -omega  0];

% Pw = [P, 0, 0]
Pw = [P, zeros(4,2)];

% Error mapping e = Cr*x + Qr*w
Cr = [1 0 0 0];
Qr = [0 -alpha 0];

%% 2) FBI solution Pi and Gamma
Pi = zeros(4,3);
Pi(1,2) = alpha;
Pi(2,3) = alpha*omega;
Pi(3,2) = -(alpha*omega^2)/(g + omega^2);
Pi(4,3) = -(alpha*omega^3)/(g + omega^2);

Gamma = [-1, -alpha*omega^2, alpha*omega];

%% 3) FBI residual check
res1 = norm(Pi*S - (A*Pi + B*Gamma + Pw), 'fro');
res2 = norm(Cr*Pi + Qr, 'fro');

fprintf('FBI residual ||Pi*S - (A*Pi+B*Gamma+Pw)||_F = %.3e\n', res1);
fprintf('FBI residual ||Cr*Pi + Qr||_F               = %.3e\n', res2);

%% 4) Pole placement for K and compute L
desired_poles = [-6, -2, -3, -5];

K = -place(A, B, desired_poles);

L = Gamma - K*Pi;

%% 5) Display controller gains
disp('--- Full-information regulator (B1) ---');
disp('State x = [s; sdot; phi; phidot], exosystem w = [d1; v1; v2]');
disp('Control law:  u(t) = K*x(t) + L*w(t)');
fprintf('K = [% .6f  % .6f  % .6f  % .6f]\n', K(1),K(2),K(3),K(4));
fprintf('L = [% .6f  % .6f  % .6f]\n', L(1),L(2),L(3));

%% Disturbance function handle
d1_of_t = @(t) Amp * sign(sin(2*pi*t/50));
