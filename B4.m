%% B4: Simulations for omega = 0.1, 1, 10
% Run linear and nonlinear closed-loop simulations
% Use the FBI explicit solution and the full-information regulator
% Generate plots and compute steady-state metrics

clear; close all; clc;

%% Common physical parameters
M = 1; F = 1; Lp = 1; g = 9.81;

alpha = 1;
Amp = 0.5;     % disturbance amplitude
T   = 50;      % disturbance period [s]

tEnd = 100;
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

%% Linearized model used for design
A = [ 0  1  0  0;
      0 -1  0  0;
      0  0  0  1;
      0  1  g  0];

B = [0; 1; 0; -1];
P = [0; 1; 0; -1];            % matched disturbance channel
Pw = [P, zeros(4,2)];         % w = [d1; v1; v2]

% outputs for plots
Cplot = [1 0 0 0;             % s
         0 0 1 0];            % phi

% regulated error definition e = s - alpha*v1
Cr = [1 0 0 0];
Qr = [0 -alpha 0];

%% Disturbance square wave
d1_of_t = @(t) Amp * square_no_toolbox(t, T);

%% Omegas to test
omegas = [0.1, 1, 10];

Results = struct([]);

for idx = 1:numel(omegas)
    omega = omegas(idx);

    fprintf('\n==================== omega = %.3g ====================\n', omega);

    %% 1) Exosystem for w = [d1; v1; v2]
    S = [0  0      0;
         0  0   omega;
         0 -omega  0];

    %% 2) FBI explicit solution Pi and Gamma
    Pi = zeros(4,3);
    Pi(1,2) = alpha;
    Pi(2,3) = alpha*omega;
    Pi(3,2) = -(alpha*omega^2)/(g + omega^2);
    Pi(4,3) = -(alpha*omega^3)/(g + omega^2);

    Gamma = [-1, -alpha*omega^2, alpha*omega];   % 1x3

    % FBI residual check
    res1 = norm(Pi*S - (A*Pi + B*Gamma + Pw), 'fro');
    res2 = norm(Cr*Pi + Qr, 'fro');
    fprintf('FBI residual1 = %.3e, residual2 = %.3e\n', res1, res2);

    %% 3) State feedback K and internal model injection gain Lgain
    if omega >= 5
        desired_poles = [-2, -3, -4, -5];
    else
        desired_poles = [-2, -3, -4, -5];
    end

    K = -place(A, B, desired_poles);   % 1x4
    Lgain = Gamma - K*Pi;              % 1x3

    %% 4) Linear closed-loop simulation with z = [x; v1; v2]
    f_lin = @(t,z) closed_loop_linear(t,z,A,B,P,K,Lgain,omega,d1_of_t);

    z0_lin = [zeros(4,1); 0; 1];
    [tL,zL] = ode45(f_lin, [0 tEnd], z0_lin, opts);

    xL  = zL(:,1:4);
    v1L = zL(:,5);
    v2L = zL(:,6);

    d1L = arrayfun(d1_of_t, tL);
    d2L = alpha*v1L;

    yL  = (Cplot*xL.').';
    sL  = yL(:,1);
    phiL= yL(:,2);
    eL  = sL - d2L;

    Wl = [d1L(:), v1L(:), v2L(:)];
    uL = xL*K.' + Wl*Lgain.';

    %% 5) Nonlinear closed-loop simulation with z = [s sdot phi phidot v1 v2]
    f_nl  = @(t,z) closed_loop_nonlinear(t,z,M,F,Lp,g,K,Lgain,omega,d1_of_t);
    z0_nl = [0;0;0;0; 0;1];

    if omega >= 5
        optsN = odeset(opts, ...
            'Events', @(t,z) stop_if_blowup(t,z), ...
            'MaxStep', 0.01);
        [tN,zN,te,ze,ie] = ode45(f_nl, [0 tEnd], z0_nl, optsN);

        if ~isempty(te)
            fprintf('Nonlinear stopped early at t=%.3f, event=%d\n', te(end), ie(end));
        end
    else
        [tN,zN] = ode45(f_nl, [0 tEnd], z0_nl, opts);
        te = []; ze = []; ie = [];
    end

    sN     = zN(:,1);
    sdotN  = zN(:,2);
    phiN   = zN(:,3);
    phidotN= zN(:,4);
    v1N    = zN(:,5);
    v2N    = zN(:,6);

    d1N = arrayfun(d1_of_t, tN);
    d2N = alpha*v1N;
    eN  = sN - d2N;

    Xn = [sN, sdotN, phiN, phidotN];
    Wn = [d1N(:), v1N(:), v2N(:)];
    muN = Xn*K.' + Wn*Lgain.';

    %% 6) Metrics over the final part of the simulation
    tSS = 0.75*tEnd;

    idxLSS = tL > tSS;
    idxNSS = tN > tSS;

    met.lin.maxE   = max(abs(eL(idxLSS)));
    met.lin.rmsE   = rms(eL(idxLSS));
    met.lin.rmsU   = rms(uL(idxLSS));
    met.lin.maxPhi = max(abs(phiL(idxLSS)));

    if any(idxNSS)
        met.nl.maxE   = max(abs(eN(idxNSS)));
        met.nl.rmsE   = rms(eN(idxNSS));
        met.nl.rmsU   = rms(muN(idxNSS));
    else
        met.nl.maxE   = max(abs(eN));
        met.nl.rmsE   = rms(eN);
        met.nl.rmsU   = rms(muN);
    end
    met.nl.maxPhi = max(abs(phiN)); % full horizon

    Results(idx).omega   = omega;
    Results(idx).K       = K;
    Results(idx).Lgain   = Lgain;
    Results(idx).metrics = met;
    Results(idx).event_t = te;
    Results(idx).event_i = ie;

    fprintf('Linear:    max|e|=%.3e, rms|e|=%.3e, rms u=%.3e, max|phi|=%.3e\n', ...
        met.lin.maxE, met.lin.rmsE, met.lin.rmsU, met.lin.maxPhi);
    fprintf('Nonlinear: max|e|=%.3e, rms|e|=%.3e, rms mu=%.3e, max|phi|=%.3e\n', ...
        met.nl.maxE, met.nl.rmsE, met.nl.rmsU, met.nl.maxPhi);

    %% 7) Plots for each omega
    figure('Color','w','Name',sprintf('B4 omega=%.3g (Linear)',omega));
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    nexttile;
    plot(tL, sL, 'LineWidth',1.5); hold on;
    plot(tL, d2L, '--', 'LineWidth',1.2);
    plot(tL, eL, ':', 'LineWidth',1.2);
    grid on; xlabel('t [s]'); ylabel('s, d_2, e');
    legend('s(t)','d_2(t)=\alpha v_1','e(t)=s-d_2','Location','best');
    title(sprintf('Linear: \\omega=%.3g',omega));

    nexttile;
    plot(tL, phiL, 'LineWidth',1.5);
    grid on; xlabel('t [s]'); ylabel('\phi(t) [rad]');
    title('Linear: \phi(t)');

    nexttile;
    plot(tL, uL, 'LineWidth',1.5); hold on;
    plot(tL, d1L, '--', 'LineWidth',1.0);
    grid on; xlabel('t [s]'); ylabel('u(t), d_1(t)');
    legend('u(t)','d_1(t)','Location','best');
    title('Linear: u(t) and d_1(t)');

    figure('Color','w','Name',sprintf('B4 omega=%.3g (Nonlinear)',omega));
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    nexttile;
    plot(tN, sN, 'LineWidth',1.5); hold on;
    plot(tN, d2N, '--', 'LineWidth',1.2);
    plot(tN, eN, ':', 'LineWidth',1.2);
    grid on; xlabel('t [s]'); ylabel('s, d_2, e');
    legend('s(t)','d_2(t)=\alpha v_1','e(t)=s-d_2','Location','best');
    title(sprintf('Nonlinear: \\omega=%.3g',omega));

    nexttile;
    plot(tN, phiN, 'LineWidth',1.5);
    grid on; xlabel('t [s]'); ylabel('\phi(t) [rad]');
    title('Nonlinear: \phi(t)');

    nexttile;
    plot(tN, muN, 'LineWidth',1.5); hold on;
    plot(tN, d1N, '--', 'LineWidth',1.0);
    grid on; xlabel('t [s]'); ylabel('\mu(t), d_1(t)');
    legend('\mu(t)','d_1(t)','Location','best');
    title('Nonlinear: \mu(t) and d_1(t)');
end

%% Quick discussion guide
fprintf('\n==================== B4 DISCUSSION GUIDE ====================\n');
fprintf(['As omega increases, the sinusoidal reference d_2(t)=alpha*v_1(t) becomes faster.\n' ...
         'The internal-model/feedforward targets that frequency, but higher omega typically\n' ...
         'requires larger control effort (terms scale with omega and omega^2).\n' ...
         'The linear model may still show small steady-state error, while the nonlinear\n' ...
         'model can deviate due to strong coupling via sin(phi), cos(phi); for large omega,\n' ...
         'phi(t) can grow and the simulation may trigger a safety event (large-angle regime).\n']);

%% ================= helper functions =================
function dz = closed_loop_linear(t,z,A,B,P,K,Lgain,omega,d1_of_t)
    x  = z(1:4);
    v1 = z(5);
    v2 = z(6);

    d1 = d1_of_t(t);
    w  = [d1; v1; v2];

    u  = K*x + Lgain*w;

    dx  = A*x + B*u + P*d1;
    dv1 =  omega*v2;
    dv2 = -omega*v1;

    dz = [dx; dv1; dv2];
end

function dz = closed_loop_nonlinear(t,z,M,F,Lp,g,K,Lgain,omega,d1_of_t)
    s      = z(1);
    sdot   = z(2);
    phi    = z(3);
    phidot = z(4);
    v1     = z(5);
    v2     = z(6);

    d1 = d1_of_t(t);
    w  = [d1; v1; v2];

    xlin = [s; sdot; phi; phidot];
    mu = K*xlin + Lgain*w;

    % nonlinear plant
    sddot   = (mu + d1 - F*sdot)/M;
    phiddot = (g/Lp)*sin(phi) - (1/Lp)*sddot*cos(phi);

    dv1 =  omega*v2;
    dv2 = -omega*v1;

    dz = [sdot; sddot; phidot; phiddot; dv1; dv2];
end

function y = square_no_toolbox(t, T)
    % +1 or -1 square wave with period T
    k = floor(t/(T/2));
    y = 1 - 2*mod(k,2); % +1,-1,+1,-1,...
end

function [value,isterminal,direction] = stop_if_blowup(~,z)
    % Stop the simulation if limits are exceeded
    phi = z(3);
    s   = z(1);

    v1 = 0.8*pi - abs(phi);
    v2 = 50000     - abs(s);
    v3 = double(all(isfinite(z))) - 0.5;

    value      = [v1; v2; v3];
    isterminal = [1;  1;  1];
    direction  = [-1; -1; -1];
end
