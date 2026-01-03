########B2

%% 0) Given parameters (B1/B2)
g     = 9.81;
alpha = 1;
omega = 0.1;

Amp = 0.5;      % d1 square wave amplitude
T   = 50;       % d1 square wave period [s]

tEnd = 200;     % simulate horizon
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);



FBI residual1 = 1.880e-19, residual2 = 0.000e+00
K = [12.23242  16.69827  93.04242  29.69827]
L = [-1.00000 -12.14767  -1.56680]
Steady-state |e(t)| max (t>150s) = 1.221e-15
Steady-state |e(t)| RMSE(root mean square error. Please note that MSE=mean squre err, RMSE =root MSE ) (t>150s) = 2.432e-16
>> 




########B3


FBI residual1 = 1.880e-19, residual2 = 0.000e+00
K (used in mu = Kx + Lw) = [12.23242  16.69827  93.04242  29.69827]
L = [-1.00000 -12.14767  -1.56680]
Nonlinear: steady max |e| (t>150s) = 2.667e-09
Nonlinear: steady RMS |e| (t>150s) = 1.282e-09








########B4w10_reB2

parameters

%% 0) Given parameters (B1/B2)
g     = 9.81;
alpha = 1;
omega = 10;

Amp = 0.5;      % d1 square wave amplitude
T   = 50;       % d1 square wave period [s]

tEnd = 10;     % simulate horizon
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);




FBI residual1 = 0.000e+00, residual2 = 0.000e+00
K = [12.23242  16.69827  93.04242  29.69827]
L = [-1.00000  -27.50205  113.46872]
Steady-state |e(t)| max (t>150s) = 
Steady-state |e(t)| RMSE(root mean square error. Please note that MSE=mean squre err, RMSE =root MSE ) (t>150s) = NaN
>> 







########B4w1_reB2


g     = 9.81;
alpha = 1;
omega = 1;

Amp = 0.5;      % d1 square wave amplitude
T   = 50;       % d1 square wave period [s]

tEnd = 100;     % simulate horizon
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);


FBI residual1 = 4.163e-17, residual2 = 0.000e+00
K = [12.23242  16.69827  93.04242  29.69827]
L = [-1.00000  -4.62535 -12.95097]
Steady-state |e(t)| max (t>150s) = 
Steady-state |e(t)| RMSE(root mean square error. Please note that MSE=mean squre err, RMSE =root MSE ) (t>150s) = NaN






############B4w1_reB3




%% Parameters (given)
M = 1;     % cart mass
F = 1;     % friction coeff
L = 1;     % pendulum length
g = 9.81;

alpha = 1;
omega = 1;  % this is the turnable parameter for B4 --> w=0.1->w =1/10/100

Amp = 0.5;   % d1 square-wave amplitude
T   = 50;    % d1 period [s]

tEnd = 100;   % this is the turnable parameter for B4
opts = odeset('RelTol',1e-9,'AbsTol',1e-11); 
% RelTol : residual tolerance 

FBI residual1 = 4.163e-17, residual2 = 0.000e+00
K (used in mu = Kx + Lw) = [12.23242  16.69827  93.04242  29.69827]
L = [-1.00000  -4.62535 -12.95097]
Nonlinear: steady max |e| (t>150s) = 
Nonlinear: steady RMS |e| (t>150s) = NaN















B4 result

==================== omega = 0.1 ====================
FBI residual1 = 1.880e-19, residual2 = 0.000e+00
Linear:    max|e|=1.351e-08, rms|e|=4.899e-09, rms u=4.383e-01, max|phi|=1.018e-03
Nonlinear: max|e|=1.496e-08, rms|e|=5.326e-09, rms mu=4.390e-01, max|phi|=2.056e-02

==================== omega = 1 ====================
FBI residual1 = 4.163e-17, residual2 = 0.000e+00
Linear:    max|e|=1.554e-15, rms|e|=4.667e-16, rms u=1.120e+00, max|phi|=9.243e-02
Nonlinear: max|e|=1.660e-03, rms|e|=1.011e-03, rms mu=1.098e+00, max|phi|=1.504e-01

==================== omega = 10 ====================
FBI residual1 = 0.000e+00, residual2 = 0.000e+00
Nonlinear stopped early at t=0.695, event=1
Linear:    max|e|=5.440e-15, rms|e|=1.962e-15, rms u=7.060e+01, max|phi|=9.067e-01
Nonlinear: max|e|=6.078e+00, rms|e|=1.795e+00, rms mu=1.576e+02, max|phi|=2.513e+00

==================== B4 DISCUSSION GUIDE ====================
As omega increases, the sinusoidal reference d_2(t)=alpha*v_1(t) becomes faster.
The internal-model/feedforward targets that frequency, but higher omega typically
requires larger control effort (terms scale with omega and omega^2).
The linear model may still show small steady-state error, while the nonlinear
model can deviate due to strong coupling via sin(phi), cos(phi); for large omega,
phi(t) can grow and the simulation may trigger a safety event (large-angle regime).





>> 
