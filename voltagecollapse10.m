% 8 Dec 2016
% Aurik Sarker & Jesse Rines
clear, close all

% Define discrete time period and time step
k = (1:1:1E4);  % sample scale (more samples = higher runtime)
dt = 1E-0;  % time step (highter time step = longer time scale)
t = dt * k; % actual time scale
b = 10;    % gain constant

% Define initial values
% rd0 = rand(1)*100;      % initial rd value, randomized
% rc0 = rand(1)*100;      % initial rc value, randomized
rd0 = 500;
rc0 = 500;
rl = 100;               % load resistance; should be in neighborhood of ep
ep = 100;               % epsilon (determines maximum power)
pc0 = 13;                % power demanded by rc

delP = 0.0005; % Change in power demanded every time interval
% Array containing multiple inital values for pd0 (power demanded by rd)
pd0P = ( 11 : delP : 11 + delP*(k(end)-1) );

samples = 1000;        % number of samples used to estimate alpha
sigma = 1;              % standard deviation of the normal distribution


% Define time arrays for resistance, voltage, power, alpha
rd = zeros(1, k(end));  % rd values over time
rc = zeros(1, k(end));  % rc values over time
v = zeros(1, k(end));   % voltage values over time
p = zeros(1, k(end));   % power values over time
pc = zeros(1, k(end));  % power at rc values over time
pd = zeros(1, k(end));  % power at rd values over time
ac = zeros(1, k(end));  % alpha at rc values over time
% ad = zeros(1, k(end));  % alpha at rd values over time

% Logical array indicating whether collapse has occurred or not
collapse = zeros(1, k(end));

% Set initial values for r, v, p
rd(1) = rd0;
rc(1) = rc0;
% v(1) = (ep * rd0 * rc0) / (rd0 * (rc0 + rd0) + rl*rc0);
v(1) = ep / (rl/rd0 + rl/rc0 + 1);
p(1) = v(1)^2 * (1/rc0 + 1/rd0);
pc(1) = v(1)^2 / rc0;
pd(1) = v(1)^2 / rd0;

ac(1) = -sign(rc(1));
% ad(1) = -sign(rd(1));


% Determine array values over time interval
for i = 2:length(k)
    pd0 = pd0P(i);
    
    % Rd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rd(i) = max(rd(i-1) + b*dt*( (v(i-1)^2 / rd(i-1) - pd0) ), 0);
    
    % Rc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Take rc to be a normal random variable to find alpha
    rcR = normrnd(rc(i-1), sigma, 1, samples);
    vcR = ep ./ (rl./rcR + rl./rd(i) + 1);
    pcR = vcR.^2 ./ rcR;
    
    % Use the change in power and resistance at rc to determine ac
    ac(i) = mean( (pcR - pc(i-1)) .* (rcR - rc(i-1)) );
    
    % dpc = pcR - pc(i-1);
    % dvc = vcR - v(i-1);
    % drc = rcR - rc(i-1);
    % ac(i) = mean( dpc./dvc .* dvc./drc + dpc./drc );
    
    % Use this new alpha to calculate rc
    rc(i) = max(rc(i-1) - b*dt * (pc(i-1) - pc0) * ac(i), 0);
    
    v(i) = ep / (rl/rd(i) + rl/rc(i) + 1);
    
    % Set p, pc, delta pc, delta rc, for next loop
    p(i) = v(i)^2 * (1/rc(i) + 1/rd(i));
    pc(i) = v(i)^2 / rc(i);
    pd(i) = v(i)^2 / rd(i);
end

% Test collapse logical
collapse = (v < 1E-4) & (rd < rd0);


maxpower = (ep^2 / (4*rl));

figure, plot(t, ad)
title('ad vs t'), xlabel('t'), ylabel('ad');

figure
set(gcf, 'Position', get(0, 'Screensize'));

subplot(321), plot(t, rd);
title('rd vs t'), xlabel('t'), ylabel('rd');

subplot(322), plot(t, rc);
title('rc vs t'), xlabel('t'), ylabel('rc');

subplot(323), plot(t, pd0P, '--');
hold on
subplot(323), plot(t, pd);
title('power at rd vs t'), xlabel('t'), ylabel('Prd');

subplot(324), plot(t, pc);
hold on
subplot(324), plot(t, pc0 .* ones(1, k(end)), '--', 'color', 'k');
title('power at rc vs t'), xlabel('t'), ylabel('Prc');

subplot(325), plot(t, v);
title('v vs t'), xlabel('t'), ylabel('v');

subplot(326), plot(t, p);
hold on
subplot(326), plot(t, maxpower * ones(1, k(end)), '--', 'color', 'k');
title('power vs t'), xlabel('t'), ylabel('power');
ylim([0 maxpower + maxpower/10]);


% start with initial conditions for rc/rd which are close to the fixed points
% (below max power), show that the resistances may converge. This is local
% stability.

% change pd0 small perturbation every time step
% rc is smart, rd is not (use old equation without alpha)