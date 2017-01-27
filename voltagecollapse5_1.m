% 27 October 2016
% Aurik Sarker
clear, close all

% Compute r_d(k+1)
% But before computing r_c(k+1)
% Pick random values for Delta r_c with nomal distribution N(0,sigma) (say 10)
% For each one of those values, compute 
% Delta P_c (use r_d(k+1) in the computation of the new voltage).
% 
% Now compute
% the expected value:  E[(Delta P_c) (Delta r_c)]
% and define \alpha = E[(Delta P_c) (Delta r_c)]/sigma


% Define discrete time period and time step
k = (1:1:1E5);
dt = 1E-2;

% Define initial values
rd0 = 30;   % initial rd value
rc0 = 50;   % initial rc value
rl = 10;    % load resistance
ep = 10;    % epsilon (determines maximum power)
pc0 = 0.5;  % power demanded by rc
b = 1;      % gain constant

samples = 10;% number of samples used to estimate rc
sigma = 5; % standard deviation of the normal distribution


% Define cell arrays for multiple initial p0
rdP = cell(1, 10);
rcP = cell(1, 10);
vP = cell(1, 10);
pP = cell(1, 10);
pcP = cell(1, 10);
aP = cell(1, 10);


% Define time arrays for rd, v, p
rd = zeros(1, k(end));
rc = zeros(1, k(end));
v = zeros(1, k(end));
p = zeros(1, k(end));
pc = zeros(1, k(end));
% Array containing multiple inital p0 values
pd0P = (0.5:.1:1.5);

% Define time arrays for dpc and drd
dpc = zeros(1, k(end));
drc = zeros(1, k(end));
a = zeros(1, k(end));

% Logical array indicating whether collapse has occurred or not
collapse = zeros(1, k(end));
collapseP = cell(1, 10);

% Set initial values for r, v, p
rd(1) = rd0;
rc(1) = rc0;
% v(1) = (ep * rd0 * rc0) / (rd0 * (rc0 + rd0) + rl*rc0);
v(1) = ep / (rl/rd0 + rl/rc0 + 1);
p(1) = v(1)^2 * (1/rc0 + 1/rd0);
pc(1) = v(1)^2 / rc0;

% Define inital delta, alpha values
dpc(1) = v(1)^2 / rc(1);
drc(1) = 1;
a(1) = -sign(dpc(1)) / sign(drc(1));

for n = 1:10
    % Set inital p0 value (from defined p0P array)
    pd0 = pd0P(n);
    
    % Determine array values over time interval
    for i = 2:length(k)
        rd(i) = max(rd(i-1) + dt*( (v(i-1)^2 / rd(i-1) - pd0) * rd(i-1)), 0);
        
        % Take rc to be a normal random variable to find alpha
        rcR = normrnd(rc(i-1), sigma, 1, samples);
        vR = ep ./ (rl/rd(i-1) + rl./rcR + 1);
        pR = vR.^2 ./ rcR;
        
        dpcR = pR - dpc(i-1);
        drcR = rcR - rc(i-1);
        a(i) = mean(dpcR .* drcR);
        
        % Use this new alpha to calculate rc
        % a(i) = -sign(dpc(i-1)) / sign(drc(i-1));
        rc(i) = max(rc(i-1) + dt * rc(i-1) * (pc(i-1) - pc0) * a(i-1), 0);
        
        % v(i) = (ep * rd(i-1) * rc(i-1)) / (rd(i-1) * (rc(i-1) + rd(i-1)) + rl*rc(i-1));
        v(i) = ep / (rl/rd(i-1) + rl/rc(i-1) + 1);
        
        % Set p, pc, delta pc, delta rc, for next loop
        p(i) = v(i)^2 * (1/rc(i) + 1/rd(i));
        pc(i) = v(i)^2 / rc(i);
        dpc(i) = pc(i) - pc(i-1);
        drc(i) = rc(i) - rc(i-1);
    end
    
    % Test collapse logical
    collapse = (v < 1E-4) & (rd < rd0);
    
    % Set the corresponding position in the cell arrays
    rdP{n} = rd;
    rcP{n} = rc;
    vP{n} = v;
    collapseP{n} = collapse; 
    pP{n} = p;
    pcP{n} = pc;
    aP{n} = a;
end

figure
set(gcf, 'Position', get(0, 'Screensize'));

% Loop through each p0
for i = 1:10
    % Define color vector outside of plot so that color matches
    c = rand(1,3);
    
    hold on
    subplot(321), plot(k, rdP{i}, 'color', c);
    title('rd vs k'), xlabel('k'), ylabel('rd');
    
    hold on
    subplot(322), plot(k, rcP{i}, 'color', c);
    title('rc vs k'), xlabel('k'), ylabel('rc');
    
    hold on
    subplot(323), plot(k, vP{i}, 'color', c);
    title('v vs k'), xlabel('k'), ylabel('v');
    
    hold on
    subplot(324), plot(k, pP{i}, 'color', c);
    title('power vs k'), xlabel('k'), ylabel('power');
    
    hold on
    subplot(325), plot(k, collapseP{i}, 'color', c);
    title('collapse vs k'), xlabel('k'), ylabel('collapse');
    ylim([-0.5 1.5]);
    
    hold on
    subplot(326), plot(k, aP{i}, 'color', c);
    title('alpha(a) vs k'), xlabel('k'), ylabel('alpha');
end