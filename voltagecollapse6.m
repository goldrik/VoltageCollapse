% 3 Nov 2016
% Aurik Sarker & Jesse Rines
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
k = (1:1:1E4);
dt = 1E-3;

% Define initial values
rd0 = rand(1)*100;   % initial rd value, randomized
rc0 = rand(1)*100;    % initial rc value, randomized
rl = 100;    % load resistance; should be in neighborhood of ep
ep = 100;    % epsilon (determines maximum power)
pc0 = .5;   % power demanded by rc
b = 1;      % gain constant

samples = 10000; % number of samples used to estimate rc
sigma = 1; % standard deviation of the normal distribution


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
        vR = ep ./ (rl/rd(i) + rl./rcR + 1);
        pcR = vR.^2 ./ rcR;
        
        dpcR = pcR - pc(i-1);
        drcR = rcR - rc(i-1);
        a(i) = mean(dpcR .* drcR);
        
        % Use this new alpha to calculate rc
        % a(i) = -sign(dpc(i-1)) / sign(drc(i-1));
        rc(i) = max(rc(i-1) + dt * rc(i-1) * (pc(i-1) - pc0) * a(i), 0);
        
        % v(i) = (ep * rd(i-1) * rc(i-1)) / (rd(i-1) * (rc(i-1) + rd(i-1)) + rl*rc(i-1));
        v(i) = ep / (rl/rd(i) + rl/rc(i) + 1);
        
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

% we need to know how much power each resistor is getting (power sharing),
% plot power of rc, power of rd, and how does that defer from the target Pc
% and Pd 

% if no voltage collapse, or if there is enough power that may be
% transferred, the power of each one should be equal to pd0 and pc0, which
% is expected because that's the demanded power

% after you demand more than the network can support, it's important to see
% how much power each is getting. so make plots which include pc and pd as
% well as rc and rd, voltage and power, alpha and collapse logical

% because the rd isnt doing anything 'smart', rd will get ?, rc will get
% remaining available power. understand if power we are getting is correct,
% and why. also alpha. 

% test many initial conditions. 
% ex. fix pc0 and pd0 to values so that there is voltage collapse
% then change the initial conditions of rc and rd (close to 0 up to
% larger than rl)
% function: imshow

% ex. repeat where both demands are controlled. compute random
% perturbations for rd and locate solutions. 

% document, write a report, motivation of this project, results

% change in voltage and power depend on random rd and rc, does this still
% work















