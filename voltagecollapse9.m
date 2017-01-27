% 30 Nov 2016
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

% Array containing multiple inital values for pd0 (power demanded by rd)
% pd0P = (.1:.1:1.0);
pd0P = (6:1:15);

samples = 1000;        % number of samples used to estimate alpha
sigma = 1;              % standard deviation of the normal distribution


% Define cell arrays
% Each cell corresponds to values from a different initial pd0 value
rdP = cell(1, 10);      % rd values for each pd0 value
rcP = cell(1, 10);      % rc values for each pd0 value
vP = cell(1, 10);       % voltage values for each pd0 value
pP = cell(1, 10);       % power values for each pd0 value
pcP = cell(1, 10);      % power at rc values for each pd0 value
pdP = cell(1, 10);      % power at rd values for each pd0 value
acP = cell(1, 10);      % alpha at rc values for each pd0 value
adP = cell(1, 10);      % alpha at rd values for each pd0 value


% Define time arrays for resistance, voltage, power, alpha
rd = zeros(1, k(end));  % rd values over time
rc = zeros(1, k(end));  % rc values over time
v = zeros(1, k(end));   % voltage values over time
p = zeros(1, k(end));   % power values over time
pc = zeros(1, k(end));  % power at rc values over time
pd = zeros(1, k(end));  % power at rd values over time
ac = zeros(1, k(end));  % alpha at rc values over time
ad = zeros(1, k(end));  % alpha at rd values over time

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
pd(1) = v(1)^2 / rd0;

% Define inital delta, alpha values
% dpc(1) = v(1)^2 / rc(1);
% dpd(1) = v(1)^2 / rd(1);
% drc(1) = 1;
% drd(1) = 1;
% ac(1) = -sign(dpc(1)) / sign(drc(1));
% ad(1) = -sign(dpd(1)) / sign(drd(1));
ac(1) = -sign(rc(1));
ad(1) = -sign(rd(1));


for n = 1:10
    % Set inital pd0 value (from defined pd0P array)
    pd0 = pd0P(n);
    
    % Determine array values over time interval
    for i = 2:length(k)
        % Rd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         rd(i) = max(rd(i-1) + b*dt*( (v(i-1)^2 / rd(i-1) - pd0) ), 0);
        
        % Take rd to be a normal random variable to find alpha
        rdR = normrnd(rd(i-1), sigma, 1, samples);
        vdR = ep ./ (rl/rc(i-1) + rl./rdR + 1);
        pdR = vdR.^2 ./ rdR;
        
        % Use the change in power and resistance at rd to determine ad
%         ad(i) = mean( (pdR - pd(i-1)) .* (rdR - rd(i-1)) );
        
        dpd = pdR - pd(i-1);
        dvd = vdR - v(i-1);
        drd = rdR - rd(i-1);
        ad(i) = mean( dpd./dvd .* dvd./drd + dpd./drd );
        
        % Use this new alpha to calculate rd
        rd(i) = max(rd(i-1) - b*dt * (pd(i-1) - pd0) * ad(i), 0);
        
        
        % Rc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Take rc to be a normal random variable to find alpha
        rcR = normrnd(rc(i-1), sigma, 1, samples);
        vcR = ep ./ (rl./rcR + rl./rd(i) + 1);
        pcR = vcR.^2 ./ rcR;
        
        % Use the change in power and resistance at rc to determine ac
        ac(i) = mean( (pcR - pc(i-1)) .* (rcR - rc(i-1)) );
        
%         dpc = pcR - pc(i-1);
%         dvc = vcR - v(i-1);
%         drc = rcR - rc(i-1);
%         ac(i) = mean( dpc./dvc .* dvc./drc + dpc./drc );
        
        % Use this new alpha to calculate rc
        rc(i) = max(rc(i-1) - b*dt * (pc(i-1) - pc0) * ac(i), 0);
        
        v(i) = ep / (rl/rd(i) + rl/rc(i) + 1);
        
        % Set p, pc, delta pc, delta rc, for next loop
        p(i) = v(i)^2 * (1/rc(i) + 1/rd(i));
        pc(i) = v(i)^2 / rc(i);
        pd(i) = v(i)^2 / rd(i);
        % dpc(i) = pc(i) - pc(i-1);
        % dpd(i) = pd(i) - pd(i-1);
        % drc(i) = rc(i) - rc(i-1);
        % drd(i) = rd(i) - rd(i-1);
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
    pdP{n} = pd;
    acP{n} = ac;
    adP{n} = ad;
end

maxpower = (ep^2 / (4*rl));

figure
set(gcf, 'Position', get(0, 'Screensize'));

hold on
subplot(324), plot(t, pc0 * ones(1, k(end)), '--', 'color', 'k');
hold on
subplot(326), plot(t, maxpower * ones(1, k(end)), '--', 'color', 'k');

% New figures
% Loop through each pd0
for i = 1:10
    % Define color vector outside of plot so that color matches
    c = rand(1,3);
    
    hold on
    subplot(321), plot(t, rdP{i}, 'color', c);
    title('rd vs t'), xlabel('k'), ylabel('rd');
    
    hold on
    subplot(322), plot(t, rcP{i}, 'color', c);
    title('rc vs t'), xlabel('k'), ylabel('rc');
    
    hold on
    subplot(323), plot(t, pd0P(i) * ones(1, k(end)), '--', 'color', c);
    hold on
    subplot(323), plot(t, pdP{i}, 'color', c);
    title('power at rd vs t'), xlabel('t'), ylabel('Prd');
    
    hold on
    subplot(324), plot(t, pcP{i}, 'color', c);
    title('power at rc vs t'), xlabel('t'), ylabel('Prc');
    
    hold on
    subplot(325), plot(t, vP{i}, 'color', c);
    title('v vs t'), xlabel('t'), ylabel('v');
    
    hold on
    subplot(326), plot(t, pP{i}, 'color', c);
    title('power vs t'), xlabel('t'), ylabel('power');
    ylim([0 maxpower + maxpower/10]);
end

% figure
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% % New figures
% % Loop through each pd0
% for i = 1:10
%     % Define color vector outside of plot so that color matches
%     c = rand(1,3);
%     
%     hold on
%     subplot(121), plot(t, acP{i}, 'color', c);
%     title('ac vs t'), xlabel('k'), ylabel('ac');
%     
%     hold on
%     subplot(122), plot(t, adP{i}, 'color', c);
%     title('ad vs t'), xlabel('k'), ylabel('ad');
% end


% test many initial conditions. 
% ex. fix pc0 and pd0 to values so that there is voltage collapse
% then change the initial conditions of rc and rd (close to 0 up to
% larger than rl)
% function: imshow

% ex. repeat where both demands are controlled. compute random
% perturbations for rd and locate solutions. 

% change in voltage and power depend on random rd and rc, does this still
% work


% calculate new alpha which is equal to d/drc (Pc(v(rc,rd),rc)) = dPc/dv * dv/drc + dPc/drc