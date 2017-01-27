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
k = (1:1:5e3);
% k = (1:1:5);
dt = 10E-1;

% Define initial values
rd0 = 40;   % initial rd value
rc0 = 60;   % initial rc value
rl = 10;    % load resistance
ep = 10;    % epsilon (determines maximum power)
g = 150;      % gain constant

samples = 5000;% number of samples used to estimate rc
sigma = 1; % standard deviation of the normal distribution
% Array containing multiple inital p0 values
pmax= ep^2/4/rl;
N=9;
Nini=1;
p0P = (pmax/N/2:pmax/N/2:pmax/2);
p0P = (1.5/N:1.5/N:1.5);
%pd0 = 1.5;
pc0 = pmax*2/3;  % power demanded by rc


% Define cell arrays for multiple initial p0
rdP = cell(1, 10);
rcP = cell(1, 10);
vP = cell(1, 10);
pP = cell(1, 10);
pcP = cell(1, 10);
pdP = cell(1, 10);
aP = cell(1, 10);
bP = cell(1, 10);


% Define time arrays for rd, v, p
rd = zeros(1, k(end));
rc = zeros(1, k(end));
v = zeros(1, k(end));
p = zeros(1, k(end));
pc = zeros(1, k(end));
pd = zeros(1, k(end));

% Define time arrays for dpc and drd
dpc = zeros(1, k(end));
drc = zeros(1, k(end));
a = zeros(1, k(end));
dpd = zeros(1, k(end));
drd = zeros(1, k(end));
b = zeros(1, k(end));
% Logical array indicating whether collapse has occurred or not
collapse = zeros(1, k(end));
collapseP = cell(1, 10);

% Set initial values for r, v, p
rd(1) = rd0*(1+randn(1)*.02);
rc(1) = rc0*(1+randn(1)*.02);
% v(1) = (ep * rd0 * rc0) / (rd0 * (rc0 + rd0) + rl*rc0);
v(1) = ep / (rl/rd(1) + rl/rc(1) + 1);
p(1) = v(1)^2 * (1/rc(1) + 1/rd(1));
pc(1) = v(1)^2 / rc(1);
pd(1) = v(1)^2 / rd(1);
% Define inital delta, alpha values
dpc(1) = v(1)^2 / rc(1);
dpd(1) = v(1)^2 / rd(1);
drc(1) = 1;
ac(1) = -0*sign(dpc(1)) / sign(drc(1));
drd(1) = 1;
ad(1) = -0*sign(dpd(1)) / sign(drd(1));
sigc=[];
sigd=[];

for n = Nini:N
    % Set inital p0 value (from defined p0P array)
    p0 = p0P(n);
    
    % Determine array values over time interval
    for i = 2:length(k)
        
        sigma = max(eps,(abs(rd(i-1))+abs(rc(i-1)))*0.01);
        sigc(i-1) = max(1e-2,1e-2*abs(pc(i-1))); 
        sigd(i-1) = max(1e-2,1e-2*abs(pd(i-1)));
        sigc(i-1) = sigma;
        sigd(i-1) = sigma;
        
        rdR = abs(normrnd(rd(i-1), sigd(i-1), 1, samples));
        rcR = abs(normrnd(rc(i-1), sigc(i-1), 1, samples));
        vdR = ep ./ (rl./rc(i-1) + rl./rdR + 1);
        vcR = ep ./ (rl./rcR + rl./rd(i-1) + 1);
        vR = ep ./ (rl./rcR + rl./rdR + 1);
        pdR = vdR.^2 ./ rdR;
        pcR = vcR.^2 ./ rcR;
        
        %Tracking total power
        pdR = vdR.^2 .*(1./ rdR+1/rc(i-1));
        pcR = vcR.^2 .*(1./ rcR+ 1/rd(i-1));
        pR = vR.^2.*(1./ rcR+ 1./rdR);
        dpR= pR - p(i-1);
        dvR = vR - v(i-1);
        
        
        
        dpdR = pdR - pd(i-1);
        dpcR = pcR - pc(i-1);
        drdR = rdR - rd(i-1);
        drcR = rcR - rc(i-1);
        
        ac(i) = mean(dpcR .* drcR)/(sigc(i-1))^2;
        ad(i) = mean(dpdR .* drdR)/(sigd(i-1))^2;
   
        ac(i) = mean(dpR .* drcR)/(sigma)^2;
        ad(i) = mean(dpR .* drdR)/(sigma)^2;
   
        
        rd(i) = max(rd(i-1) + 0.01*dt*( (v(i-1)^2 / rd(i-1) - p0) * rd(i-1)), 0);
        gd= 1 * (ad(i));
        rd(i) = max(rd(i-1) - max(g)*dt  * (pd(i-1) - p0) * gd, eps);

        % Take rc to be a normal random variable to find alpha
        
        
        % Use this new alpha to calculate rc
        gc = 1*(ac(i));
        rc(i) = max(rc(i-1) - max(g) * dt * (pc(i-1) - pc0) * gc, eps);
        %rc(i) = max(rc(i-1) - 0.01*dt * rc(i-1) * (pc(i-1) - pc0) * sign(a(i)), 0);
        
        % v(i) = (ep * rd(i-1) * rc(i-1)) / (rd(i-1) * (rc(i-1) + rd(i-1)) + rl*rc(i-1));
        v(i) = ep / (rl/rd(i) + rl/rc(i) + 1);
        
        % Set p, pc, delta pc, delta rc, for next loop
        p(i) = v(i)^2 * (1/rc(i) + 1/rd(i));
        pc(i) = v(i)^2 / rc(i);
        dpc(i) = pc(i) - pc(i-1);
        drc(i) = rc(i) - rc(i-1);
        pd(i) = v(i)^2 / rd(i);
        dpd(i) = pd(i) - pd(i-1);
        drd(i) = rd(i) - rd(i-1);
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
    aP{n} = ac;
    bP{n} = ad;
end

figure
set(gcf, 'Position', get(0, 'Screensize'));

% Loop through each p0
for i = Nini:N
    % Define color vector outside of plot so that color matches
    c = rand(1,3);
    
    hold on
    subplot(421), plot(k*dt, rdP{i}, 'color', c, 'LineWidth',2);
    title('rd vs k'), xlabel('k'), ylabel('rd');
    
    hold on
    subplot(422), plot(k*dt, rcP{i}, 'color', c, 'LineWidth',2);
    title('rc vs k'), xlabel('k'), ylabel('rc');
    
    hold on
    subplot(423), plot(k*dt, vP{i}, 'color', c, 'LineWidth',2);
    title('v vs k'), xlabel('k'), ylabel('v');
    
    hold on
    subplot(424), plot(k*dt, pP{i}, 'color', c, 'LineWidth',2);
    title('power vs k'), xlabel('k'), ylabel('power');
    
%     hold on
%     subplot(425), plot(k*dt, collapseP{i}, 'color', c, 'LineWidth',2);
%     title('collapse vs k'), xlabel('k'), ylabel('collapse');
%     ylim([-0.5 1.5]);
%     
    hold on
    subplot(425), plot(k*dt, bP{i}, 'color', c, 'LineWidth',2);
    title('beta(b) vs k'), xlabel('k'), ylabel('beta');
    
    hold on
    subplot(426), plot(k*dt, aP{i}, 'color', c, 'LineWidth',2);
    title('alpha(a) vs k'), xlabel('k'), ylabel('alpha');
   
    hold on
    subplot(427), plot(k*dt, pdP{i}, 'color', c, 'LineWidth',2);
    title('power of rd vs k'), xlabel('k'), ylabel('power');
    
    
    hold on
    subplot(428), plot(k*dt, pcP{i}, 'color', c, 'LineWidth',2);
    title('power of rc vs k'), xlabel('k'), ylabel('power');
    
end

return

x = 0.01:.1:100;
[Rd,Rc]=meshgrid(x);
V = ep./(rl./Rd+rl./Rc + 1);
Pd = V.^2./Rd;
Pc = V.^2./Rc;
P = Pd + Pc;

delta=10;
figure
set(gcf, 'Position', get(0, 'Screensize'));

subplot(241)
surf(Rd(1:delta:end,1:delta:end),Rc(1:delta:end,1:delta:end),Pd(1:delta:end,1:delta:end))
title('Pd')
xlabel('Rd')
ylabel('Rc')
subplot(242)
surf(Rd(1:delta:end,1:delta:end),Rc(1:delta:end,1:delta:end),Pc(1:delta:end,1:delta:end))
xlabel('Rd')
ylabel('Rc')
title('Pc')
subplot(243)
surf(Rd(1:delta:end,1:delta:end),Rc(1:delta:end,1:delta:end),P(1:delta:end,1:delta:end))
title('P')
xlabel('Rd')
ylabel('Rc')
subplot(244)
surf(Rd(1:delta:end,1:delta:end),Rc(1:delta:end,1:delta:end),V(1:delta:end,1:delta:end))
title('V')
xlabel('Rd')
ylabel('Rc')
subplot(245)
contour(Rd,Rc,Pd)
title('Pd')
xlabel('Rd')
ylabel('Rc')
subplot(246)
contour(Rd,Rc,Pc)
xlabel('Rd')
ylabel('Rc')
title('Pc')
subplot(247)
contour(Rd,Rc,P)
title('P')
xlabel('Rd')
ylabel('Rc')
subplot(248)
contour(Rd,Rc,V)
title('V')
xlabel('Rd')
ylabel('Rc')

