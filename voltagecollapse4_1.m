% 20 October 2016
% Aurik Sarker
clear

% - compute the power consumed by rc(k), i.e. Pc(k)= v(k)^2/rc_(k).
% - compute  alpha(k) = (Pc(k)-Pc(k-1))/(rc(k) - rc(k-1))
%     update rc(k+1) = rc(k+1) + dt*kc*alpha(k)( Pc(k) - Pc_0 )

% Define discrete time period and time step
k = (1:1:1E5);
dt = 1E-4;

% Define starting values
rd0 = 30;
rc0 = 50;
rl = 25;
ep = 10;
pc0 = 0.5;
b = 1;

% Define time arrays for rd, v, p
rd = zeros(1, k(end));
rc = zeros(1, k(end));
v = zeros(1, k(end));
p = zeros(1, k(end));
pc = zeros(1, k(end));
% Different p0
p0P = (.1:.1:1);

% Define time arrays for dpc and drd
dpc = zeros(1, k(end));
drc = zeros(1, k(end));
a = zeros(1, k(end));

% Define cell arrays for multiple starting p0
rdP = cell(1, 10);
rcP = cell(1, 10);
vP = cell(1, 10);
pP = cell(1, 10);
pcP = cell(1, 10);
aP = cell(1, 10);

% Logical array indicating whether collapse has occurred or not
collapse = zeros(1, k(end));
collapseP = cell(1, 10);

% Set zero values for rd and vd
rd(1) = rd0;
rc(1) = rc0;
v(1) = (ep * rd0 * rc0) / (rd0 * (rc0 + rd0) + rl*rc0);
p(1) = v(1)^2 * (1/rc0 + 1/rd0);
pc(1) = v(1)^2 / rc0;

for n = 1:10
    p0 = p0P(n);
    
    dpc(1) = v(1)^2 / rc(1);
    drc(1) = 1;
    a(1) = -dpc(1) / drc(1);
    
    for i = 2:length(k)
        a(i) = -sign(dpc(i-1)) / sign(drc(i-1));
        rc(i) = max(rc(i-1) + dt * rc(i-1) * (pc(i-1) - pc0) * a(i-1), 0);
        rd(i) = max(rd(i-1) + dt*( (v(i-1)^2 / rd(i-1) - p0) * rd(i-1)), 0);
        v(i) = (ep * rd(i-1) * rc(i-1)) / (rd(i-1) * (rc(i-1) + rd(i-1)) + rl*rc(i-1));
        
        p(i) = v(i)^2 * (1/rc(i) + 1/rd(i));
        pc(i) = v(i)^2 / rc(i);
        dpc(i) = pc(i) - pc(i-1);
        drc(i) = rc(i) - rc(i-1);        
    end
    
    collapse = (v < 1E-4) & (rd < rd0);
    
    % Set the corresponding position in the cell array
    rdP{n} = rd;
    rcP{n} = rc;
    vP{n} = v;
    collapseP{n} = collapse; 
    pP{n} = p;
    pcP{n} = pc;
    aP{n} = a;
end

figure(2)
% Loop through each p0
for i = 1:10
    hold on
    
    % Define color vector outside of plot so that color matches 
    c = rand(1,3);
    
%     subplot(411), plot(k, rdP{i}, 'color', c)
%     hold on
%     title('Rd vs k'), xlabel('k'), ylabel('Rd');
%     
%     subplot(412), plot(k, vP{i}, 'color', c)
%     hold on
%     title('v vs k'), xlabel('k'), ylabel('vRd');
%     
%     subplot(413), plot(k, collapseP{i}, 'color', c)
%     hold on
%     title('collapse vs k'), xlabel('k'), ylabel('collapse');
%     ylim([-0.5 1.5]);
%     
%     subplot(414), plot(k, pP{i}, 'color', c)
%     hold on
%     title('p vs k'), xlabel('k'), ylabel('p');

    subplot(321), plot(k, rdP{i}, 'color', c);
    hold on
    title('rd vs k'), xlabel('k'), ylabel('rd');
    
    subplot(322), plot(k, vP{i}, 'color', c);
    hold on
    title('v vs k'), xlabel('k'), ylabel('v');
    
    subplot(323), plot(k, collapseP{i}, 'color', c);
    hold on
    title('collapse vs k'), xlabel('k'), ylabel('collapse');
    ylim([-0.5 1.5]);
    
    subplot(324), plot(k, pP{i}, 'color', c);
    hold on
    title('power vs k'), xlabel('k'), ylabel('power');
    
    subplot(325), plot(k, aP{i}, 'color', c);
    hold on
    title('alpha vs k'), xlabel('k'), ylabel('alpha');
    ylim([-1.5 1.5]);
    
    subplot(326), plot(k, rcP{i}, 'color', c);
    hold on
    title('rc vs k'), xlabel('k'), ylabel('rc');
end