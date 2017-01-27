% 13 October 2016
% Aurik Sarker
clear, close all

% Define discrete time period and time step
k = (1:1:1E5);
dt = 1E-4;

% Define starting values
rd0 = 30;
rc = 100;
rl = 30;
ep = 9.75;

% Define time arrays for rd, v, p
rd = zeros(1, k(end));
v = zeros(1, k(end));
p = zeros(1, k(end));
% Different p0
p0P = (.1:.1:1);

% Define cell arrays for multiple starting p0
rdP = cell(1, 10);
vP = cell(1, 10);
pP = cell(1, 10);

% Logical array indicating whether collapse has occurred or not
collapse = zeros(1, k(end));
collapseP = cell(1, 10);

% Set zero values for rd and vd
rd(1) = rd0;
v(1) = (ep * rd(1) * rc) / (rd(1) * (rc + rd(1)) + rl*rc);

for n = 1:10
    p0 = p0P(n);
    
    for i = 2:length(k)
        rd(i) = max(rd(i-1) + dt*( (v(i-1)^2 / rd(i-1) - p0) * rd(i-1)), 0);
        v(i) = (ep * rd(i-1) * rc) / (rd(i-1) * (rc + rd(i-1)) + rl*rc);
    end
    
    collapse = (v < 1E-4) & (rd < rd0);
    p = v.^2 .* (1./rc + 1./rd);
    
    % Set the corresponding position in the cell array
    rdP{n} = rd;
    vP{n} = v;
    collapseP{n} = collapse; 
    pP{n} = p;
end

figure
% Loop through each p0
for i = 1:10
    hold on
    
    % Define color vector outside of plot so that color matches 
    c = (1:1:3);            
    c = rand(1,3);         
    
    subplot(411), plot(k, rdP{i}, 'color', c)
    hold on
    title('Rd vs k'), xlabel('k'), ylabel('Rd');
    
    subplot(412), plot(k, vP{i}, 'color', c)
    hold on
    title('v vs k'), xlabel('k'), ylabel('vRd');
    
    subplot(413), plot(k, collapseP{i}, 'color', c)
    hold on
    title('collapse vs k'), xlabel('k'), ylabel('collapse');
    ylim([-0.5 1.5]);
    
    subplot(414), plot(k, pP{i}, 'color', c)
    hold on
    title('p vs k'), xlabel('k'), ylabel('p');
end