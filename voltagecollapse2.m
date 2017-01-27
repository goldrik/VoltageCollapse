% 13 October 2016
% Jesse Rines
% sweeping P0's
% Plotted to 10E6 (k)


clear, close all
k = (1:1:1E6);

P0 = 0;           % Power Demanded 
r0 = 30;          % Rd initial resistance 
dt = 0.0001;      % timestep of calculation
rc = 100;         % Resistance under control (here constant)
rl = 30;          % Load Resistance
e = 9.75;         % Generator voltage output
P0_last = 1;      % in discrete P0 sweep, this is highest P0

% defining rd, v, and p vectors, size same as k, initializing to zero
rd = k;           % resistance connected to P0
rd = rd.*0;       
v = k;            % voltage at node after Rl
v = v.*0;
p = k;            % instantaneous power of circuit
p = p.*0;

% defining initial values of rd, v
rd(1) = r0;
v(1) = (e * rd(1) * rc) / (rd(1) * (rc + rd(1)) + rl*rc);


% defining cells with cell index for each rd, v at each P0 iteration
rd_P0 = {.1:.1:P0_last};        % Cell rd; space for each P0
v_P0 = {.1:.1:P0_last};         % Cell v; space for each P0
p_P0 = {.1:.1:P0_last};         % Cell p; space for each P0

% initializing boolean vector collapse for all k , P0
collapse = k;                   % vector     
collapse_P0 = {.1:.1:P0_last};  % Cell collapse; space for each P0

% loop through P0 and calculate for all k within
for P0 = (1:1:P0_last*10);      % loop through each P0 , 1-->P0_last
    
    % loop through each k; compute rd(k), v(k)
    for i = 1:(length(k)-1)     % loop variable i 
        rd(i+1) = max(rd(i) + dt*(( v(i)^2 / rd(i) - P0/10) * rd(i)), 0);
        v(i+1) = (e * rd(i) * rc) / (rd(i) * (rc + rd(i)) + rl*rc);
    end
    
    % collapse TRUE if voltage < 0.0001 AND rd < r0
    % rd < r0 term to avoid true at very low k
    collapse = lt(v, 0.0001) .* lt(rd, r0);
    
    % compute instantaneous power based on v, rc, rd
    p = (v.^2).*(1./rc + 1./rd);
    
    % assign rd, v, p, collapse to respective 
    % P0 index in appropriate cell
    rd_P0{P0} = rd;
    v_P0{P0} = v;
    collapse_P0{P0} = collapse; 
    p_P0{P0} = p;
end

% create figure
figure

% loop through each P0 when plotting 
for j = (1:1:P0_last*10)      % j is loop variable
    
    % define color vector c outside of plot; color is same
    % regardless of graph for each P0
    c = (1:1:3);            
    c = rand(1,3);         
    
    % Plot Rd_P0, v_P0, collapse_P0, p_P0:
    
    % rd_P0 block
    subplot(4,1,1), plot(k, rd_P0{j}, 'color', c);
    hold on;
    title('Rd as fx of k');
    xlabel('k');
    ylabel('Rd (Ohms)');
    
    % v_P0 block
    subplot(4,1,2), plot(k, v_P0{j}, 'color', c);
    hold on;
    title('V as fx of k');
    xlabel('k');
    ylabel('V (Volts)');
    
    % collapse_P0 block
    subplot(4,1,3), plot(k, collapse_P0{j}, 'color', c);
    hold on;
    title('Collapse True/False');
    xlabel('k');
    ylabel('Collapse True/False');
    ylim([-0.5 1.5]);
    
    % p_P0 block 
    subplot(4,1,4), plot(k, p_P0{j}, 'color', c);
    hold on;
    title('Instantaneous Power as fx of k');
    xlabel('k');
    ylabel('Power (W)');
end



