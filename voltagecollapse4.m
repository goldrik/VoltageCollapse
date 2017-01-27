% 14 October 2016
% Jesse Rines
% sweeping P0's
% Reactive collapse prevention
% Plotted to 10E6 (k)


clear, close all
k = (1:1:1E5);

P0 = 0;           % Power Demanded 
rd_0 = 30;        % Rd initial resistance 
dt = 0.0001;       % timestep of calculation
rc_0 = 50;         % Resistance under control (here constant)
rl = 25;          % Load Resistance
e = 10;         % Generator voltage output
P_c = 0.5;        % Power demaded at control resistance
P0_last = 1;      % in discrete P0 sweep, this is highest P0
beta = 1;         % a gain constant 

% defining rd, rc, v, and p and p_c vectors, size same as k, initializing to zero
rd = k;           % resistance connected to P0
rd = rd.*0;       
v = k;            % voltage at node after Rl
v = v.*0;
p = k;            % instantaneous power of circuit
p = p.*0;
p_c = k;          % power at control resistance 
p_c = p_c.*0;
rc = k;           % control resistance
rc = rc.*0;


% defining initial values of rd, v, rc
rd(1) = rd_0;
rc(1) = rc_0;
v(1) = (e * rd(1) * rc_0) / (rd(1) * (rc_0 + rd(1)) + rl*rc_0);


% defining cells with cell index for each rd, v at each P0 iteration
rd_P0 = {.1:.1:P0_last};        % Cell rd; space for each P0
rc_P0 = {.1:.1:P0_last};        % Cell rc; space for each P0
v_P0 = {.1:.1:P0_last};         % Cell v; space for each P0
p_P0 = {.1:.1:P0_last};         % Cell p; space for each P0
p_c_P0 = {.1:.1:P0_last};       % Cell p_c; space for each P0
alpha_P0 = {.1:.1:P0_last};

% initializing boolean vector collapse for all k , P0
collapse = k;                   % vector     
collapse_P0 = {.1:.1:P0_last};  % Cell collapse; space for each P0

% Initialize del_p_c and del_r_d, each functions w/ input k
del_p_c = k;        % change in power at controlled resistor 
del_p_c = 0.*del_p_c;
del_rc = k;         % change in resistance at dependant resistor 
del_rc = 0*del_rc; 
alpha = k;
alpha = alpha.*0;

% loop through P0 and calculate for all k within
for P0 = (1:1:P0_last*10);      % loop through each P0 , 1-->P0_last
    
    del_p_c(1) = (v(1)^2/rc_0);
    del_rc(1) = 1;
    alpha(1) = -1*del_p_c(1)/del_rc(1);
    
    % loop through each k; compute rd(k), v(k)
    for i = 1:(length(k)-1)     % loop variable i 
        alpha(i+1) = -1*sign(del_p_c(i))*sign(del_rc(i));
        rc(i+1) = max(rc(i) + dt*rc(i)*(v(i)^2/rc(i) - P_c)*alpha(i), 0);
        rd(i+1) = max(rd(i) + dt*((v(i)^2/rd(i) - P0/10)*rd(i)), 0);
        v(i+1) = (e * rd(i) * rc(i)) / (rd(i) * (rc(i) + rd(i)) + rl*rc(i));
        del_p_c(i+1) = (v(i+1)^2/rc(i+1)) - (v(i)^2/rc(i));
        del_rc(i+1) = rc(i+1) - rc(i);
    end
    
    % collapse TRUE if voltage < 0.0001 AND rd < r0
    % rd < r0 term to avoid true at very low k
    collapse = lt(v, 0.0001) .* lt(rd, rd_0);
    
    % compute instantaneous power based on v, rc, rd
    p = (v.^2).*(1./rc + 1./rd);
    
    % compute power at rc 
    p_c = (v.^2)./rc;
    
    % assign rd, v, p, collapse to respective 
    % P0 index in appropriate cell
    rd_P0{P0} = rd;
    v_P0{P0} = v;
    collapse_P0{P0} = collapse; 
    p_P0{P0} = p;
    p_c_P0 = p_c;      
    alpha_P0{P0} = alpha; 
    rc_P0{P0} = rc;

end

figure(2)
% loop through each P0 when plotting 
for j = (1:1:P0_last*10)      % j is loop variable
    
    % define color vector c outside of plot; color is same
    % regardless of graph for each P0
    c = (1:1:3);            
    c = rand(1,3);         
    
    % Plot Rd_P0, v_P0, collapse_P0, p_P0:
    % rd_P0 block
    subplot(3,2,1), plot(k, rd_P0{j}, 'color', c);
    hold on;
    title('Rd as fx of k');
    xlabel('k');
    ylabel('Rd (Ohms)');
    
    % v_P0 block
    subplot(3,2,2), plot(k, v_P0{j}, 'color', c);
    hold on;
    title('V as fx of k');
    xlabel('k');
    ylabel('V (Volts)');
    
    % collapse_P0 block
    subplot(3,2,3), plot(k, collapse_P0{j}, 'color', c);
    hold on;
    title('Collapse True/False');
    xlabel('k');
    ylabel('Collapse True/False');
    ylim([-0.5 1.5]);
    
    % p_P0 block 
    subplot(3,2,4), plot(k, p_P0{j}, 'color', c);
    hold on;
    title('Instantaneous Power as fx of k');
    xlabel('k');
    ylabel('Power (W)');
    
    % alpha block
    subplot(3,2,5), plot(k, alpha_P0{j}, 'color', c);
    hold on;
    title('alpha(k)');
    xlabel('k');
    ylabel('alpha(k)');
    
    % rc block
    subplot(3,2,6), plot(k, rc_P0{j}, 'color', c);
    hold on;
    title('rc(k)');
    xlabel('k');
    ylabel('rc(k)');
end


