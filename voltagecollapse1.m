% 03 October 2016
clear, close all
k = (1:1:100000);

% P0 = 32400;
P0 = 10;
r0 = 300;
dt = 0.00001;
rc = 300;
rl = 30;
e = 9.75;

rd = (1:1:100000);
rd = rd.*0;


v = (1:1:100000);
v = v.*0;
v(1) = (e * rd(1) * rc) / (rd(1) * (rc + rd(1)) + rl*rc);

rd(1) = r0;
for i = 1:99999
    rd(i+1) = max( rd(i) + dt*( ( v(i)^2 / rd(i) - P0) * rd(i)) , 0);
    v(i+1) = (e * rd(i) * rc) / (rd(i) * (rc + rd(i)) + rl*rc);
end

figure
plot(k, rd), title('');
figure 
plot(k, v), title('');

