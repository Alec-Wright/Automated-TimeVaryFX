close all
clear all

chirps = 100;
spc = 20;

fs = 44100;

sps = round(fs*spc/1000);

a = -0.9;

delhf = (1 - a^2)/(1 + 2*a*cos(pi) + a^2);
dellf = (1 - a^2)/(1 + 2*a*cos(0) + a^2);
grpdel = max(delhf,dellf);

out = filter([a, 1], [1, a], [1; zeros(1300,1)]);

in = filter([a, 1], [1, a], flip(out));

M = 64;

x = [0];

for ch = 0:chirps-1
   x(ch*sps + 1) = 1; 
end

x(end+1:end+(round(M*grpdel)) + 50) = 0;

subplot(3,1,1)
plot(x)

for n = 1:M
    x = filter([a, 1], [1, a], x);
end

subplot(3,1,2)
plot(x)

y = flip(x);

for n = 1:M
    y = filter([a, 1], [1, a], y);
end

subplot(3,1,3)
plot(y)