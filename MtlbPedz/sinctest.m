rng default

t = 1:10;
x = randn(size(t))';
ts = linspace(0,15,60000);
[Ts,T] = ndgrid(ts,t);
y = sinc(Ts - T)*x;

plot(t,x,'o',ts,y)
xlabel Time, ylabel Signal
legend('Sampled','Interpolated','Location','SouthWest')
legend boxoff