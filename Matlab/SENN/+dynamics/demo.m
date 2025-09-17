x = -2e-2:1e-4:2e-2;

for a = [-10 -5 5 10]
    plot(x, dynamics.vtrap(x, a)); hold on
    plot(x, x ./ (1-exp(a*x)))
end
hold off