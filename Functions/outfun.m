function stop = outfun(x, optimValues, state)
stop = false;
hold on;
plot(-log(abs(x(3))),abs(x(4)/x(1)),'.');
drawnow

