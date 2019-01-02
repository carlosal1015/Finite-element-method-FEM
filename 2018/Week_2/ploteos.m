fplot(@(x) exp(x)*sin(x), [-10 10], 'Linewidth',2);
hold on
fplot(@(x) 0.5*x ,  [-10 10] ,'--or');
hold off