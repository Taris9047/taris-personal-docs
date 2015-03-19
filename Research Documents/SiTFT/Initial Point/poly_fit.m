data = load('SiTFTTRANSDK.txt');
xdata = data(:,1);
ydata = data(:,2);

p = polyfit(xdata, ydata, 2);

disp(p);

x = linspace(min(xdata),max(xdata),101);
y = polyval(p,x);

plot(x, y, '-', xdata, ydata, '--')
legend('Fitted polynomial', 'Original Data')
print -deps poly_fit_plot.eps

quit

