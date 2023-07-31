

syms x y 
y = linspace(-4, 5, 100);
x = zeros(size(y));
x = 0*y-5;
plot(x,y);
line_label = 'x <=  0';
%text(-1, 0, line_label, 'HorizontalAlignment', 'left');
hold on;




x = linspace(-5, 1, 100);
y = (10-x)/3;
plot(x,y);
%line_label = 'x <=  0';
%text(-1, 0, line_label, 'HorizontalAlignment', 'left');
hold on;

%return


x = linspace(-1, 1, 100);
y = (3*x+3)/2;
plot(x,y);
%line_label = 'x+y <= 0';
%text(-3, 0, line_label, 'HorizontalAlignment', 'left');
hold on;


x = linspace(-5, -1, 100);
y = (x+1);
plot(x,y);
%line_label = 'x+y <= 0';
%text(-3, 0, line_label, 'HorizontalAlignment', 'left');
hold on;

x = linspace(-5, 1, 100);
y = (7*x+11)/6;
%plot(x,y);
%line_label = 'x+y <= 0';
%text(-3, 0, line_label, 'HorizontalAlignment', 'left');
hold on;

x = linspace(-5, -1, 100);
y = (-5*x-5)/4;
plot(x,y);
%line_label = 'x+y <= 0';
%text(-3, 0, line_label, 'HorizontalAlignment', 'left');
hold on;
