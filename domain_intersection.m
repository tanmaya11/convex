x = linspace(-6, 6, 100);
y = -x;
plot(x,y);
%line_label = 'x+y <= 0';
%text(-3, 0, line_label, 'HorizontalAlignment', 'left');
hold on;

x = linspace(-6, 6, 100);
%y = zeros(size(x));
y = x+1;
plot(x,y);
%line_label = 'x <=  0';
%text(-1, 0, line_label, 'HorizontalAlignment', 'left');
hold on;

syms x y 
y = linspace(-6, 6, 100);
x = zeros(size(y));
%x = y+1
plot(x,y);
line_label = 'x <=  0';
%text(-1, 0, line_label, 'HorizontalAlignment', 'left');
hold on;

syms x y 
y = linspace(-6, 6, 100);
%x = zeros(size(y));
x = 0*y+1
plot(x,y);
line_label = 'x <=  0';
%text(-1, 0, line_label, 'HorizontalAlignment', 'left');
hold on;


x = linspace(-5, 5, 100);
%y = zeros(size(x));
y = 2-x;
plot(x,y);
hold on;

ezplot('8*x-(x+y)^2')
hold on

syms x y 
x = linspace(-6, 6, 100);
%x = zeros(size(y));
y = 0*x+2
plot(x,y);
line_label = 'x <=  0';
%text(-1, 0, line_label, 'HorizontalAlignment', 'left');
hold on;