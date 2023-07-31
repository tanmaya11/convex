syms x y 
x = linspace(-6, -4, 100);

y = 4*x/5  + 9;
plot(x,y);
hold on;

y = 3*x  + 20;
plot(x,y);
hold on;


x = linspace(-2, 0, 100);
y = -2*x/3  - 2/3;
plot(x,y);
hold on;

y = 4*x/5  + 4/5;
plot(x,y);
hold on;


x = linspace(0, 2, 100);
y = 3*x  ;
plot(x,y);
hold on;

y = -2*x/3  + 11/3;
plot(x,y);
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