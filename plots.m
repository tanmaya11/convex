figure
%vertices_ineq1 = [-5, -4; 0,-4; 2, 0;2,1;1,3;-5,5];
 vertices_ineq1 = [-4, -3; -3, -3; -1, 1];
% vertices_ineq1 = [0, 0; 2, 0; 2, 1; 1,1;0,0;2,0;1,1];
%vertices_ineq2 = [6, 3; 3, 5; 2, 1];
%return

% Plot the inequalities
%figure;
 xlim([-5, 0]); % Set x-axis range from 0 to 6
 ylim([-2, 2]); % Set y-axis range from 0 to 6
fill(vertices_ineq1(:, 1), vertices_ineq1(:, 2), 'b', 'FaceAlpha', .5);
 return
%hold on;
syms x y 
z = x*y
% %z = (x^2-y^2+x*y)
 fsurf(x,y,z,[0 2 0 1],'r');
 set(h,'edgecolor','none','facecolor',[.1 .9 .1])
 xlim([0, 2]); % Set x-axis range from 0 to 6
 ylim([0, 1]); % Set y-axis range from 0 to 6
 hold on
 return
figure
 syms x y 
 z = -5*y -4*x -20+40
% %z = (x^2-y^2+x*y)
 h=fsurf(x,y,z,[-7 5 -5 5],'r');
 set(h,'edgecolor','none','facecolor',[.1 .9 .1])
 xlim([-5, 0]); % Set x-axis range from 0 to 6
 ylim([-5, 5]); % Set y-axis range from 0 to 6
 hold on
 vertices_ineq1 = [-5, -4, -5; 0, -4, -5; -5, 5, -5];
%vertices_ineq2 = [6, 3; 3, 5; 2, 1];
%return

% Plot the inequalities
%figure;
fill(vertices_ineq1(:, 1), vertices_ineq1(:, 2), 'b', 'FaceAlpha', .5);
hold on;
 return
 
 z = -5*y -4*x -20
  fsurf(-5,0,0,[-7 5 -5 5],'b');
  hold on
  fsurf(0,-4,0,[-7 5 -5 5],'b');
  hold on
  
  %fsurf(x,x-5,0,[-5 5 -5 5],'b');
  hold on

 return
% 

figure
vertices_ineq1 = [2, 1; 4, 0; 6, 3];
vertices_ineq2 = [6, 3; 3, 5; 2, 1];
%return

% Plot the inequalities
%figure;
fill(vertices_ineq1(:, 1), vertices_ineq1(:, 2), 'b', 'FaceAlpha', 0.5);
hold on;
fill(vertices_ineq2(:, 1), vertices_ineq2(:, 2), 'r', 'FaceAlpha', 0.5);
axis equal;
return

syms x y 
%z = (15*x^2-3*x*y+10*y^2+90*x-65*y+75)/(3*x-2*y+25);
%z = 2*y^2/(y+2-x);
%z = 2*y+x-2;
%z = x*y;
%z = (90*x - 65*y - 3*x*y + 15*x^2 + 10*y^2 + 75)/(3*x - 2*y + 25)
%z = (30*x - 30*y - x*y + 5*x^2 + 5*y^2 + 25)/(x - y + 10)
figure
z = x*y
z = (8*x + 6*y - 4*x*y - 2*x^2 + 2*y^2 - 8)/(y - 2*x + 3)
z = - 4*x - 5*y - 20
z = (155*x - 5*y + 4*x*y + 35*x^2 + 5*y^2 - 100)/(7*x - y + 40)
z = (2*x^2-y^2-x*y)
ezsurf(x,y,z,[-5 5 -5 5]);
xlim([-5, 5]); % Set x-axis range from 0 to 6
ylim([-4, 4]); % Set y-axis range from 0 to 6
return


z = -exp(-x^2 ) + exp(x+0.25)
z = (336*x+672*y-12432)/sqrt(3*(-x-2*y+37))
z = (336*x+672*y-12432)/sqrt(3*(-x-2*y+37))
%ezsurf(x,y,z,[0 2 0 2]);
%ezsurf(x,y,z,[-5 1 -4 5]);
ezsurf(x,y,z,[-5 5 -5 5]);
return

vertices_ineq1 = [0, 0; 1, 0; 0, 1];
vertices_ineq2 = [-1, -1; 0, 1; 1, -1];

% Plot the inequalities
figure;
fill(vertices_ineq1(:, 1), vertices_ineq1(:, 2), 'b', 'FaceAlpha', 0.5);
hold on;
fill(vertices_ineq2(:, 1), vertices_ineq2(:, 2), 'r', 'FaceAlpha', 0.5);
axis equal;