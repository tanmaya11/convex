syms x y 
%z = (15*x^2-3*x*y+10*y^2+90*x-65*y+75)/(3*x-2*y+25);
%z = 2*y^2/(y+2-x);
%z = 2*y+x-2;
%z = x*y;
%z = (90*x - 65*y - 3*x*y + 15*x^2 + 10*y^2 + 75)/(3*x - 2*y + 25)
%z = (30*x - 30*y - x*y + 5*x^2 + 5*y^2 + 25)/(x - y + 10)
z = x*y
ezsurf(x,y,z,[-5 5 -5 5]);
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