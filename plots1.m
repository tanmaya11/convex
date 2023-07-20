syms x y 
%z = (15*x^2-3*x*y+10*y^2+90*x-65*y+75)/(3*x-2*y+25);
%z = 2*y^2/(y+2-x);

z = 2*y+x-2;
%z = x*y;
ezsurf(x,y,z,[0 2 0 2]);
%ezsurf(x,y,z,[-5 1 0 5]);


%z = x^2+y^2+2*x*y-4*x;
%ezsurf(x,y,z,[0 5 0 5]);
%ezsurf(x,y,z,[-5 -5 0 0]);
hold on

return 
% triangular domain with P(x,y)
P1 = [-5,5];
P2 = [1,3];
P3 = [-1,0];
plot([P1(1);P2(1);P3(1);P1(1)],[P1(2);P2(2);P3(2);P1(2)]);

P1 = [0,0];
P2 = [1,0];
P3 = [0,1];
%plot([P1(1);P2(1);P3(1);P1(1)],[P1(2);P2(2);P3(2);P1(2)]);
