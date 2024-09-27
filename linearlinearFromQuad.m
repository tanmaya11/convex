% quad - linear case of stp 1 
% rational function with 0/0 at vertex
x= sym('x');
y= sym('y');

a= sym('a');
b= sym('b');
m= sym('m');
q= sym('q');


dl= sym('dl');
du= sym('du');


xv1 = sym('x_1');
yv1 = sym('y_1');

fv1 = xv1 * yv1 %sym('fv1')
xv2 = sym('x_2');
yv2 = sym('y_2');

fv2 = xv2 * yv2 %sym('fv2')

xv3 = sym('x_3');
yv3 = sym('y_3');

fv3 = xv3 * yv3 %sym('fv3')


eta(1) = fv1 - a*xv1 - b*yv1;
eta(2) = fv2 - a*xv2 - b*yv2;
eta(3) = fv3 - a*xv3 - b*yv3;

%2-3 is convex edge
% taking obj1 as 1-2 and obj2 as 1-3

n = 0;


obj = eta(1) + (a*x + b*y);
eq = eta(1) - eta(2);

av = solve(eq,a)

obj0 = subs(obj, [a],[av]);
objfacts = coeffs(obj0,b);
if size(objfacts,2) == 1
  eta0 = 0;
  eta1 = objfacts(1);  
else
  eta0 = objfacts(1);
  eta1 = objfacts(2);  
end

    

                   
c = a+m*b - dl
c =  simplifyFraction(subs(c,[a],[av]))
%c = simplifyFraction(subs(c,m,(yv2-yv3)/(xv2-xv3)))
[d,terms] = coeffs(c,b)
bv = simplifyFraction(solve(c,b)         );
%bv - simplifyFraction(subs(bv,m,(yv2-yv3)/(xv2-xv3)))
lb= -inf;
ub = bv;

obj0 = subs(obj,a,av);
obj1 = simplifyFraction(subs(obj0,b,bv))
subs(obj1,[x,y],[xv1,yv1])
subs(obj1,[x,y],[xv2,yv2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('checking')

%%%%%
subs(eta(1),[xv1,yv1],[2,0])

eta(3)
subs(eta(3),[xv3,yv3],[1,1])
% eta(1) = -2*a
% eta(3) = 1 - b - a

%%%%%


obj = eta(1) + (a*x + b*y)
eq = eta(1) - eta(3)

av = solve(eq,a)
%av = subs(av,[xv1,yv1],[2,0])
%subs(av,[xv3,yv3],[1,1])


obj0 = subs(obj, [a],[av]);
objfacts = coeffs(obj0,b);
if size(objfacts,2) == 1
  eta0 = 0;
  eta1 = objfacts(1);  
else
  eta0 = objfacts(1);
  eta1 = objfacts(2);  
end


c =  du - (a+m*b) ;
%c = a+b-2


c =  simplifyFraction(subs(c,[a],[av]));
d = coeffs(c,b);
bv = simplifyFraction(solve(c,b)         );
ub= bv;
lb = -inf;

obj0 = subs(obj,a,av);
obj2 = simplifyFraction(subs(obj0,b,bv))

% obj2 = subs(obj2,[xv1,yv1],[2,0])
% obj2 = subs(obj2,[xv3,yv3],[1,1])
% obj2 = subs(obj2,[m],[1])
% obj2 = subs(obj2,[du],[2])
% obj2 = subs(obj2,[x,y],[0,0])
% return
% 


simplifyFraction(subs(obj1,[x,y],[xv2,yv2]))
simplifyFraction(subs(obj1,[x,y],[xv1,yv1]))

simplifyFraction(subs(obj2,[x,y],[xv1,yv1]))
simplifyFraction(subs(obj2,[x,y],[xv3,yv3]))




obj1


g0 =simplifyFraction(subs(obj1,m,(yv2-yv3)/(xv2-xv3)))
g0 = simplifyFraction(subs(g0,[x,y],[xv3,yv3]))
obj2
g0 = simplifyFraction(subs(obj2,m,(yv2-yv3)/(xv2-xv3)))
g0 = simplifyFraction(subs(g0,[x,y],[xv2,yv2]))





% We have du >= 0 as convex edge f'x >= 0
% and x2-x3 < 0
% So g0 < xy 
% so skip
return
% obj2
% 
% simplifyFraction(subs(obj2,[x,y],[xv2,yv2]))

g = subs(obj2,[xv1,yv1],[2,0])
g = subs(g,[xv2,yv2],[0,0])
g = subs(g,[xv3,yv3],[1,1])
g = subs(g,[m],[1])
g = subs(g,[du],[2])

subs(g,[x,y],[0,0])




