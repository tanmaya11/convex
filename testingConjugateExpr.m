x = sym('x');
y = sym('y');

s1 = sym('s1');
s2 = sym('s2');


a = sym('a');
b = sym('b');
c = sym('c');
d = sym('d');
e = sym('e');
f0 = sym('f0');

af = sym('af');
bf = sym('bf');
cf = sym('cf');
df0 = sym('df0');
ef = sym('ef');
ff = sym('ff');
gf = sym('gf');
hf = sym('hf');
jf = sym('jf');


% % For step 4 Fract
% edge = -(a*x+b*y)^2 + d*x + e*y +f0;
% f = sqrt(af*x+bf*y+cf)+(df0*x+ef*y+gf)
% %f = subs(f,bf,sqrt(4*af*cf));
% conj = conjugateExpr(edge,f,x,y)
% return
% %%%%%%%%%%

% For step 4
%edge = x^2 + b*x*y + (b^2/4)*y^2 + d*x + e*y +f0;
%edge = subs(edge,b,sqrt(4*a*c));
edge = (af*x+bf*y)^2 + d*x + e*y +f0;
f = jf*(af*x+bf*y)^2+(af*x-bf*y) +cf
%f = subs(f,bf,sqrt(4*af*cf));
conj = conjugateExpr(edge,f,x,y)
[num,den] = numden(conj)
 [coef,terms] =coeffs(num,[x,y])
 [coef,terms] =coeffs(den,[x,y])
return
%%%%%%%%%%
% Checking example3
f = (-2*x^2+2*y^2-4*x*y+8*x+6*y-8)/(y-2*x+3)
edge = x-y/2-2
conj = conjugateExpr(edge,f,x,y)
return
%%%%%%%%%%

%%%%%%%%%%
% % Checking example2
% f = (35*x^2+5*y^2+4*x*y+155*x-5*y-100)/(7*x-y+40)
% edge = x-y/7-4/7
% conj = conjugateExpr(edge,f,x,y)
% return
%%%%%%%%%%
%%%%%%%%%%
% Checking example
f = (2*y^2)/(y-x+2)
edge = y-x
conj = conjugateExpr(edge,f,x,y)
return
%%%%%%%%%%







% For step 2
 disp('step2')
 edge = d*x + e*y +f0;

 %f = (af*x+bf*y)^2 /(df0*x+ef*y+ff)  + (gf*x+hf*y+jf)
 % imposing psi21 = -m psi22
% df0 = -mef
% df0 = -d/e * ef

 f = (af*x+bf*y+cf)^2 /((-d/e)*ef*x+ef*y+ff)  + (gf*x+hf*y+jf)
 f = simplifyFraction(f);
 conj = conjugateExpr(edge,f,x,y)

