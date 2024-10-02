% quad - quad verification
x = sym('x')
y = sym('y')

mh = sym('mh')
mw = sym('mw')
qh = sym('qh')
qw = sym('qw')
a = sym('a')
b = sym('b')

etah = - (a+mh*b-qh)^2/(4*mh) - b*qh;
etaw = - (a+mw*b-qw)^2/(4*mw) - b*qw;
eq = simplifyFraction(etah-etaw)
[c,t] = coeffs(eq,a)
s = solve(eq,a)
[c,t] = coeffs(s(1),b)

alpha1 = simplifyFraction(c(1))
alpha0 = simplifyFraction(c(2))

obj1 = etah + a*x + b*y
obj1 = simplifyFraction(subs(obj1, a, alpha1*b+alpha0))

obj2 = etaw + a*x + b*y
obj2 = simplifyFraction(subs(obj2, a, alpha1*b+alpha0))

disp("Check 1 : obj1 == obj2")
simplifyFraction(obj1 - obj2)

[c,t] = coeffs(obj1,b)

psi0 = c(3)
psi1 = c(2)/2
psi2 = -c(1)

disp('function')
f0 = simplifyFraction(psi1^2/psi2+psi0)
[c0,t0] = coeffs(f0,[x,y])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('verifying grad = 0')
bv = simplifyFraction(psi1/psi2);
simplifyFraction(-2*psi2*bv+2*psi1)

%%%%%%%

%%%%%%%

disp('verify common point')
eq1 = y - mh*x - qh
eq2 = y - mw*x - qw
cpt = solve([eq1,eq2],[x,y])
f1 = subs(f0,[x,y],[cpt.x,cpt.y])
simplifyFraction(cpt.x*cpt.y - f1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('verifying psi1, psi2, psi0')
disp('uncomment')
 psi2L = (alpha1+mh)^2/(4*mh)
 psi1L = -((alpha1+mh)*(alpha0-qh)/(2*mh)) + y - qh + alpha1*x
 psi0L = -(alpha0-qh)^2/(4*mh)+alpha0*x
% 
  simplifyFraction (psi2-psi2L)
  simplifyFraction (psi1-psi1L)
  simplifyFraction (psi0-psi0L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Example 11')
mhv = 1
mwv = 1/2
qhv = -1 
qwv = 0

% h v1-v3
% w v1-v2
v1 = [2,1]
v2 = [0,0]
v3 = [1,0]

f0 = subs(f0, [mh,mw],[mhv,mwv])
f0 = simplifyFraction(subs(f0, [qh,qw],[qhv,qwv]))
double(subs(f0, [x,y],v1))
double(subs(f0, [x,y],v2))
double(subs(f0, [x,y],v3))
double(subs(f0, [x,y],[0.6,0.3]))


disp('bounds')
fn = x*y
fn = subs(fn,y,mh*x+qh)
dfn = diff(fn,x)
lb = subs(dfn,[x,mh,qh],[v3(1),mhv,qhv])
ub = subs(dfn,[x,mh,qh],[v1(1),mhv,qhv])
dl = a+mh*b-lb
dl = subs(dl,a,alpha1*b+alpha0)
dl = subs(dl,[mh,qh],[mhv,qhv])
dl = subs(dl,[mw,qw],[mwv,qwv])
dl = solve(dl,b)

du = a+mh*b-ub
du = subs(du,a,alpha1*b+alpha0)
du = subs(du,[mh,qh],[mhv,qhv])
du = subs(du,[mw,qw],[mwv,qwv])
du = solve(du,b)

psi1v = subs(psi1,[mh,mw,qh,qw],[mhv,mwv,qhv,qwv])
psi2v = subs(psi2,[mh,mw,qh,qw],[mhv,mwv,qhv,qwv])

ineq1 = simplifyFraction(psi1v/psi2v - dl)

ineq1 = simplifyFraction(psi1v/psi2v - du)

disp('obj')
f1 = subs(obj1,b,dl)
f1 = subs(f1,[mh,qh],[mhv,qhv])
f1 = subs(f1,[mw,qw],[mwv,qwv])

double(subs(f1, [x,y],v1))
double(subs(f1, [x,y],v2))
double(subs(f1, [x,y],v3))
double(subs(f1, [x,y],[0.6,0.3]))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = x*y
fn = subs(fn,y,mw*x+qw)
dfn = diff(fn,x)

lb = subs(dfn,[x,mw,qw],[v2(1),mwv,qwv])
ub = subs(dfn,[x,mw,qw],[v1(1),mwv,qwv])
dl = a+mw*b-lb
dl = subs(dl,a,alpha1*b+alpha0)
dl = subs(dl,[mh,qh],[mhv,qhv])
dl = subs(dl,[mw,qw],[mwv,qwv])
dl = solve(dl,b)

du = a+mw*b-ub
du = subs(du,a,alpha1*b+alpha0)
du = subs(du,[mh,qh],[mhv,qhv])
du = subs(du,[mw,qw],[mwv,qwv])
du = solve(du,b)

psi1v = subs(psi1,[mh,mw,qh,qw],[mhv,mwv,qhv,qwv])
psi2v = subs(psi2,[mh,mw,qh,qw],[mhv,mwv,qhv,qwv])

ineq1 = simplifyFraction(psi1v/psi2v - dl)

ineq1 = simplifyFraction(psi1v/psi2v - du)


return
%%%%%%%%%%%%%%%%%%%%%%%%%%
double(subs(f2,[x,y],[0,0]))