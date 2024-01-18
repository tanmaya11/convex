m = sym('m')
q = sym('q')

psi01 = sym('psi01');
psi02 = sym('psi02');
psi03 = sym('psi03');
psi11 = sym('psi11');
psi12 = sym('psi12');
psi13 = sym('psi13');
psi21 = sym('psi21');
psi22 = sym('psi22');
psi23 = sym('psi23');

% psi01 = 4;
% psi02 = -8;
% psi03 = 14;
% psi11 = 3/2;
% psi12 = -1/2;
% psi13 = -5/2;
% psi21 = 1/4;
% psi22 = -1/8;
% psi23 = -3/8;
% m = 2;q = -4;

t0 = (-psi01-m*psi02)/(2*(psi11+m*psi12));
t1 = 1/(2*(psi11+m*psi12));
t2 = m/(2*(psi11+m*psi12));
gamma10 = t1*(psi23+q*psi22)/(psi11+m*psi12);
gamma01 = t2*(psi23+q*psi22)/(psi11+m*psi12);
gamma00 = (t0*(psi23+q*psi22)-psi13-q*(psi12))/(psi11+m*psi12);
%gamma01 = sym('gamma01');
%gamma10 = sym('gamma10');

zeta11 = -(psi11*gamma10+m*psi12*gamma10)^2/(psi23+q*psi22) + gamma10
zeta12 = -(2*(psi11*gamma01+m*psi12*gamma01)*(psi11*gamma10+m*psi12*gamma10))/(psi23+q*psi22) + gamma01 + m *gamma10
zeta22 = -(psi11*gamma01+m*psi12*gamma01)^2/(psi23+q*psi22) + gamma01 *m

% zeta11 = sym('zeta11');
% zeta12 = sym('zeta12');
% zeta22 = sym('zeta22');
                 
D = zeta12^2 - 4*zeta11*zeta22
simplifyFraction(D)
return
D1 = D * 8 *(psi11+m*psi12)^4
simplify(D)
D0 = simplifyFraction(D)
x1 = sym('x1')
y1 = sym('y1')
f1 = sym('f1')
D0 = D0 * 8 *(psi11+m*psi12)^4
D0 = subs(D0,psi22,-1/(4*m*(q-y1+m*x1)))
D0 = subs(D0,psi23,(y1-m*x1)/(4*m*(q-y1+m*x1)))
simplifyFraction(D0)
return
D1 = D0 * 8 *(psi11+m*psi12)^4
simplifyFraction(D1/m)

x = sym('x')
y = sym('y')
x1 = sym('x1')
y1 = sym('y1')
f1 = sym('f1')


% f = ((psi11 * x + psi12*y + psi13 )^2 /(psi21 * x + psi22*y + psi23 )) + psi01 * x + psi02*y + psi03 
% g(1) = diff(f,x)
% g(2) = diff(f,y)
% H(1,1) = diff(g(1),x)
% H(1,2) = diff(g(1),y)
% H(2,1) = diff(g(2),x)
% H(2,2) = diff(g(2),y)
% simplifyFraction(H(1,2)-H(2,1))
% simplifyFraction(H(1,1)*H(2,2)-H(1,2)^2)
% simplifyFraction(H)
% H

a = sym('a')
b = sym('b')
z = sym('z')

etah = (a+m*b-q)^2/(4*m) - b*q
etaw = f1 -a*x1 - b*y1
s = etah - etaw
s = subs(s,a,z-m*b)
bv = solve(s ,b)
obj = etah + a*x+b*y
obj0 = subs(obj,a,z-m*b)
obj0 = subs(obj0,b,bv)

%a = z-m*b
obj0 = subs(obj0,a,z-m*b)
simplifyFraction(obj0)
z0 = coeffs(obj0,z)
psi2 = -z0(end)
psi2 = simplifyFraction(psi2)
coeffs(psi2,x)
coeffs(psi2,y)
          %obj0 = obj0.subsVarsPartial([a],[av]);
          %obj0 = obj0.subsVarsPartial([b],[bv]);

