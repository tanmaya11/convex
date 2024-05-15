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
%return
s1 = sym('s1')
s2 = sym('s2')

zeta00 = (psi21 + m*psi22)^2  ;
delta1 = -2*(psi13*psi21 - psi11*psi23 + m*psi13*psi22 - m*psi12*psi23 - q*psi11*psi22 + q*psi12*psi21)*(psi01*psi21 + psi11^2 + m*(psi12^2*m + psi01*psi22 + psi02*psi21 + 2*psi11*psi12 + psi02*psi22*m))  ;
delta0 = 2*psi11*psi23*(psi11+m*psi12) -2*psi13*psi21*(psi11+m*psi12)-psi03*psi21*(psi21+m*psi22)+ psi01*psi23*(psi21+m*psi22) -2*m*psi13*psi22*(psi11+m*psi12) + 2*m*psi12*psi23*(psi11+m*psi12) - m * psi03*psi22 * (psi21+m*psi22) + m * psi02*psi23*(psi21+m*psi22 ) + 2*q*psi11*psi22*(psi11+m*psi12)- 2*q*psi12*psi21* (psi11+m*psi12) + q*psi01*psi22*(psi21+m*psi22)-q*psi02*psi21*(psi21+m*psi22);
si1 = 2*(psi21 + m*psi22) * (psi13*psi21 - psi11*psi23 + m*psi13*psi22 - m*psi12*psi23 - q*psi11*psi22 + q*psi12*psi21)*s1 + 2*m*(m*psi22 + psi21) * (psi13*psi21 - psi11*psi23 + m*psi13*psi22 - m* psi12*psi23 - q*psi11*psi22 + q*psi12*psi21)*s2 + delta1;
si1_2 = -(psi21 + m*psi22)*s1 -m *(psi21+m*psi22)*s2+psi01*psi21 +psi11^2 + m*(m*psi12^2 + psi01*psi22+psi02*psi21+2*psi11*psi12+m*psi02*psi22);
si0 = (-psi23*(psi21+m*psi22)-q*psi22*(psi21+m*psi22))*s1 + (q*psi21*(psi21+m*psi22)-m*psi23*(psi21+m*psi22))*s2 + delta0;  

ff = simplifyFraction((si1 / (zeta00 * sqrt(si1_2))) + si0)


a = sym('a')
b = sym('b')
c = sym('c')

fl = a*s1+b*s2+c


% fractional form = linear
temp = simplifyFraction(((si1 / (zeta00 * sqrt(si1_2))))^2 -  (fl - si0)^2)

[cx,tx] = coeffs(temp,[s1,s2]);
tx
coeffs(cx(3),s1)
coeffs(cx(3),s2)
%- 8*(b - psi21^2*q + m^2*psi22*psi23 + m*psi21*psi23 - m*psi21*psi22*q)*(a + psi21*psi23 + m*psi22^2*q + m*psi22*psi23 + psi21*psi22*q) - (2*(a + psi21*psi23 + m*psi22^2*q + m*psi22*psi23 + psi21*psi22*q)^2*(2*psi11^2*psi23^2 + 2*psi13^2*psi21^2 + a*psi03*psi21^5 + c*psi21^4*psi23 + psi03*psi21^6*psi23 - psi01*psi21^5*psi23^2 + 2*m^2*psi12^2*psi23^2 + 2*m^2*psi13^2*psi22^2 - 2*psi11^2*psi21^4*psi23^2 + 2*psi11^2*psi22^2*q^2 + 2*psi12^2*psi21^2*q^2 + a*c*psi21^3 - 4*psi11*psi13*psi21*psi23 - m^5*psi01*psi22^5*psi23^2 - m^6*psi02*psi22^5*psi23^2 - m^5*psi01*psi22^7*q^2 - psi01*psi21^5*psi22^2*q^2 - a*psi01*psi21^4*psi23 + 2*a*psi11*psi13*psi21^4 + a*psi02*psi21^5*q + c*psi21^4*psi22*q + 4*m*psi11*psi12*psi23^2 + 4*m*psi13^2*psi21*psi22 + 2*psi11*psi13*psi21^5*psi23 + 4*psi12*psi13*psi21^2*q + psi02*psi21^6*psi23*q + psi03*psi21^6*psi22*q + 4*psi11^2*psi22*psi23*q - 2*m^2*psi12^2*psi21^4*psi23^2 - 2*m^4*psi11^2*psi22^4*psi23^2 - 2*m^6*psi12^2*psi22^4*psi23^2 - 2*m^4*psi11^2*psi22^6*q^2 - 2*psi11^2*psi21^4*psi22^2*q^2 + a*c*m^3*psi22^3 + a*m^5*psi03*psi22^5 + c*m^4*psi22^4*psi23 + c*m^4*psi22^5*q - 2*a*psi11^2*psi21^3*psi23 - m*psi02*psi21^5*psi23^2 + m^6*psi03*psi22^6*psi23 + m^6*psi03*psi22^7*q + psi02*psi21^6*psi22*q^2 - 10*m^2*psi01*psi21^3*psi22^2*psi23^2 - 10*m^3*psi01*psi21^2*psi22^3*psi23^2 - 10*m^3*psi02*psi21^3*psi22^2*psi23^2 - 10*m^4*psi02*psi21^2*psi22^3*psi23^2 - 8*m^3*psi11^2*psi21*psi22^3*psi23^2 - 8*m^3*psi12^2*psi21^3*psi22*psi23^2 - 8*m^5*psi12^2*psi21*psi22^3*psi23^2 - 10*m^2*psi01*psi21^3*psi22^4*q^2 + 10*m^2*psi02*psi21^4*psi22^3*q^2 - 10*m^3*psi01*psi21^2*psi22^5*q^2 + 10*m^3*psi02*psi21^3*psi22^4*q^2 + 5*m^4*psi02*psi21^2*psi22^5*q^2 - 8*m*psi11^2*psi21^3*psi22^3*q^2 - 8*m^3*psi11^2*psi21*psi22^5*q^2 + 2*m^5*psi12^2*psi21*psi22^5*q^2 + 3*a*c*m^2*psi21*psi22^2 - a*m^4*psi01*psi22^4*psi23 + 5*a*m^4*psi03*psi21*psi22^4 + 2*a*m^4*psi11*psi13*psi22^4 - a*m^5*psi02*psi22^4*psi23 + 2*a*m^5*psi12*psi13*psi22^4 - a*m^4*psi01*psi22^5*q + 2*a*m*psi12^2*psi21^4*q + 4*c*m^3*psi21*psi22^3*psi23 + 4*c*m*psi21^3*psi22^2*q + 4*c*m^3*psi21*psi22^4*q - 2*a*psi11^2*psi21^3*psi22*q - 5*m*psi01*psi21^4*psi22*psi23^2 - 4*m*psi11*psi12*psi21^4*psi23^2 + 6*m^5*psi03*psi21*psi22^5*psi23 + 2*m^5*psi11*psi13*psi22^5*psi23 + 2*m^6*psi12*psi13*psi22^5*psi23 + 6*m*psi03*psi21^5*psi22^2*q - 2*m^5*psi01*psi22^6*psi23*q + 6*m^5*psi03*psi21*psi22^6*q + 2*m^5*psi11*psi13*psi22^6*q - m^6*psi02*psi22^6*psi23*q + 2*m^6*psi12*psi13*psi22^6*q + 2*m*psi12^2*psi21^5*psi23*q + 2*psi11*psi12*psi21^5*psi22*q^2 - 4*psi11^2*psi21^4*psi22*psi23*q - 12*m^2*psi11^2*psi21^2*psi22^2*psi23^2 - 12*m^4*psi12^2*psi21^2*psi22^2*psi23^2 - 12*m^2*psi11^2*psi21^2*psi22^4*q^2 + 8*m^2*psi12^2*psi21^4*psi22^2*q^2 + 12*m^3*psi12^2*psi21^3*psi22^3*q^2 + 8*m^4*psi12^2*psi21^2*psi22^4*q^2 - 4*m*psi11*psi13*psi22*psi23 - 4*m*psi12*psi13*psi21*psi23 - 4*psi11*psi12*psi21*psi23*q - 4*psi11*psi13*psi21*psi22*q + 10*a*m^2*psi03*psi21^3*psi22^2 + 10*a*m^3*psi03*psi21^2*psi22^3 - 2*a*m^2*psi12^2*psi21^3*psi23 - 2*a*m^3*psi11^2*psi22^3*psi23 - 2*a*m^5*psi12^2*psi22^3*psi23 - 2*a*m^3*psi11^2*psi22^4*q + 6*c*m^2*psi21^2*psi22^2*psi23 + 6*c*m^2*psi21^2*psi22^3*q - 5*m^2*psi02*psi21^4*psi22*psi23^2 + 15*m^2*psi03*psi21^4*psi22^2*psi23 - 5*m^4*psi01*psi21*psi22^4*psi23^2 + 20*m^3*psi03*psi21^3*psi22^3*psi23 + 15*m^4*psi03*psi21^2*psi22^4*psi23 - 5*m^5*psi02*psi21*psi22^4*psi23^2 - 4*m^5*psi11*psi12*psi22^4*psi23^2 - 8*m*psi11^2*psi21^3*psi22*psi23^2 - 5*m*psi01*psi21^4*psi22^3*q^2 + 5*m*psi02*psi21^5*psi22^2*q^2 + 15*m^2*psi03*psi21^4*psi22^3*q + 20*m^3*psi03*psi21^3*psi22^4*q - 5*m^4*psi01*psi21*psi22^6*q^2 + 15*m^4*psi03*psi21^2*psi22^5*q + m^5*psi02*psi21*psi22^6*q^2 - 2*m^5*psi11*psi12*psi22^6*q^2 + 2*m*psi12^2*psi21^5*psi22*q^2 - 4*m^4*psi11^2*psi22^5*psi23*q - 2*m^6*psi12^2*psi22^5*psi23*q + 3*a*c*m*psi21^2*psi22 - a*m*psi02*psi21^4*psi23 + 5*a*m*psi03*psi21^4*psi22 + 2*a*m*psi12*psi13*psi21^4 + 4*c*m*psi21^3*psi22*psi23 - a*psi01*psi21^4*psi22*q + 2*a*psi11*psi12*psi21^4*q - 4*m^2*psi12*psi13*psi22*psi23 + 6*m*psi03*psi21^5*psi22*psi23 + 2*m*psi12*psi13*psi21^5*psi23 - 4*m*psi11*psi13*psi22^2*q - 4*m*psi12^2*psi21*psi23*q - 4*psi11*psi12*psi21*psi22*q^2 - 2*psi01*psi21^5*psi22*psi23*q + 2*psi11*psi12*psi21^5*psi23*q + 2*psi11*psi13*psi21^5*psi22*q - 6*a*m^2*psi01*psi21^2*psi22^2*psi23 + 12*a*m^2*psi11*psi13*psi21^2*psi22^2 - 6*a*m^3*psi02*psi21^2*psi22^2*psi23 + 12*a*m^3*psi12*psi13*psi21^2*psi22^2 - 6*a*m^2*psi11^2*psi21*psi22^2*psi23 - 6*a*m^3*psi12^2*psi21^2*psi22*psi23 - 6*a*m^4*psi12^2*psi21*psi22^2*psi23 - 6*a*m^2*psi01*psi21^2*psi22^3*q + 6*a*m^2*psi02*psi21^3*psi22^2*q + 4*a*m^3*psi02*psi21^2*psi22^3*q - 6*a*m*psi11^2*psi21^2*psi22^2*q - 6*a*m^2*psi11^2*psi21*psi22^3*q + 6*a*m^2*psi12^2*psi21^3*psi22*q + 2*a*m^4*psi12^2*psi21*psi22^3*q - 16*m^2*psi11*psi12*psi21^3*psi22*psi23^2 + 20*m^2*psi11*psi13*psi21^3*psi22^2*psi23 + 20*m^3*psi11*psi13*psi21^2*psi22^3*psi23 - 16*m^4*psi11*psi12*psi21*psi22^3*psi23^2 + 20*m^3*psi12*psi13*psi21^3*psi22^2*psi23 + 20*m^4*psi12*psi13*psi21^2*psi22^3*psi23 + 6*m*psi11*psi12*psi21^4*psi22^2*q^2 - 20*m^2*psi01*psi21^3*psi22^3*psi23*q + 20*m^2*psi11*psi13*psi21^3*psi22^3*q + 5*m^2*psi02*psi21^4*psi22^2*psi23*q + 10*m^2*psi12*psi13*psi21^4*psi22^2*q - 20*m^3*psi01*psi21^2*psi22^4*psi23*q + 20*m^3*psi11*psi13*psi21^2*psi22^4*q + 20*m^3*psi12*psi13*psi21^3*psi22^3*q - 6*m^4*psi11*psi12*psi21*psi22^5*q^2 - 5*m^4*psi02*psi21^2*psi22^4*psi23*q + 20*m^4*psi12*psi13*psi21^2*psi22^4*q - 16*m*psi11^2*psi21^3*psi22^2*psi23*q + 6*m^2*psi12^2*psi21^4*psi22*psi23*q - 16*m^3*psi11^2*psi21*psi22^4*psi23*q - 6*m^5*psi12^2*psi21*psi22^4*psi23*q - 4*a*m*psi01*psi21^3*psi22*psi23 - 4*a*m*psi11*psi12*psi21^3*psi23 + 8*a*m*psi11*psi13*psi21^3*psi22 + 4*a*m*psi02*psi21^4*psi22*q + 10*m*psi11*psi13*psi21^4*psi22*psi23 + 4*m*psi02*psi21^5*psi22*psi23*q + 2*m*psi12*psi13*psi21^5*psi22*q + 6*a*m^3*psi12^2*psi21^2*psi22^2*q - 24*m^3*psi11*psi12*psi21^2*psi22^2*psi23^2 + 4*m^2*psi11*psi12*psi21^3*psi22^3*q^2 - 4*m^3*psi11*psi12*psi21^2*psi22^4*q^2 - 24*m^2*psi11^2*psi21^2*psi22^3*psi23*q + 4*m^3*psi12^2*psi21^3*psi22^2*psi23*q - 4*m^4*psi12^2*psi21^2*psi22^3*psi23*q - 4*a*m^2*psi02*psi21^3*psi22*psi23 + 8*a*m^2*psi12*psi13*psi21^3*psi22 - 4*a*m^3*psi01*psi21*psi22^3*psi23 + 8*a*m^3*psi11*psi13*psi21*psi22^3 - 4*a*m^4*psi02*psi21*psi22^3*psi23 - 4*a*m^4*psi11*psi12*psi22^3*psi23 + 8*a*m^4*psi12*psi13*psi21*psi22^3 - 6*a*m*psi11^2*psi21^2*psi22*psi23 - 4*a*m*psi01*psi21^3*psi22^2*q - 4*a*m^3*psi01*psi21*psi22^4*q + a*m^4*psi02*psi21*psi22^4*q - 2*a*m^4*psi11*psi12*psi22^4*q + 10*m^2*psi12*psi13*psi21^4*psi22*psi23 + 10*m^4*psi11*psi13*psi21*psi22^4*psi23 + 10*m^5*psi12*psi13*psi21*psi22^4*psi23 - 10*m*psi01*psi21^4*psi22^2*psi23*q + 10*m*psi11*psi13*psi21^4*psi22^2*q - 10*m^4*psi01*psi21*psi22^5*psi23*q + 10*m^4*psi11*psi13*psi21*psi22^5*q - 4*m^5*psi02*psi21*psi22^5*psi23*q - 6*m^5*psi11*psi12*psi22^5*psi23*q + 10*m^5*psi12*psi13*psi21*psi22^5*q + 4*m*psi11*psi12*psi22*psi23*q + 4*m*psi12*psi13*psi21*psi22*q - 12*a*m^2*psi11*psi12*psi21^2*psi22*psi23 - 12*a*m^3*psi11*psi12*psi21*psi22^2*psi23 - 4*a*m^3*psi11*psi12*psi21*psi22^3*q - 22*m^4*psi11*psi12*psi21*psi22^4*psi23*q - 12*m^2*psi11*psi12*psi21^3*psi22^2*psi23*q - 28*m^3*psi11*psi12*psi21^2*psi22^3*psi23*q + 4*a*m*psi11*psi12*psi21^3*psi22*q + 2*m*psi11*psi12*psi21^4*psi22*psi23*q))/(psi21 + m*psi22)^3

simplifyFraction(cx(2)^2 - 4*cx(1)*cx(4))

%ff-fl




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

