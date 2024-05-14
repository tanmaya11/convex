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

t0 = (-psi01-m*psi02)/(2*(psi11+m*psi12));
t1 = 1/(2*(psi11+m*psi12));
t2 = m/(2*(psi11+m*psi12));
gamma10 = t1*(psi23+q*psi22)/(psi11+m*psi12);
gamma01 = t2*(psi23+q*psi22)/(psi11+m*psi12);
gamma00 = (t0*(psi23+q*psi22)-psi13-q*(psi12))/(psi11+m*psi12);

zeta11 = -(psi11*gamma10+m*psi12*gamma10)^2/(psi23+q*psi22) + gamma10
zeta12 = -(2*(psi11*gamma01+m*psi12*gamma01)*(psi11*gamma10+m*psi12*gamma10))/(psi23+q*psi22) + gamma01 + m *gamma10
zeta22 = -(psi11*gamma01+m*psi12*gamma01)^2/(psi23+q*psi22) + gamma01 *m

zeta10 = -2*(psi11*gamma01+m*psi12*gamma10)*(psi13+psi11*gamma00+psi12*(q+m*gamma00))/(psi23+q*psi22) - m*psi02*gamma10 + gamma00 - psi01*gamma10;
zeta01 = -(2*(psi11*gamma01+m*psi12*gamma01)*(psi13+psi11*gamma00+psi12*(q+m*gamma00)))/((psi23+q*psi22)) - m*psi02*gamma01 - psi01*gamma01 + m*gamma00+q;
zeta00 = -(psi13+psi11*gamma00 +psi12*(q+m*gamma00))^2/(psi23+q*psi22) -psi03 - psi01*gamma00 - psi02*(q+m*gamma00);

D = zeta12^2 - 4*zeta11*zeta22
simplifyFraction(D)

s1 = sym('s1')
s2 = sym('s2')

disp('quad')
fq = simplifyFraction(zeta11*s1^2 + zeta12*s1*s2 + zeta22*s2^2   + zeta10*s1 + zeta01*s2 + zeta00)

[cx,tx] = coeffs(fq,[s1,s2])

delta = simplifyFraction(cx(2)^2 - 4*cx(1)*cx(4))
return

mb = sym('mb')
qb = sym('qb')

psi01b = sym('psi01b');
psi02b = sym('psi02b');
psi03b = sym('psi03b');
psi11b = sym('psi11b');
psi12b = sym('psi12b');
psi13b = sym('psi13b');
psi21b = sym('psi21b');
psi22b = sym('psi22b');
psi23b = sym('psi23b');

t0 = (-psi01b-mb*psi02b)/(2*(psi11b+mb*psi12b));
t1 = 1/(2*(psi11b+mb*psi12b));
t2 = mb/(2*(psi11b+mb*psi12b));
gamma10b = t1*(psi23b+qb*psi22b)/(psi11b+mb*psi12b);
gamma01b = t2*(psi23b+qb*psi22b)/(psi11b+mb*psi12b);
gamma00b = (t0*(psi23b+qb*psi22b)-psi13b-qb*(psi12b))/(psi11b+mb*psi12b);

zeta11 = -(psi11b*gamma10b+mb*psi12b*gamma10b)^2/(psi23b+qb*psi22b) + gamma10b
zeta12 = -(2*(psi11b*gamma01b+mb*psi12b*gamma01b)*(psi11b*gamma10b+mb*psi12b*gamma10b))/(psi23b+qb*psi22b) + gamma01b + mb *gamma10b
zeta22 = -(psi11b*gamma01b+mb*psi12b*gamma01b)^2/(psi23b+qb*psi22b) + gamma01b *mb

zeta10 = -2*(psi11b*gamma01b+mb*psi12b*gamma10b)*(psi13b+psi11b*gamma00b+psi12b*(qb+mb*gamma00b))/(psi23b+qb*psi22b) - mb*psi02b*gamma10b + gamma00b - psi01b*gamma10b;
zeta01 = -(2*(psi11b*gamma01b+mb*psi12b*gamma01b)*(psi13b+psi11b*gamma00b+psi12b*(qb+mb*gamma00b)))/((psi23b+qb*psi22b)) - mb*psi02b*gamma01b - psi01b*gamma01b + mb*gamma00b+qb;
zeta00 = -(psi13b+psi11b*gamma00b +psi12b*(qb+mb*gamma00b))^2/(psi23b+qb*psi22b) -psi03b - psi01b*gamma00b - psi02b*(qb+mb*gamma00b);

D = zeta12^2 - 4*zeta11*zeta22
simplifyFraction(D)

s1 = sym('s1')
s2 = sym('s2')

disp('quad')
fqb = simplifyFraction(zeta11*s1^2 + zeta12*s1*s2 + zeta22*s2^2   + zeta10*s1 + zeta01*s2 + zeta00)

%fq - fqb

%[cx,tx] = coeffs(fq-fqb,[s1,s2])
%delta = simplifyFraction(cx(2)^2 - 4*cx(1)*cx(4))
