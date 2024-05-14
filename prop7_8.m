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

s1 = sym('s1')
s2 = sym('s2')

% Fractional part 
zeta00 = (psi21 + m*psi22)^2  ;
delta1 = -2*(psi13*psi21 - psi11*psi23 + m*psi13*psi22 - m*psi12*psi23 - q*psi11*psi22 + q*psi12*psi21)*(psi01*psi21 + psi11^2 ...
    + m*(psi12^2*m + psi01*psi22 + psi02*psi21 + 2*psi11*psi12 + psi02*psi22*m))  ;
delta0 = 2*psi11*psi23*(psi11+m*psi12) -2*psi13*psi21*(psi11+m*psi12)-psi03*psi21*(psi21+m*psi22)+ psi01*psi23*(psi21+m*psi22) ...
    -2*m*psi13*psi22*(psi11+m*psi12) + 2*m*psi12*psi23*(psi11+m*psi12) - m * psi03*psi22 * (psi21+m*psi22) + m * psi02*psi23*(psi21+m*psi22 ) ...
    + 2*q*psi11*psi22*(psi11+m*psi12)- 2*q*psi12*psi21* (psi11+m*psi12) + q*psi01*psi22*(psi21+m*psi22)-q*psi02*psi21*(psi21+m*psi22);
si1 = 2*(psi21 + m*psi22) * (psi13*psi21 - psi11*psi23 + m*psi13*psi22 - m*psi12*psi23 - q*psi11*psi22 + q*psi12*psi21)*s1 + ...
    2*m*(m*psi22 + psi21) * (psi13*psi21 - psi11*psi23 + m*psi13*psi22 - m* psi12*psi23 - q*psi11*psi22 + q*psi12*psi21)*s2 + delta1;
si1_2 = -(psi21 + m*psi22)*s1 -m *(psi21+m*psi22)*s2+psi01*psi21 +psi11^2 + m*(m*psi12^2 + psi01*psi22+psi02*psi21+2*psi11*psi12+m*psi02*psi22);
si0 = (-psi23*(psi21+m*psi22)-q*psi22*(psi21+m*psi22))*s1 + (q*psi21*(psi21+m*psi22)-m*psi23*(psi21+m*psi22))*s2 + delta0;  

a = sym('a')
b = sym('b')
c = sym('c')

% Linear
fl = a*s1+b*s2+c

% fractional form = linear
disp('frac-lin')

temp = simplifyFraction(((si1^2 / (zeta00^2 * si1_2))) -  (fl - si0)^2)
[cx,tx] = coeffs(temp,[s1,s2])
simplifyFraction(cx(2)^2 - 4*cx(1)*cx(4))


