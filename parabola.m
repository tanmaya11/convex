x = sym('x')
y = sym('y')
a = sym('a')
b = sym('b')
c = sym('c')
d = sym('d')
e = sym('e')
f = sym('f')
a1 = sym('a1')
b1 = sym('b1')
c1 = sym('c1')
d1 = sym('d1')
e1 = sym('e1')
f1 = sym('f1')

p1 = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
p2 = a1*x^2 + b1*y^2 + c1*x*y + d1*x + e1*y + f1

% c^2-4*a*b=0
% c = 2*sqrt(a*b)

q1 = subs(p1,c,2*sqrt(a*b))

q2 = subs(p2,c1,2*sqrt(a1*b1))

dp = q1-q2

[cx,tx] = coeffs(dp,[x,y])

simplifyFraction(cx(2)^2 - 4*cx(1)*cx(4))