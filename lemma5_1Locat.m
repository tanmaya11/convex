a = sym('a')
b = sym('b')
mh = sym('mh')
mw = sym('mw')
qh = sym('qh')
qw = sym('qw')
x = sym('x')
y = sym('y')
alpha1 = sym('alpha1')
alpha0 = sym('alpha0')


eq = (1/mw-1/mh)*a^2 - 2*(qw/mw-qh/mh)*a+(1/mw)*(mw*b+qw)^2-+(1/mh)*(mh*b+qh)^2

av = solve(eq,a)

simplifyFraction(av(1))


disp('objective')

obj = - (a+mh*b-qh)^2/(4*mh) -b*qh + a*x + b * y

obj = subs(obj,a,alpha1*b+alpha0)
obj = expand(obj)

c = coeffs(obj,b)

si1 = -(alpha1+mh)*(alpha0-qh)/(2*mh) + y - qh + alpha1*x

simplify(c(2) - si1)