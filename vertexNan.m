% quad - linear case of stp 1 
% rational function with 0/0 at vertex

a= sym('a')
b= sym('b')
m= sym('m')
q= sym('q')
fv = sym('fv')
xv = sym('xv')
yv = sym('yv')

etah = -(a+m*b-q)^2/(4*m)-b*q
etaw = fv - a*xv - b*yv

eq = etah - etaw

z = sym('z')

% z = a+mb

eq2 = subs(eq, a, z-m*b)

bv = solve(eq2,b)

obj = etaw + a*x+ b*y

obj = subs(obj, a, z-m*b)

obj = subs(obj, b, bv)

cz = coeffs(obj,z)

cz2 = simplify(cz(2) )

subs(cz2, [x,y], [xv,yv])