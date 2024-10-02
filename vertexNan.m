% quad - linear case of stp 1 
% rational function with 0/0 at vertex

a= sym('a')
b= sym('b')
m= sym('m')
q= sym('q')
fv = sym('fv')
xv = sym('x_1')
yv = sym('y_1')
fv = xv*yv

dl= sym('dl')
du= sym('du')

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

[cz,terms] = coeffs(obj,z)

cz2 = simplify(cz(1) )

cz0 = simplify(cz(3) )

cz2 = simplify(cz(2) )

psi2 = -simplify(cz(1) )
psi1 = simplify(cz(2) )/2
psi0 = simplify(cz(3) )


% vertex 0/0 for obj1
subs(psi2, [x,y], [xv,yv])
subs(psi1, [x,y], [xv,yv])

obj1 = simplifyFraction(psi1^2/psi2 + psi0)
obj2 = simplifyFraction(-psi2*dl^2+2*dl*psi1+psi0)
obj2 = simplifyFraction(-psi2*du^2+2*du*psi1+psi0)

return

disp('verification for region')


xv2 = sym('x_2')
xv3 = sym('x_3')
yv2 = sym('y_2')
yv3 = sym('y_3')



edge = simplifyFraction((2*m*xv2+q) * psi2 - psi1)
subs(edge,[x,y],[xv1,yv1])
e1=subs(edge,[x,y],[xv2,yv2])
e1=simplifyFraction(subs(e1,[q],[yv2-m*xv2]))

disp('verification for region end')
return

e1=subs(edge,[x,y],[xv3,yv3])
e1=simplifyFraction(subs(e1,[q],[yv3-m*xv3]))



e1 = simplifyFraction(subs(edge,[m],[(yv2-xv2)/(yv1-xv1)]))
%return
e1 = subs(edge,[xv2],[0])
e1 = subs(e1,[xv1,yv1],[2,0])
e1 = subs(e1,[m,q],[1,0])


e1 = subs(edge,[xv2],[1])
e1 = subs(e1,[xv1,yv1],[2,0])
e1 = subs(e1,[m,q],[1,0])
