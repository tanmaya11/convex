% quad - linear case of stp 1 
% rational function with 0/0 at vertex

a= sym('a')
b= sym('b')
m= sym('m')
q= sym('q')
fv = sym('fv')
xv = sym('x_1')
yv = sym('y_1')
x = sym('x')
y = sym('y')
fv = xv*yv

dl= sym('dl')
du= sym('du')

[obj1, obj2, obj3] = quadLinear(a, b, m ,q, fv, xv, yv, x, y, dl, du)
