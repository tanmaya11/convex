x = sym('x');
y = sym('y');
x_1 = sym('x_1');
y_1 = sym('y_1');
m = sym('m');
q = sym('q');
t = sym('t');
s1 = sym('s_1')
s2 = sym('s_2')

psi2 = (y - y_1 - m*x + m*x_1)/(4*m*(q - y_1 + m*x_1))
 
psi1 = (q*y - q*y_1 + m*q*x - m*q*x_1 - 2*m*x*y_1 + 2*m*x_1*y)/(4*m*(q - y_1 + m*x_1))
 
psi0 = (q^2*y_1 - q^2*y + m*q^2*x - m*q^2*x_1 + 4*m*q*x_1*y_1 - 4*m*x_1*y*y_1 + 4*m^2*x*x_1*y_1)/(4*m*(q - y_1 + m*x_1))

vars = [x,y]
cpsi2 = coeffs (psi2, vars)
cpsi0 = coeffs (psi0, vars)
cpsi1 = coeffs (psi1, vars)
              
vs1 = s1 - (2*cpsi1(1)*t - cpsi2(1)*t^2 + cpsi0(1));
vs2 = s2 - (2*cpsi1(2)*t - cpsi2(2)*t^2 + cpsi0(2));

vt = solve (cpsi2(2)*vs1 - cpsi2(1)*vs2, t );    %  cpsi2(2)*vs1 - cpsi2(1)*vs2  cancels t^2 term - then solve for t 
crs = simplifyFraction(subs(vs1,t, vt))

[c,t] = coeffs(crs,[s1,s2])

simplifyFraction(c(2)^2-4*c(1)*c(4))
