% quad - quad verification

mh = sym('mh')
mw = sym('mw')
qh = sym('qh')
qw = sym('qw')
a = sym('a')
b = sym('b')

alpha0 = simplifyFraction((qw*mh-qh*mw+sqrt(mh*mw)*(qh-qw))/(mh-mw)    )
%alpha0 = subs(alpha0,[mh,mw],[1,1/2])
%alpha0 = subs(alpha0,[qh,qw],[-1,0])

alpha1 = -sqrt(mh*mw)
%alpha1 = subs(alpha1,[mh,mw],[1,1/2])

obj = -(a+mh*b-qh)^2/(4*mh) - b*mh+a*x+b*y
obj = subs(obj,a,alpha1*b+alpha0)

[c,t] = coeffs(obj,b)

psi0 = c(3)
psi1 = c(2)/2
psi2 = -c(1)

f0 = simplifyFraction(psi1^2/psi2+psi0)

% verifying grad = 0
bv = simplifyFraction(psi1/psi2)
simplifyFraction(-2*psi2*bv+2*psi1)
%%%%%%%

%return
% Example 1

% mh = 1
% mw = 1/2
% qh = -1 
% qw = 0


% psi1 formula doesnt verify
% psi2 = (alpha1+mh)^2/(4*mh)
% psi1 = -((alpha1+mh)*(alpha0-qh)/(2*mh)) + y - qh + alpha1*x
% psi0 = -(alpha0-qh)^2/(4*mh)+alpha0*x

% simplifyFraction (c(1)+psi2)
% simplifyFraction (c(2)-2*psi1)
% simplifyFraction (c(3)-psi0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% check bounds
f = x*y
f = subs (f,y,mh*x+qh)
dfdx = diff(f,x)
xv1 = sym('xv1')
subs(dfdx,x,xv1)
disp('bound')
subs(dfdx,x,1)
%return

boundEq = a + mh* b - subs(dfdx,x,xv1)
boundEq = subs(boundEq, a, alpha1*b+alpha0)
lb = solve(boundEq,b)

ineql = lb * psi2 - psi1    % >= 0
% when ineq <= 0
% obj -psi2*lb^2 + 2*lb*psi1+psi0


boundEq = a + mh* b - subs(dfdx,x,xv2)
boundEq = subs(boundEq, a, alpha1*b+alpha0)
ub = solve(boundEq,b)

inequ = ub * psi2 - psi1    % >= 0
% when ineq <= 0
% obj -psi2*ub^2 + 2*ub*psi1+psi0



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


simplifyFraction(c(2))
simplifyFraction(simplifyFraction(2*psi1))

%return

psi1 = c(2)/2




f = -psi2*bv^2+2*bv*psi1+psi0
f = simplifyFraction(subs(f,[mh,mw,qh,qw],[1,1/2,-1,0]))

double(2- subs(f,[x,y],[2,1]))
double(0- subs(f,[x,y],[1,0]))
double(0- subs(f,[x,y],[0,0]))
return
% 
% f1 = 5*x - 4*y - 3*2^(1/2)*x + 4*2^(1/2)*y + 3*2^(1/2) - 4
% 
% disp('f-f0')
% simplifyFraction(f-f0)
% 
% %return
% 
% f = subs(f,[mh,mw],[1,1/2])
% f = subs(f,[qh,qw],[-1,0])
% f = simplifyFraction(f)
% disp('f-f1')
% simplifyFraction(f-f1)
% 
% double(simplifyFraction(subs(f,[x,y],[0,0])))
% double(simplifyFraction(subs(f,[x,y],[1,0])))
% double(simplifyFraction(subs(f,[x,y],[2,1]))) -2
% double(simplifyFraction(subs(f,[x,y],[0.146447,0])))
% double(simplifyFraction(subs(f,[x,y],[1.146447,0.146447]))) - 1.146447*0.146447
% double(simplifyFraction(subs(f,[x,y],[0.792893,0.396447]))) - 0.792893*0.396447
% double(simplifyFraction(subs(f,[x,y],[0.085786,0.042893]))) - 0.085786*0.042893
% double(simplifyFraction(subs(f,[x,y],[0.3964,0])))

% check bounds
f = x*y
f = subs (f,y,mh*x+qh)
dfdx = diff(f,x)
xv1 = sym('xv1')
subs(dfdx,x,xv1)
disp('bound')
subs(dfdx,x,1)
%return
b
boundEq = a + mh* b - subs(dfdx,x,xv1)
boundEq = subs(boundEq, a, alpha1*b+alpha0)
lb = solve(boundEq,b)

ineq = lb * psi2 - psi1
ineq = subs(ineq,[mh,mw],[1,1/2])
ineq = subs(ineq,[qh,qw],[-1,0])
ineq = subs(ineq,[xv1],[0])
disp('ineq1')
ineq = simplifyFraction(ineq)



ineq = lb * psi2 - psi1
ineq = subs(ineq,[mh,mw],[1,1/2])
ineq = subs(ineq,[qh,qw],[-1,0])
ineq = subs(ineq,[xv1],[2])
disp('ineq2')
ineq = simplifyFraction(ineq)

lb = subs(lb,[xv1,mh,mw,qh,qw],[2,1,1/2,-1,0])
double(lb)
f2 = -psi2*lb^2 + 2*lb*psi1+psi0
f2 = simplifyFraction(subs(f2,[mh,mw,qh,qw],[1,1/2,-1,0]))
double(subs(f2,[x,y],[0,0]))