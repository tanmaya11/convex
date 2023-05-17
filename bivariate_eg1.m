function bivariate_eg1 

  % f = xy defined on {(0,0),(2,0),(2,1),(1.1)}
  x=sym('x');
  y=sym('y');
  f=functionF(x*y);

  d=domain([0,0;2,0;2,1;1,1],x,y);

  % d.nVertices
  % 4
  % d.vx
  %   1
  %   2
  %   2
  %   0
  % d.vy
  %   1
  %   1
  %   0
  %   0

  %for i = 1:size(d.ineqs,2)
  %    d.ineqs(i).print
  %end
  %y - 1
  %x - 2
  %-y
  %y - x

  %d.E
  % 4     1

  %d.mE
  %     1

  %d.cE
  %  0

  %d.V
  %2     3

  p=plq_1piece(d,f);
  %p.f.print
  %x*y
 
  %p.d
  % domain with properties:

  %  nVertices: 4
  %         vx: [4×1 double]
  %         vy: [4×1 double]
  %      ineqs: [1×4 functionF]
  %         nE: 1
  %          E: [4 1]
  %         mE: 1
  %         cE: 0
  %         nV: 2
  %          V: [2 3]
  
  p=p.convexEnvelope;
  disp("out")
  %getEtaFunctions
  %Printing list
  %etaV
  % 2 - b - 2*a
  % -2*a
  %Printing list
  %etaE
  % 0
  % -(a + b)^2/4
  % 1 - b - a
  %Printing list
  %etaR
  %a + b
  %0
  %2

  %feasible pairs
  %ix = 1     1     1
  %jx = 1     2     2
  %vix = 1     1     1
  %vjx = 0     0     0
  %ixd = 3     3     2
  %jxd = 0     0     0


  %solve
  % i = 1
  % yintercept = 0
  % m = 1
  % c = f: 2 - b - a
  % etah 
  % f: 1 - b - a
  % etaw = 
  % f: 2 - b - 2*a
  % degreeh = 1
  % degreew = 1

  %i = 2
  %yintercept = 0
  %m = 1
  %c = f: 2 - b - a
  % etah = 
  %f: 1 - b - a
  %etaw = 
  %f: -2*a
  %degreeh = 1
  %degreew = 1

  %i = 3
  %yintercept = 0
  %m = 1
  %c = f: 2 - b - a
  %c = a + b - 2
  
   
  %etah = f: -(a + b)^2/4
  %etaw = f: -2*a
  %degreeh = 2
  %degreew = 1
  
  % tests linear-linear
  % i=1, j=1
  % av = 1
  %obj0 = a*x - b - a + b*y + 1
  %     = x - b + b*y
  %objfacts = [x, y - 1]
  % f : 1-b    from etaR
  %lb 1
  %ub Inf
  %etaK
  % f :  a - b + 1
  %  2 - b
  %lb 2
  %ub Inf
  % lb 1 2
  % ub Inf Inf
  %mlb = 2
  %mub=Inf
  % b = lb
  % x + 2*y - 2
  % y - 1

  % tests linear-linear
  % i=1
  % j=2
  % obj0 = b*y - 2*b + x*(b - 1) + 2
  % objFacts = [2 - x, x + y - 2]
  % 3 - 2*b
  % lb = 1.5
  % etak 
  % a-1
  %f: b - 2
  %lb =     1.5000      -Inf
  %ub =     Inf     2
  %mlb = 1.5
  %mub = 2
  %x/2 + (3*y)/2 - 1
  %x + 2*y - 2
  % x + y - 2, 2 - y - x

  %i=1
  %j=2
  % tests quad-linear
  %eq : 2*a - (a + b)^2/4
  % av = z - b
  % bv = - z^2/8 + z
  % lb = 0
  % ub = 2
  % etak = f: 2 - b - 2*a
  % c = f: 2*a + b - (a + b)^2/4 - 2
  % -(z - 4)^2/8
  % check this
 
  % obj0 = f: a*x + b*y - (a + b)^2/4
  % obj0 = f: (x*z^2)/8 - z^2/4 + y*(- z^2/8 + z)
  % objfacts = [y, x/8 - y/8 - 1/4]
  % psi0 = 0
  % psi1 = y/2
  % psi2 = y/8 - x/8 + 1/4
  % f0 = (2*y^2)/(y - x + 2)
  % f0 = x/2 + (3*y)/2 - 1
  
 
  
  disp("envelope")
  size(p.envf)
  for i=1:size(p.envf,2) 
    disp('Function')  
    p.envf(i).print
    disp('Domain')
    p.envd(i).print

  end
  return
  %
  %f = functionF(x-1);
  %g = functionF(x-2);
  %isSubset(f,g)
  %return
disp("b4")
  p = p.maxEnvelope;
  disp("Max")
  for i=1:size(p.envf,2) 
    disp('Function')  
    p.envf(i).print
    disp('Domain')
    p.envd(i).print
  end


end


