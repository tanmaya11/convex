function bivariate_eg3

  % f = xy defined on {(0,0),(2,0),(2,1),(1.1)}
  x=sym('x');
  y=sym('y');
  f=functionF(x*y);

  d=domain([0,0;1,2;2,2;2,1],x,y);
  % E
  % 1     2
  % 4     1
  % mE
  % 2.0000    0.5000
  % V 
  % 3

  %ineqs
  %y - 2*x
  %y - 2
  %x - 2
  %x/2 - y

  
  p=plq_1piece(d,f);
  
 
  p=p.convexEnvelope;
  %etaV
  %4-2a-2b

  %etaE
  % 0    -(a/2 + b)*(a/4 + b/2)   2 - 2*b - a
  % 0    -(2*a + b)^2/8           2 - b - 2*a
 
  %etaR
  %a + 2*b  0  5
  %a + b/2  0  5/2

  %feasible pairs
  % (ix,jx) feasible pair
  % (vix,vjx) 0 : from V, 1 : from E
  % ixd,jxd : if from E - gives region no - 1,2,3
  % ix =     1     2     1     1     1     1     1     1     1     1     1
  % jx =     1     1     2     2     2     2     2     2     2     2     2
  %vix =     1     1     1     1     1     1     1     1     1     1     1
  %vjx =     0     0     1     1     1     1     1     1     1     1     1
  %ixd =     3     3     1     1     1     2     2     2     3     3     3
  %jxd =     0     0     1     2     3     1     2     3     1     2     3

  %solve
  %i = 1

  %etah = f: 2 - 2*b - a
  %mh = 2
  %qh = 0
  %etaw = f: 4 - 2*b - 2*a
  % lin-lin
  %av = 2
  %obj0 = f: a*x - 2*b - a + b*y + 2
  %obj0 = f: 2*x - 2*b + b*y
  %objfacts = [2*x, y - 2]
  %lb =[1.5,0,2]
  %ub =[Inf,Inf,Inf]
  %mlb = 2
  %mub = Inf
  %envfs0  2*x + 2*y - 4

  %i = 2
  %etah = f: 2 - b - 2*a
  %mh = 0.5000
  %qh = 0
  %etaw = f: 4 - 2*b - 2*a
  %lin-lin
  %f0 = f: b - 2
  %a = a
  %av = Empty sym: 0-by-1
  %f0 = f: b - 2
  %a = b
  %av = 2
  %obj0 = f: a*x - b - 2*a + b*y + 2
  %obj0 = f: 2*y - 2*a + a*x
  %objfacts = [2*y, x - 2]
  %lb = 1.5000         0    4.0000
  %ub = Inf   Inf     4
  %mlb =    4
  %mub =    4
  %4*x + 2*y - 8

  % i = 3
  % etah =    f: 0
  % etaw =    f: 0

  %i = 4
  %etah = f: 0
  %mh = 2
  %qh =0
  %etaw = f: -(2*a + b)^2/8
  %mw = 0.5000
  %qw = 0

  %i = 5
  %etah =f: 0
  %etaw = f: 2 - b - 2*a

  %i= 6
  %etaw = 0

  %i = 7
  %etah = f: -(a/2 + b)*(a/4 + b/2)
  %mh = 2
  %qh = 0
  %etaw = f: -(2*a + b)^2/8
  %mw = 0.5000
  %qw = 0
  % quad-quad
  %mh /= mw
  %alpha1 =  1
  %alpha0 = 0
  %lb = [0, -5]
  %ub = [5, 0]


  for i=1:size(p.envf,2) 
    disp('Function')  
    p.envf(i).print
    disp('Domain')
    p.envd(i).print
  end
return
  p = p.maxEnvelope;
  for i=1:size(p.envf,2) 
    disp('Function')  
    p.envf(i).print
    disp('Domain')
    p.envd(i).print
  end


end


