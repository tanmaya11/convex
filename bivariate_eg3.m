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

  for i=1:size(p.envf,2) 
    disp('Function')  
    p.envf(i).print
    disp('Domain')
    p.envd(i).print
  end

  p = p.maxEnvelope;
  for i=1:size(p.envf,2) 
    disp('Function')  
    p.envf(i).print
    disp('Domain')
    p.envd(i).print
  end


end


