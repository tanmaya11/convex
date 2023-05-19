function bivariate_eg2

  % f = xy defined on {(0,0),(2,0),(2,1),(1.1)}
  x=sym('x');
  y=sym('y');
  f=functionF(x*y);

  d=domain([-5,5;1,3;-1,0;-5,-4],x,y);
  %return
  % nVertices =   4
  %d.vx
  %  -5
  %   1
  %  -1
  %  -5

  %d.vy
  %   5
  %   3
  %   0
  %  -4

 % E =
 %   2     3
 %   3     4

 % m =  1.5000    1.0000
 % V = 1
 % q =  1.5000    1.0000

 %x/3 + y - 10/3
 %(3*x)/2 - y + 3/2
 %x - y + 1
 %- x - 5
  
  
  p=plq_1piece(d,f);

  p=p.convexEnvelope;

  %etaV
  %5*a - 5*b - 25

  %etaE
  % a
  % a/2 - (3*b)/4 - (a*b)/2 - a^2/6 - (3*b^2)/8 - 3/8
  % 3 - 3*b - a
  % 5*a + 4*b + 20
  % - a^2/4 - (a*b)/2 + a/2 - b^2/4 - b/2 - 1/4
  % a

  %etaR
  %a + (3*b)/2
  %-3/2
  %9/2
  %a + b
  %-9
  %-1
 
  %feasiblePairs
  %ix =   1     2     1     1     1     1     1     1     1     1     1
  %jx =   1     1     2     2     2     2     2     2     2     2     2
  %vix =      1     1     1     1     1     1     1     1     1     1     1
  %vjx =      0     0     1     1     1     1     1     1     1     1     1
  %ixd =     2     3     1     1     1     2     2     2     3     3     3
  %jxd =     0     0     1     2     3     1     2     3     1     2     3
  disp("Envelope")
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


