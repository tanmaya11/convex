function bivariate_eg2

  % f = xy defined on {(0,0),(2,0),(2,1),(1.1)}
  x=sym('x');
  y=sym('y');
  f=functionF(x*y);

  d=domain([-5,5;1,3;-1,0;-5,-4],x,y);
  
  p=plq_1piece(d,f);

  p=p.convexEnvelope;

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


