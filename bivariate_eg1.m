function bivariate_eg1 

  % f = xy defined on {(0,0),(2,0),(2,1),(1.1)}
  x=sym('x');
  y=sym('y');
  f=functionF(x*y);

  d=domain([0,0;2,0;2,1;1,1],x,y);
  
  p=plq_1piece(d,f);

  p=p.convexEnvelope;
  disp("envelope")
  for i=1:size(p.envf,2) 
    disp('Function')  
    p.envf(i).print
    disp('Domain')
    p.envd(i).print
  end
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


