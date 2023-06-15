function plq_eg
  x=sym('x');
  y=sym('y');
  f=functionF(x*y);

  
  

  
  d1=domain([0,0;2,0;2,1;1,1],x,y);
  d2=domain([-5,5;1,3;-1,0;-5,-4],x,y);
  d3=domain([0,0;1,2;2,2;2,1],x,y);

  p(1)=plq_1piece(d1,f);
  p(2)=plq_1piece(d2,f);
  p(3)=plq_1piece(d3,f);


  PS = plq(p);
  convex_PS = PS.convexEnvelope();
  
  convex_PS.print
  
  conjugate_PS = convex_PS.conjugate();
  return
  max_PS = conjugate_PS.maximum();
  
  
end
