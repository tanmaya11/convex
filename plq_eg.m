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
  %return
  
  convex_PS = PS.convexEnvelope();
  %return
  
%  convex_PS.print
%  return
  conjugate_PS = convex_PS.conjugate();
%  conjugate_PS.print
% return
%disp ("Check intersection")
  [f2, r2] = conjugate_PS.intersectionConjugateDomain;
  %r2(5).print
%  return
  disp("Start maximum")
  [maxf,maxd] = conjugate_PS.maximum(f2,r2);
  
  for i = 1:size(maxf,2)
      disp(i)
      maxf(i).print
      maxd(i).print
  end
  
  
end
