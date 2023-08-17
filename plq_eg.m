function plq_eg
  x=sym('x');
  y=sym('y');
  f=functionF(x*y);
  
  


  
  

  
  d1=domain([0,0;2,0;2,1;1,1],x,y);
  d2=domain([-5,5;1,3;-1,0;-5,-4],x,y);
  d3=domain([-5,-4;0,0;1,1;1,3],x,y);
  %d3=domain([0,0;1,2;2,2;2,1],x,y);
  p(1)=plq_1piece(d1,f);
  p(2)=plq_1piece(d2,f);
  p(3)=plq_1piece(d3,f);
  

  PS = plq(p);
  %return
  
  convex_PS = PS.convexEnvelope();
  disp("convex done")
  convex_PS.plot;
%  return
  
  convex_PS.print
  return
  conjugate_PS = convex_PS.conjugate();
%   disp("conjugate done")
%  conjugate_PS.plot;
%  return
%  conjugate_PS.print
% return
%disp ("Check intersection")
  %[f2, r2] = conjugate_PS.intersectionConjugateDomain;
  mconjugate_PS = conjugate_PS.intersectionConjugateDomain;
  %r2(5).print
 % return
  disp("Start maximum")
  %[maxf,maxd] = mconjugate_PS.maximum %(f2,r2);
  mconjugate_PS = mconjugate_PS.maximum %(f2,r2);
  mconjugate_PS.print;
  return
%  for i = 1:size(maxf,2)
%      disp(i)
%      maxf(i).print
%      maxd(i).print
%  end
%  maxd.plotL
%return  
% merge wont work due to union vs intersection
  %[maxf,maxd] = conjugate_PS.merge(maxf,maxd);
 % return

  for i = 1:size(mconjugate_PS.pieces(1).maxf,2)
      disp(i)
      mconjugate_PS.pieces(1).maxf(i).print
      mconjugate_PS.pieces(1).maxd(i).print
  end
  %maxd.plotL

  
end
