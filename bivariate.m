% f = xy defined on {(0,0),(2,0),(2,1),(1.1)}
x=sym('x');
y=sym('y');
f=functionF(x*y);

d=domain([0,0;2,0;2,1;1,1]);
p=plq_1piece(d,f);

p=p.convexEnvelope;
%return

for i=1:size(p.envf,2) 
  p.envf(i)
  p.envd(i)
end

%p = p.maxEnvelope(x,y)
