x=sym('x');
y=sym('y');
f=functionF(x*y);
d=domain([0,0;2,0;2,1;1,1]);
d = d.getEdges;
p=plq(d,f);

p = p.convexEnvelope1(1,x,y);


for i=1:size(p.envf,2) 
  p.envf(i)
  p.envd(i)
end

p = p.maxEnvelope(x,y)
