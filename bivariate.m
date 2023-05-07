x=sym('x');
y=sym('y');
f=functionF(x*y);
%f2=functionF(x*y);
%f3 = f-f2
d=domain([0,0;2,0;2,1;1,1]);
d = d.getEdges;
p=plq(d,f);

p.convexEnvelope1(1,x,y)
