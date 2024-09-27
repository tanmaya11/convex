% quad - linear case of stp 1 
% rational function with 0/0 at vertex

a= sym('a')
b= sym('b')
m= sym('m')
q= sym('q')
fv1 = sym('fv1')
xv1 = sym('x_1')
yv1 = sym('y_1')

fv2 = sym('fv2')
xv2 = sym('x_2')
yv2 = sym('y_2')

fv3 = sym('fv3')
xv3 = sym('x_3')
yv3 = sym('y_3')


eta(1) = fv1 - a*xv1 - b*yv1
eta(2) = fv2 - a*xv2 - b*yv2
eta(3) = fv3 - a*xv3 - b*yv3

n = 0

for i = 1:3
    obj = eta(i) + (a*x + b*y)
    for j = i+1:3
        eq = eta(i) - eta(j)

        av = solve(eq,a);

        obj0 = subs(obj, [a],[av]);
        objfacts = coeffs(obj0,b);
        if size(objfacts,2) == 1
          eta0 = 0;
          eta1 = objfacts(1);  
    
        else
          eta0 = objfacts(1);
          eta1 = objfacts(2);  
        end

    

        for k=1:3
            if k == i
                continue
            end
            if k == j
                continue
            end
                   
            c = eta(i) - eta(k)
            c =  simplifyFraction(subs(c,[a],[av]))
            s = coeffs(c,b)
            bv = simplifyFraction(solve(c,b)         )
            lb= -inf
            ub = bv
            obj0 = subs(obj,a,av)
            obj0 = simplifyFraction(subs(obj0,b,bv))

            obj0 = simplifyFraction(subs(obj0,fv1,xv1*yv1))
            obj0 = simplifyFraction(subs(obj0,fv2,xv2*yv2))
            obj0 = simplifyFraction(subs(obj0,fv3,xv3*yv3))
            n = n + 1
            f(n) = obj0
            coeffs(obj0,[x,y])
        end
        
   end
            
end

f(1)
f(2)
f(3)



f(1) - f(2)
f(1) - f(2)
f(2) - f(3)

[c,t] = coeffs(f(1),[x,y])
f(1)
[num,den] = numden(f(1))
c = den*c
c(1)
c(2)
c(3)

subs(c(1),[xv1,xv2,xv3,yv1,yv2,yv3],[1,-1,-1,-1,-1,1])/subs(den,[xv1,xv2,xv3,yv1,yv2,yv3],[1,-1,-1,-1,-1,1])
subs(c(2),[xv1,xv2,xv3,yv1,yv2,yv3],[1,-1,-1,-1,-1,1])/subs(den,[xv1,xv2,xv3,yv1,yv2,yv3],[1,-1,-1,-1,-1,1])
subs(c(3),[xv1,xv2,xv3,yv1,yv2,yv3],[1,-1,-1,-1,-1,1])/subs(den,[xv1,xv2,xv3,yv1,yv2,yv3],[1,-1,-1,-1,-1,1])


subs(c(1),[xv1,xv2,xv3,yv1,yv2,yv3],[0,2,1,2,0,2])/subs(den,[xv1,xv2,xv3,yv1,yv2,yv3],[0,2,1,2,0,2])
subs(c(2),[xv1,xv2,xv3,yv1,yv2,yv3],[0,2,1,2,0,2])/subs(den,[xv1,xv2,xv3,yv1,yv2,yv3],[0,2,1,2,0,2])
subs(c(3),[xv1,xv2,xv3,yv1,yv2,yv3],[0,2,1,2,0,2])/subs(den,[xv1,xv2,xv3,yv1,yv2,yv3],[0,2,1,2,0,2])




































