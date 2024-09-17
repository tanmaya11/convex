% conjugate F*(s) = sup(s.x-f(x))  restricted to edge
% edge : edge
% f : f
% make appropriate substitutions to clean this

function conj = conjugateExpr(edge,f,x,y)
    edge
    f

    lambda = sym('lambda');
    infs = sym('inf');

    s1 = sym('s1');
    s2 = sym('s2');

    %edge
    %f
    dedge = sym.empty();
    dedge(1) = diff(edge,x)
    dedge(2) = diff(edge,y)

    df = sym.empty();
    df(1) = diff(f,x)
    df(2) = diff(f,y)
    
    eq1 = s1 -  df(1) - lambda*dedge(1)
    eq2 = s2 -  df(2) - lambda*dedge(2)
    eq3 = simplifyFraction(eq1*dedge(2)-eq2*dedge(1))
    eq3 = simplify(eq3) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % edge
    % [co2, va2] = coeffs(simplifyFraction(edge),[x,y,s1,s2])
    %  [co, va] = coeffs(simplifyFraction(eq3),[x,y,s1,s2])
    % 
    %  co4 = co2 * co(2)
    %  co3 = co * co2(2)
    % 
    %  simplifyFraction(co3(1)-co4(1))
    %  simplifyFraction(co3(2)-co4(2))
    %  simplifyFraction(co3(6)-co4(4))
    %  return
    % a1 = sym('a1');
    % b1 = sym('b1');
    % c1 = sym('c1');
    % d1 = sym('d1');
    % e1 = sym('e1');
    % f1 = sym('f1');
    % g1 = sym('g1');
    % h1 = sym('h1');
    % i1 = sym('i1');
    % j1 = sym('j1');
    % 
    % eq3 = a1*x^2 + b1*x*y + c1*y^2 + d1*x+ e1*y + f1 + g1*s1*x + h1*s2*y + i1*s1 + j1*s2;

    xyl = solve([eq1,eq2,edge],[x,y,lambda])
    size(xyl.x)
    if size(xyl.x,1) > 1
        xyl.x = xyl.x(2)
        xyl.y = xyl.y(2)
        xyl.lambda = xyl.lambda(2)
    end

    simplifyFraction(xyl.x)
    simplifyFraction(xyl.y)
    xyl.lambda

   
    disp('df subs')
    % simplifyFraction(subs(df(1),[x,y],[xyl.x,xyl.y]))
    % simplifyFraction(subs(df(2),[x,y],[xyl.x,xyl.y]))
    % simplifyFraction(subs(dedge(1),[x,y],[xyl.x,xyl.y]))
    % simplifyFraction(subs(dedge(2),[x,y],[xyl.x,xyl.y]))
    % simplifyFraction(s1*xyl.x + s2*xyl.y)
    % simplifyFraction(subs(f,[x,y],[xyl.x,xyl.y]))
    % simplifyFraction(s1*xyl.x + s2*xyl.y-subs(f,[x,y],[xyl.x,xyl.y]))
    % 
    % simplifyFraction(subs(df(1),[x,y],[xyl.x,xyl.y])- xyl.lambda*subs(dedge(1),[x,y],[xyl.x,xyl.y]))
    % simplifyFraction(subs(edge,[x,y],[xyl.x,xyl.y]))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % xy = solve([eq3,edge],[x,y]);
     % % size(xy.x)
     % % 
     % xy.x 
     % % 
     % xy.y

    % simplifyFraction(xy.x(1)-xy.x(2))
    % return
   
    conj = infs;
    
    if isempty(xyl)
        return
    elseif isempty(xyl.x) | isempty(xyl.y) | isempty(xyl.lambda)
        return
    end
    
    %conj = s1*xy.x + s2*xy.y - subs(f,[x,y],[xy.x,xy.y]);
    conj = s1*xyl.x + s2*xyl.y - subs(f,[x,y],[xyl.x,xyl.y]);
    conj = simplifyFraction(conj);
    %[co, va] = coeffs(conj,[s1,s2])
    conj = subs(conj,[s1,s2],[x,y])
    
    
end