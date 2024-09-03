% conjugate F*(s) = sup(s.x-f(x))  restricted to edge
% edge : edge
% f : f
% make appropriate substitutions to clean this

function [conj] = conjugateExpr(edge,f,x,y)

    lambda = sym('lambda');
    infs = sym('inf');

    s1 = sym('s1');
    s2 = sym('s2');

    dedge = sym.empty();
    dedge(1) = diff(edge,x)
    dedge(2) = diff(edge,y)

    df = sym.empty();
    df(1) = diff(f,x)
    df(2) = diff(f,y)
    
    eq1 = s1 -  df(1) - lambda*dedge(1)
    eq2 = s2 -  df(2) - lambda*dedge(2)
    eq3 = eq1*dedge(2)-eq2*dedge(1)
    edge
    
    xy = solve([eq3,edge],[x,y])
    conj = infs;
    
    if isempty(xy)
        return
    elseif isempty(xy.x) | isempty(xy.y)
        return
    end
    
    conj = s1*xy.x + s2*xy.y - subs(f,[x,y],[xy.x,xy.y]);
    conj = simplifyFraction(conj);
    conj = subs(conj,[s1,s2],[x,y]);
    
    
end