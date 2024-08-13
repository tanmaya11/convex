function [conj] = conjugateExpr(q,eq,x,y)

    lambda = sym('lambda');
    infs = sym('inf');

    s1 = sym('s1');
    s2 = sym('s2');

    del = sym.empty();
    del(1) = diff(q,x);
    del(2) = diff(q,y);
    del1 = sym.empty();
    del1(1) = diff(eq,x);
    del1(2) = diff(eq,y);
    eq1 = s1 -  del1(1) - lambda*del(1);
    eq2 = s2 -  del1(2) - lambda*del(2);
    eq3 = eq1*del(2)-eq2*del(1);
    
    xy = solve([eq3,q],[x,y])
    conj = infs;
    
    if isempty(xy)
        return
    elseif isempty(xy.x) | isempty(xy.y)
        return
    end
    
    conj = s1*xy.x + s2*xy.y - subs(eq,[x,y],[xy.x,xy.y]);
    conj = simplifyFraction(conj);
    %eq1 = subs(eq1,[x,y],[xy.x,xy.y]);
    %l = solve(eq1,lambda);
    
end