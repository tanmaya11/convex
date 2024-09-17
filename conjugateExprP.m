% conjugate F*(s) = sup(s.x-f(x))  restricted to edge
% edge : edge
% f : f
% make appropriate substitutions to clean this

function conj = conjugateExpr(edge,f,x,y)

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


    
    %edge
     [co, va] = coeffs(simplifyFraction(eq3),[x,y])
    % 
    % a1 = simplifyFraction(co(2)^2-4*co(1)*co(3))
    % subs(a1,b,sqrt(4*a*c))
    % 
    % [co2, va2] = coeffs(simplifyFraction(edge),[x,y])
    xy = solve([eq3,edge],[x,y]);
    % size(xy.x)
    % 
    %  xy.x 
    % 
    % xy.y
   
    conj = infs;
    
    if isempty(xy)
        return
    elseif isempty(xy.x) | isempty(xy.y)
        return
    end
    
    conj = s1*xy.x(1) + s2*xy.y(1) - subs(f,[x,y],[xy.x(1),xy.y(1)]);
    conj = simplifyFraction(conj);
    conj = subs(conj,[s1,s2],[x,y]);
    
    
end