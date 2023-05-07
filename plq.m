classdef plq
    properties
        f=functionF.empty();
        d=domain.empty();
        envf=functionF.empty();
        envd=domain.empty();
    end

    methods
        function obj = plq(d,f)
            if nargin > 0
              obj.f = [obj.f,f];
              obj.d = [obj.d,d];
            end 
        end

        function convexEnvelope1 (obj,i1,x,y)
            a=sym('a');
            b=sym('b');
            obj.d(i1).V    
            [etaV, etaE, etaR] =  getEtaFunctions (obj,i1,x,y,a,b)
            % put a check that eta are only polynomials 

             % (vix, vjx) is the pair :  1 for edge else 0 
            [ix,jx,vix, vjx] = feasiblePairs (obj, i1, etaR,a,b)
            solve (obj, ix,jx,vix, vjx,etaV, etaE, etaR,a, b, x, y);
        end 

        function solve (obj, ix,jx,vix, vjx,etaV, etaE, etaR, a, b, x, y)
            for i=1:size(ix,2)
                disp("i")
                disp(i)
                nc = 0;
                    
                if (vix(i)==1)
                    etah = etaE(ix(i));
                    degreeh = polynomialDegree(etah);
                    nc = nc + 1;
                    c(nc) = etaR(ix(i))
                else
                    etah = etaV(ix(i));
                    degreeh = polynomialDegree(etah);
                end
                if (vjx(i)==1)
                    etaw = etaE(jx(i));
                    degreew = polynomialDegree(etaw);
                    nc = nc + 1;
                    c(nc) = etaR(jx(i))
                else
                    etaw = etaV(jx(i));
                    degreew = polynomialDegree(etaw);
                end
                disp("nc")
                disp(nc)
                disp(c(nc))
                disp("eta")
                etah
                etaw
                if (degreeh==1 & degreew==1)
                    disp("here")
                    obj = etah + a*x+b*y
                    bv = solve(etah==etaw,a);
                    obj = subs(obj,a, bv)
                    for k = 1:nc
                      c(k) = subs(c(k),a, bv);
                    end
                    disp("etaV")
                    size(etaV,2)
                    etaV
                    for j=1:size(etaV,2)
                        j
                        etaV(j)
                        if ((etaV(j) == etah)  )
                            continue;
                        end
                        if ( (etaV(j) == etaw) )
                            continue;
                        end
                        etak = etaV(j)
                        nc = nc+1;
                        disp("c")
                        c(nc) = etah - etak  % <= 0   easier for substitution
                        c(nc) = subs(c(nc),a, bv)
                        
                    end

                    % put etaE loop
                    disp("here end")
                    
                end 
                for k = 1:nc
                    disp("c")
                    c(k)
                end
                    
                        
            end     
        end

        
        function [ix,jx,vix, vjx] = feasiblePairs (obj, i1, etaR,a,b)
            n = 0;
            for i = 1:size(etaR,2)
                for j = 1:size(obj.d(i1).V,2)
                    vertex = [obj.d(i1).getVertex(obj.d(i1).V(j))];
                    s = subs(etaR(i),[a,b],vertex);
                    
                    % temperory fix
                    % piecewise stored as inequality , others as lhs, where
                    % lhs <= 0 

                    if (symType(s)=="equation")
                    
                        if (s )
                          n = n + 1;
                          ix(n) = i;
                          jx(n) = j;
                          vix(n) = 1;
                          vjx(n) = 0;
                        end
                    
                    else
                        if (s <= 0)
                          n = n + 1;
                          ix(n) = i;
                          jx(n) = j;
                          vix(n) = 1;
                          vjx(n) = 0;
                        end
                    
                    end
                    
                end
            end
            for i = 1:size(etaR,2)
                for j = i+1:size(etaR)
                    if (logical(eval(etaR(i) + etaR(j))))
                        
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 1;
                    end
                end
            end
            
        end 
        
        function [etaV, etaE, etaR] = getEtaFunctions (obj,i1,x,y,a,b)
            f=sym('f');
            eta=sym('eta');
            m=sym('m');
            c=sym('c');
            f2 = functionF(a*x+b*y);
            eta = obj.f(i1) - f2;
            
            edgey = m * x + c;
            
            % Eta for Edges
            for i = 1:obj.d(i1).nE
                xv1 = obj.d(i1).vx(obj.d(i1).E(i,1));
                yv1 = obj.d(i1).vy(obj.d(i1).E(i,1));
                

                fT = obj.f(i1).subsF(x,y,xv1,yv1);
                etaT = subs(eta,[x,y,f],[xv1,yv1,fT]);
                
                etaE(3*(i-1)+1) = etaT;

                


                xv2 = obj.d(i1).vx(obj.d(i1).E(i,2))
                yv2 = obj.d(i1).vy(obj.d(i1).E(i,2))
                
                fT = obj.f(i1).subsF(x,y,xv2,yv2);
                etaT = subs(eta,[x,y,f],[xv2,yv2,fT]);
                
                etaE(3*(i-1)+3) = etaT;



                m = obj.d(i1).mE(i);
                c = obj.d(i1).cE(i);
                yT = subs(edgey);
                
                etaT = subs(eta,y,yT);

                %dfdx = diff (etaT,x)
                xp = solve(diff (etaT,x),x);

                etaM = simplify(subs(etaT,x,xp));
                etaE(3*(i-1)+2) = etaM;
                
                
                f1 = obj.f(i1).subsF(x,y,x,yT);
                df1 = diff(f1,x);
                % store only <= in order to simplify later  (>= better but currently no multiplications involved)
                etaR(3*(i-1)+1) = a + m*b - subs(df1,[x,y],[xv1,yv1]); % <= 0
                %etaR(3*(i-1)+2) =  subs(df1,[x,y],[xv1,yv1]) <= a + m*b <= subs(df1,[x,y],[xv2,yv2]);
                etaR(3*(i-1)+2) =  piecewise( -a - m*b  + subs(df1,[x,y],[xv1,yv1])<=0,  (a + m*b - subs(df1,[x,y],[xv2,yv2])<=0));
                %etaR(3*(i-1)+2) =  piecewise( -a - m*b  + subs(df1,[x,y],[xv1,yv1]),  (a + m*b - subs(df1,[x,y],[xv2,yv2])));
                etaR(3*(i-1)+3) = -a - m*b + subs(df1,[x,y],[xv2,yv2]); % <=0

            end

            % Eta for Vertices
            for i = 1:obj.d(i1).nV
                xv = obj.d(i1).vx(obj.d(i1).V(i));
                yv = obj.d(i1).vy(obj.d(i1).V(i));
                fT = obj.f(i1).subsF(x,y,xv,yv);
                f=fT;
                x=xv;
                y=yv;
                etaT = subs(eta);
                etaV(i) = etaT; 
            end
            
        end
    end
end