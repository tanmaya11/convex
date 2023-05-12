classdef plq_1piece
    properties
        f;
        d;
        envf=functionF.empty();
        envd=functionF.empty();
    end

    methods
        function obj = plq_1piece(d,f)
            % put checks for type of f and d
            if nargin > 0
              obj.f = f;
              obj.d = d;
            end 
        end


        % fix for multiple intersections
        % remove outside polytope ineqs
        function obj = maxEnvelope (obj)
            envdT = [];
            envfT = [];
            for i = 1:size(obj.envd,2)
                l(i)=1;
            end
            for i = 1:size(obj.envd,2)
                for j = i+1:size(obj.envd,2)
                    if (obj.envd(i) == obj.envd(j))
                      envdT = [envdT,obj.envd(i)];
                      %pointwise_max(obj.envf(i), obj.envf(j), obj.d.vx, obj.d.vy, obj.d.ineqs, [obj.envd(j)])
                      envfT = [envfT, pointwise_max(obj.envf(i), obj.envf(j), obj.d.vx, obj.d.vy, obj.d.ineqs, [obj.envd(j)])]     ;
                      l(i) = 0;
                      l(j) = 0;
                      disp("in1")
                    
                    end
                end
              
              
            end
            for i = 1:size(obj.envd,2)
              if (l(i) == 1)
                    disp("in2")
                    envdT = [envdT,obj.envd(i)];
                    envfT = [envfT, obj.envf(i)];     
              end
            end
            envfT
            obj.envf = envfT;
            obj.envd = envdT;
            
        end

        function f = max(obj,i,j,x,y)
            for i = 1:obj.d.nVertices
                if (subs(obj.envd(i).f,[x,y],[obj.d.vx(i),obj.d.vy(i)]))    % change this to vertexindomain
                    f0 = obj.envf(i)-obj.envf(j)
                end
            end
            f = functionF
        end

        function obj = convexEnvelope(obj)
            for i = 1:size(obj.f,2)
              vars = obj.f(i).getVars;
              if (size(vars,2)==2)
                  x = vars(1);
                  y = vars(2);
              else
                  disp("not bivariate in 'plq.m'")
                  return
              end
              obj = convexEnvelope1 (obj,x,y);
            end
            return
            
        end

        function obj = convexEnvelope1 (obj,x,y)
            a=sym('a');
            b=sym('b');
              
            % etaV : eta functions corresponding to set obj.d.V
            % etaE : eta functions corresponding to set obj.d.E
            % etaR : domain of etaE - stored as etaR(i,1:3) : [function,lb,ub] => lb <= function <= ub 
            disp("getEtaFunctions")
            [etaV, etaE, etaR] =  getEtaFunctions (obj,x,y,a,b);

            %etaV.printL();
            %etaE.printL();
            %etaR.printL();
            
            % put a check that eta are only polynomials 

             % (vix, vjx) is the pair :  1 for edge else 0 
            disp("feasiblePairs") 
            [ix,jx,vix, vjx, ixd, jxd] = feasiblePairs (obj,etaR, a,b);
            disp("solve")
            [envfs, envds] = solve (obj, ix,jx,vix, vjx,ixd, jxd, etaV, etaE, etaR,a, b, x, y);
            for i = 1:size(envfs,2)
                obj.envf = [obj.envf, envfs(i)];
                obj.envd = [obj.envd, functionF(envds(i))];
            end
            for i = 1:size(obj.envd,2)
              obj.envd(i) = removeDenominator (obj.envd(i),x,y);
            end 
           
            % put code for max 
        end 

                 

        function [etah, m, nc, c] = edgeInfoInSolve(obj, etaE, etaR, ix, ixd, nc, c)
           etah = etaE(ix,ixd);
           if size(obj.d.mE,2) == 1
             m = obj.d.mE;
           else  
             m = obj.d.mE(ix);
           end
           if (ixd==1)
             nc = nc + 1;
             c(nc) = etaR(ix,1)-etaR(ix,3);
           end
           if (ixd==3)
             nc = nc + 1;
             c(nc) = etaR(ix,3)-etaR(ix,1);
           end
           if (ixd==2)
             nc = nc + 1;
             c(nc) = etaR(ix,3)-etaR(ix,1);
             nc = nc + 1;
             c(nc) = etaR(ix,1)-etaR(ix,3);
           end
           
        end
                
        function [envfs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, nc, c, x, y, a, b, envfs, envds) 
            f0 = etah-etaw;
            av = f0.solve(a);
            if (isempty(av))
                lSol = false;
                return;
            end
            lSol = true;
            obj0 = obj0.subsVarsPartial([a],[av]);
            objfacts = coeffs(obj0.f,b);
            for k = 1:nc
               c(k) = c(k).subsVarsPartial([a],[av])  ;
            end
            for j=1:size(etaV,2)
               if (lV(j)) 
                 continue;
               end
               etak = etaV(j);
               nc = nc+1;
               c(nc) = etah - etak;  % <= 0   easier for substitution
               c(nc) = c(nc).subsVarsPartial([a],[av]);
            end
            for j=1:size(etaE,1)
               if (lE(j))
                 continue;
               end
               for k = 1:size(etaE,2) 
                  etak = etaE(j,k);
                  nc = nc+1;
                  c(nc) = etah - etak;  % <= 0   easier for substitution
                  c(nc) = c(nc).subsVarsPartial([a],[av]);
                  
              end
            end
            nbl = 0;
            nbr = 0;
            nzr = 0;
            for j = 1:nc
                 bv = c(j).solve(b);
                 for k = 1:size(bv)
                    %c(j).subsVarsPartial([b],[bv(k)]);
                    %c(j).subsVarsPartial([b],[bv(k)]) == functionF(0)
                    if (c(j).subsVarsPartial([b],[bv(k)]) == functionF(0))
                      nzr = nzr+1;
                      bz(nzr)=bv(k);
                    elseif (c(j).subsVarsPartial([b],[bv(k)]) < functionF(0))
                      nbl = nbl+1;
                      bl(nbl)=bv(k);
                    else
                      nbr = nbr+1;
                      br(nbr)=bv(k);
                    end
                 end
            end
                if nbl > 0 & nbr > 0
                  blv = max(bl);
                  brv = min(br);
                  if (blv > brv)
                    disp('infeasible');
                  else
                    envfs = [envfs, obj0.subsVarsPartial([b],[blv])];
                    envds = [envds, objfacts(2) ];
                    envfs = [envfs, obj0.subsVarsPartial([b],[brv])];
                    envds = [envds, -objfacts(2) ];
                    disp('feasible')    
                  end
                elseif nbl > 0 
                  blv = max(bl);
                  brv = min(bz);
                  if (blv > brv)
                    disp('infeasible')
                  else
                    envfs = [envfs, obj0.subsVarsPartial([b],[blv])];
                    envds = [envds, objfacts(2) ];
                    envfs = [envfs, obj0.subsVarsPartial([b],[brv])];
                    envds = [envds, -objfacts(2) ];
                    disp('feasible')    
                  end
                elseif nbr > 0 
                  blv = max(bz);
                  brv = min(br);
                  if (blv > brv)
                    disp('infeasible')
                  else
                    envfs = [envfs, obj0.subsVarsPartial([b],[blv])];
                    envds = [envds, objfacts(2) ];
                    envfs = [envfs, obj0.subsVarsPartial([b],[brv])];
                    envds = [envds, -objfacts(2) ];
                    disp('feasible')    
                  end
                else  
                  blv = min(bz);
                  brv = max(bz);
                  if (blv > brv)
                    disp('infeasible')
                  else
                    envfs = [envfs, obj0.subsVarsPartial([b],[blv])];
                    envds = [envds, objfacts(2) ];
                    envfs = [envfs, obj0.subsVarsPartial([b],[brv])];
                    envds = [envds, -objfacts(2) ];
                    disp('feasible')    
                  end
                end
        
        end

        function [envfs, envds] = solve (obj, ix,jx,vix, vjx, ixd, jxd, etaV, etaE, etaR, a, b, x, y)
            envfs = [];

            envds = [];
            disp("solve")
            
            for i=1:size(ix,2)

              for j = 1:size(etaV,2)
                lV(j) = false;
              end
              for j = 1:size(etaE,1)
                lE(j) = false;
              end
            
                
                nc = 0;
                if (vix(i)==1)
                    c(1)=etaR(1,1);
                    [etah, mh, nc, c] = edgeInfoInSolve(obj, etaE, etaR, ix(i), ixd(i), nc, c);
                    lE(ix(i)) = true;
                else
                    etah = etaV(ix(i));
                    lV(ix(i)) = true;
                end
                if (vjx(i)==1)
                    [etaw, mw, nc, c] = edgeInfoInSolve(obj, etaE, etaR, jx(i), jxd(i), nc, c);
                    lE(jx(i)) = true;
                else
                    etaw = etaV(jx(i));
                    lV(jx(i)) = true;
                end
                etah
                etaw
                degreeh = polynomialDegree(etah.f)
                degreew = polynomialDegree(etaw.f)    
                if (degreeh==1 & degreew==1)
                    obj0 = etah + functionF(a*x+b*y);
            
                    [envfs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, nc, c, x, y, a, b, envfs, envds) 
                    if (~lSol)
                        disp("lSol")
                        [envfs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, nc, c, x, y, b, a, envfs, envds) 
                    end
                end 
                
                if (degreeh==2 & degreew==1)
                    disp("one quad")
                    z = sym('z');
                    av = z - mh*b;
                    eq = etah-etaw;
                    eq = eq.subsVarsPartial([a],[av]);
                    bv = eq.solve(b);
                    
                    for j=1:size(etaV,2)
                        if ((etaV(j) == etah)  )
                            continue;
                        end
                        if ( (etaV(j) == etaw) )
                            continue;
                        end
                        etak = etaV(j);
                        nc = nc+1;
                        c(nc) = etah - etak;  % <= 0   easier for substitution
                        c(nc) = c(nc).subsVarsPartial([a],[av]);
                        c(nc) = c(nc).subsVarsPartial([b],[bv]);
                    end
                    obj0 = etah + functionF(a*x+b*y);
                    obj0 = obj0.subsVarsPartial([a],[av]);
                    obj0 = obj0.subsVarsPartial([b],[bv]);
                    deg = polynomialDegree(obj0.f,z);
                    objfacts = coeffs(obj0.f,z);
                    
                    if (deg == size(objfacts,2)+1)
                        psi0 = objfacts(1);
                        psi1 = objfacts(2)/2;
                        psi2 = -objfacts(3);
                    elseif (deg == size(objfacts,2))
                        psi0=0;
                        psi1 = objfacts(1)/2;
                        psi2 = -objfacts(2);
                    else
                        psi0=0;if (degreeh==2 & degreew==1)
                    disp("one quad")
                    z = sym('z');
                    av = z - mh*b;
                    eq = etah-etaw;
                    eq = eq.subsVarsPartial([a],[av]);
                    bv = eq.solve(b);
                    
                    for j=1:size(etaV,2)
                        if ((etaV(j) == etah)  )
                            continue;
                        end
                        if ( (etaV(j) == etaw) )
                            continue;
                        end
                        etak = etaV(j);
                        nc = nc+1;
                        c(nc) = etah - etak;  % <= 0   easier for substitution
                        c(nc) = c(nc).subsVarsPartial([a],[av]);
                        c(nc) = c(nc).subsVarsPartial([b],[bv]);
                    end
                    obj0 = etah + functionF(a*x+b*y);
                    obj0 = obj0.subsVarsPartial([a],[av]);
                    obj0 = obj0.subsVarsPartial([b],[bv]);
                    deg = polynomialDegree(obj0.f,z);
                    objfacts = coeffs(obj0.f,z);
                    
                    if (deg == size(objfacts,2)+1)
                        psi0 = objfacts(1);
                        psi1 = objfacts(2)/2;
                        psi2 = -objfacts(3);
                    elseif (deg == size(objfacts,2))
                        psi0=0;
                        psi1 = objfacts(1)/2;
                        psi2 = -objfacts(2);
                    else
                        psi0=0;
                        psi1=0;
                        psi2 = -objfacts(1);
                    end
                    if positivePsi (obj,-psi2,x,y)==1
                       disp('+ to be implemented')
                    else
                       f0 = simplify(psi1^2/psi2+psi0);
                       r0 = -simplify(etaR(ix(i),3).f*psi2-psi1);
                       
                       envfs = [envfs, f0];
                       envds = [envds, r0];
                       f0 = -etaR(ix(i),3).f^2*psi2 +2*etaR(ix(i),3).f*psi1+psi0;
                       r0 = simplify(etaR(ix(i),3).f*psi2-psi1);
                       envfs = [envfs, f0];
                       envds = [envds, r0];
                            
                    end
                  end
                
                        psi1=0;
                        psi2 = -objfacts(1);
                    end
                    if positivePsi (obj,-psi2,x,y)==1
                       disp('+ to be implemented')
                    else
                       f0 = simplify(psi1^2/psi2+psi0);
                       r0 = -simplify(etaR(ix(i),3).f*psi2-psi1);
                       
                       envfs = [envfs, f0];
                       envds = [envds, r0];
                       f0 = -etaR(ix(i),3).f^2*psi2 +2*etaR(ix(i),3).f*psi1+psi0;
                       r0 = simplify(etaR(ix(i),3).f*psi2-psi1);
                       envfs = [envfs, f0];
                       envds = [envds, r0];
                            
                    end
                end
                

                
            end     
        end

        function l=positivePsi (obj,p,x,y)
            l = 0;
            for i = 1: obj.d.nVertices
               if (subs(p,[x,y],[obj.d.vx(i),obj.d.vy(i)]) < 0)
                   return
               end
            end
            l = 1;
        end

        %feasible pairs
        % (ix,jx) feasible pair
        % (vix,vjx) 0 : from V, 1 : from E
        % ixd,jxd : if from E - gives region no - 1,2,3
        function [ix,jx,vix, vjx, ixd, jxd] = feasiblePairs (obj, etaR,a,b)
            n = 0;
            for i = 1:size(etaR,1)
                for j = 1:size(obj.d.V,2)
                    vertex = [obj.d.getVertex(obj.d.V(j))];
                    s = etaR(i,1).subsVarsPartial ([a,b],vertex);
                    
                    if (s <= etaR(i,2))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 1;
                      jxd(n) = 0;
                    end
                    if (s >= etaR(i,3))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 3;
                      jxd(n) = 0;
                    end
                    if ((etaR(i,2)<=s)&(s<=etaR(i,3)))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 2;
                      jxd(n) = 0;
                    end
                    
                end
            end
            for i = 1:size(etaR,1)
                for j = i+1:size(etaR,1)
                      for i1 = 1:3
                          for j1=1:3
                            n = n + 1;
                            ix(n) = i;
                            jx(n) = j;
                            vix(n) = 1;
                            vjx(n) = 1;
                            ixd(n) = i1;
                            jxd(n) = j1;
                          end
                      end
                end
            end
        end 
        
        function [etaV, etaE, etaR] = getEtaFunctions (obj,x,y,a,b)
            eta = obj.f - functionF(a*x+b*y);
            
            % Eta for Edges
            for i = 1:obj.d.nE
                xv1 = obj.d.vx(obj.d.E(i,1));
                yv1 = obj.d.vy(obj.d.E(i,1));
                etaE(i,1) = eta.subsVarsPartial([x,y],[xv1,yv1]);
            
                xv2 = obj.d.vx(obj.d.E(i,2));
                yv2 = obj.d.vy(obj.d.E(i,2));
                etaE(i,3) = eta.subsVarsPartial([x,y],[xv2,yv2]);
            
                edgey = obj.d.mE(i) * x + obj.d.cE(i);
                etaT = eta.subsVarsPartial([y],[edgey]);
                df = etaT.dfdx(x);
                xp = df.solve(x);
                etaE(i,2) = etaT.subsVarsPartial([x],[xp]); 
                
                

                f1 = obj.f.subsVarsPartial([y],[edgey]);
                df1 = f1.dfdx(x);

                
                etaR(i,1) = functionF(a+obj.d.mE(i)*b);
                df1 = functionF(a+obj.d.mE(i)*b);

                b1 = df1.subsVarsPartial ([a,b],[xv1,yv1]);
                b2 = df1.subsVarsPartial ([a,b],[xv2,yv2]);
                if (double(b1.f) < double(b2.f))
                  etaR(i,2) = b1;  %subs(df1,[x,y],[xv1,yv1]);
                  etaR(i,3) = b2; %subs(df1,[x,y],[xv2,yv2]);
                else
                  etaR(i,3) = b1;  %subs(df1,[x,y],[xv1,yv1]);
                  etaR(i,2) = b2; %subs(df1,[x,y],[xv2,yv2]);
                  t = etaE(i,1);
                  etaE(i,1) = etaE(i,3);
                  etaE(i,3) = t;
                end
            end

            % Eta for Vertices
            for i = 1:obj.d.nV
                xv = obj.d.vx(obj.d.V(i));
                yv = obj.d.vy(obj.d.V(i));
                etaV(i) = eta.subsVarsPartial([x,y],[xv,yv]);
            end
            
        end
    end
end