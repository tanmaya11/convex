classdef plq_1piece
    properties
        f;
        d;
        envf=functionF.empty();

        % fix envd to be intersection of ineqs
        % envd=functionF.empty();
        envd = region.empty();
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
            %envfT
            obj.envf = envfT;
            obj.envd = envdT;
            return

            % to be implemented
            envfT=[]
            ndTs = []
            envdTs=[]
            n = 0
            for i = 1:size(obj.envd,2)
                for j = i+1:size(obj.envd,2)
                    z0 = solve (obj.envd(i).f<=0,obj.envd(j).f<0);
                    if isempty (z0.x) | isempty (z0.y) 
                        obj.envd(i).print
                        obj.envd(j).print
                        disp('no intersection')
                    else
                        isSubset = isAlways((obj.envd(i).f<=0) <= (obj.envd(j).f<0))
                        isSubset = isAlways((obj.envd(j).f<=0) <= (obj.envd(i).f<0))
                        obj.envd(i).print
                        obj.envd(j).print
                        disp(' intersection')
                        solve (obj.envd(i).f<=0,obj.envd(j).f<0)
                        n = n + 1
                        envfT(n) = pointwise_max(obj.envf(i), obj.envf(j), obj.d.vx, obj.d.vy, obj.d.ineqs, [obj.envd(j)])
                       
                    end
                end
            end
            
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
            %return
            % put a check that eta are only polynomials 

             % (vix, vjx) is the pair :  1 for edge else 0 
            disp("feasiblePairs") 
            [ix,jx,vix, vjx, ixd, jxd] = feasiblePairs (obj,etaR, a,b);
            %ix=[2]
            %jx=[1]
            %vix=[1]
            %vjx=[0]
            %ixd=[2]
            %jxd=[0]
            %return
            disp("solve")
            [envfs, envds] = solve (obj, ix,jx,vix, vjx,ixd, jxd, etaV, etaE, etaR,a, b, x, y);
            
            for i = 1:size(envfs,2)
                obj.envf = [obj.envf, envfs(i)];
                obj.envd = [obj.envd, envds(i).removeDenominator];
                
            end
             
            % put code for max 
        end 

        function getIntersections (obj)
            for i=1:size(obj.envd,2)
                for j=i+1:size(obj.envd,2)
                    if isempty(solve(obj.envd(i).f <0, obj.envd(j).f <0))
                        continue
                    end
                    i
                    j
                    obj.envd(i).print
                    obj.envd(j).print
                end
            end
        end

        function [etah, m, yintercept] = edgeInfoInSolve(obj, etaE, etaR, ix, ixd) %, nc, c)
           etah = etaE(ix,ixd);
           if size(obj.d.mE,2) == 1
             m = obj.d.mE;
             yintercept = obj.d.cE;
           else  
             m = obj.d.mE(ix);
             yintercept = obj.d.cE(ix);
           end
           
        end

        function [nb,lb, ub] = getBound1 (obj,linfeasible, c,b,nb,lb,ub)
            
          c1 =  c.solve(b);
          linfeasible = isempty(c1);
          if (linfeasible) 
                      return;
          end
          s = coeffs(c.f,b);
          if (size(c1,1) > 1)
              if(s(end) > 0)
              nb = nb+1;
              lb(nb)=c1(1);
              ub(nb)=c1(2);
              end
          else
            nb = nb+1;
            
              if s(end) < 0
                ub(nb)= inf;
                lb(nb) = c1(1);
              else
                lb(nb)= -inf;
                ub(nb) = c1(1);
              end
            
          end
          
        end

            function [nb, lb, ub] = getBoundsLinear (obj, linfeasible, etaR,ixd, a,b,av, nb, lb, ub)
            if ixd == 1
              nc = 1  ;
              c(1) = etaR(1)-etaR(2);
            elseif ixd == 2
                nc = 2;
              c(1) = etaR(2)-etaR(1);
              c(2) = etaR(1)-etaR(3);
              
            elseif ixd == 3
                nc = 1;
              c(1) = etaR(3)-etaR(1);
            end
            for i = 1:nc
              c(i) = c(i).subsVarsPartial([a],[av]);
              [nb,lb, ub] = getBound1 (obj,linfeasible, c(i),b,nb,lb,ub);
              if (linfeasible) 
                      return;
                  end
            end
        end
    
        function [envfs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRix, ixd, etaRjx, jxd, x, y, a, b, envfs, envds) 
            f0 = etah-etaw;
            linfeasible = false;
            av = f0.solve(a);
            if (isempty(av))
                lSol = false;
                return;
            end
            lSol = true;
            obj0 = obj0.subsVarsPartial([a],[av]);
            objfacts = coeffs(obj0.f,b);
            
            nb = 0 ;
            lb = [];
            ub = [];
            if (ixd > 0)
               [nb, lb, ub] = getBoundsLinear (obj, linfeasible, etaRix,ixd, a,b,av, nb, lb, ub);
               if (linfeasible) 
                      return;
                  end
            end
            if (jxd > 0)
              [nb, lb, ub] = getBoundsLinear (obj, linfeasible, etaRjx,jxd, a,b,av, nb, lb, ub);
              if (linfeasible)  
                      return;
                  end
            end
            for j=1:size(etaV,2)
               if (lV(j)) 
                 continue;
               end
               etak = etaV(j);
               c = etah - etak ; % <= 0   easier for substitution
               c = c.subsVarsPartial([a],[av]);
               [nb,lb, ub] = getBound1 (obj,linfeasible, c,b,nb,lb,ub);
               if (linfeasible) 
                      return;
                  end
            end
            
            for j=1:size(etaE,1)
               if (lE(j))
                 continue;
               end
               for k = 1:size(etaE,2) 
                  etak = etaE(j,k);
                  c = etah - etak;  % <= 0   easier for substitution
                  c = c.subsVarsPartial([a],[av]);
                  if (isZero(c)) 
                      continue
                  end
                  [nb,lb, ub] = getBound1 (obj,linfeasible, c,b,nb,lb,ub);
                  if (linfeasible) 
                      return;
                  end
                end
            end
            
            mlb = max(lb);
            mub = min(ub);
            
            %for i = 1:size(lb,2)
            if (mlb > mub)
                disp('infeasible')
            else
              if (mub == inf)
                envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                envds = [envds, region(objfacts(2)) ];
              elseif (mlb == -inf)
                envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                envds = [envds, region(-objfacts(2)) ];
              elseif (mlb == mub)
                envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                envds = [envds, region(objfacts(2)) ];
                
              else
                envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                envds = [envds, region(objfacts(2)) ];
                envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                envds = [envds, region(-objfacts(2)) ];
              end
              
            %end
            end
                    
            return
            
        end

        
        
        function [envfs, envds] = solveQuadLinear1 (obj, m, a, b, x, y, etah, etaw, etaR, etaV, lV, etaE, lE, ix, envfs, envds)
          disp("one quad")
          z = sym('z');
          %av = z - m*b;
          %z = a+m*b;
          %lSolve = true;
          av = solve(z-a-m*b,a);
          eq = etah-etaw
          eq = eq.subsVarsPartial([a],[av])
          bv = eq.solve(b)

          if (isempty(bv))
              eq
              zv = solve(eq,z)
              
              obj0 = etah + functionF(a*x+b*y)
              obj0 = obj0.subsVarsPartial([a],[av])
              obj0 = obj0.subsVarsPartial([z],[zv])
          
              disp ("Fix this")
           %   lSolve = false;
              return
          end
          nb = 1;
          lb(nb) = etaR(ix,2).f;
          ub(nb) = etaR(ix,3).f;

            

          % solve these to get bounds
          for j=1:size(etaV,2)
            if (lV(j))
              continue;
            end
            etak = etaV(j);
            c = etah - etak;  % <= 0   easier for substitution
            c = c.subsVarsPartial([a],[av]);
            c = c.subsVarsPartial([b],[bv]);
            if isNegativeSqr(c,z)
              disp ('negative sqr')
              continue
            end
            z0 = c.solve(z)
            % check this part
            % as z is always satisfied
            if isreal(z0(1))
              nb = nb+1;
              lb(nb) = min(z0);
              ub(nb) = max(z0);
            end
          end
          lb
          ub
          for j=1:size(etaE,1)
            if (lE(j))
              continue;
            end
            for k = 1:size(etaE,2)
              etak = etaE(j,k);
              %nc = nc+1;
              c = etah - etak;  % <= 0   easier for substitution""
              c = c.subsVarsPartial([a],[av]);
              c = c.subsVarsPartial([b],[bv]);
              %c.print
              if isNegativeSqr(c,z)
                  disp ('negative sqr')
                  continue
                  
              end 
              cz = coeffs(c.f,z);
              z1 = roots(cz);
              %z0 = c.solve(z)
              %disp('size')
              %size(z0)
              nz = 0;
              for iz = 1:size(z1)
                  if (~isreal(z1(iz)))
                      continue
                  end
                  nz = nz + 1;
                  z0(nz) = z1(iz);
              end
              if nz > 0
                nb = nb+1;
                lb(nb) = min(z0);
                ub(nb) = max(z0);
              end
              
            end
          end
                   
          %lb
          %ub
          obj0 = etah + functionF(a*x+b*y);
          obj0 = obj0.subsVarsPartial([a],[av]);
          obj0 = obj0.subsVarsPartial([b],[bv])
          deg = polynomialDegree(obj0.f,z);
          if (deg ~= 2) 
              disp("check degree of polynomial in z in quad-lin")
              return
          end
          objfacts = coeffs(obj0.f,z)
          
          if (deg+1 == size(objfacts,2))
            psi0 = objfacts(1);
            psi1 = objfacts(2)/2;
            psi2 = -objfacts(3);
          elseif (deg == size(objfacts,2))
            z0 = simplify(obj0.f - objfacts(end)*z^2);
            %polynomialDegree(z0,z)
            if polynomialDegree(z0,z) == 1
              psi0=0;
              psi1 = objfacts(1)/2;
            else
              psi0=objfacts(1);
              psi1 =0;

            end

            psi2 = -objfacts(2);
          else
            psi0=0;
            psi1=0;
            psi2 = -objfacts(1);
          end

          psi0
          psi1
          psi2
          lb
          ub

          mlb = max(lb);
          mub = min(ub);
            
            %for i = 1:size(lb,2)
            
          if positivePsi (obj,-psi2,x,y)==1
            disp("-psi2 +ve")
            psi1
            psi2
            f0 = simplify(psi1/psi2);
            f1 = etaR(ix,3).f;
            
            r0 = simplify(psi1 - etaR(ix,3).f*psi2);
            envfs = [envfs, functionF(f0)];
            envds = [envds, region(r0)];
            r0 = -simplify(psi1 - etaR(ix,3).f*psi2);
            envfs = [envfs, functionF(f1)];
            envds = [envds, region(r0)];
          else
          
            
            for i = 1:size(lb,2)
                mlb = lb(i)
                mub = ub(i)
            if (mlb >= mub)
                disp('infeasible')
            else
          

            %for i = 1:nb
              f0 = simplify(psi1^2/psi2+psi0)
              r0 = -simplify(mub*psi2-psi1);
              r1 = simplify(mlb*psi2-psi1);
              % put in r1
              disp("in quad-lin")
              f0
              r0
              size(envds)
              envfs = [envfs, functionF(f0)];
              envds = [envds, region(r0)];
              size(envds)
              envds(1).print
              envfs = [envfs, functionF(f0)];
              envds = [envds, region(r1)];
              
              
              f0 = -mub^2*psi2 +2*mub*psi1+psi0 ;
              r0 = simplify(mub*psi2-psi1);
              envfs = [envfs, functionF(f0)];
              envds = [envds, region(r0)];


              f0 = -mlb^2*psi2 +2*mub*psi1+psi0 ;
              r0 = simplify(-mlb*psi2+psi1);
              envfs = [envfs, functionF(f0)];
              envds = [envds, region(r0)];
              size(envds)
            end
            end
            
            %
            %envfs = [envfs, f0];
            %envds = [envds, r0];
            %r0 = simplify(etaR(ix,3).f*psi2-psi1);
            %envfs = [envfs, f1];
            %envds = [envds, r0];
          end
        end
           
        function [envfs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0, alpha1, mh, qh, ix, jx, etaR, etaV, lV, etaE, lE, envfs, envds)
            disp("quadquad")
          av = alpha1*b + alpha0;
          
          obj0 = etah + functionF(a*x+b*y);
                    
          obj0 = obj0.subsVarsPartial([a],[av]);
          
          nb = 1;
          
          c = etaR(ix,1) - etaR(ix,2);
          c = c.subsVarsPartial([a],[av])
          coef = coeffs(c.f)
          if (coef(end) > 0)
            lb(nb) = c.solve(b);
            c = etaR(ix,1) - etaR(ix,3);
            c = c.subsVarsPartial([a],[av]);
            ub(nb) = c.solve(b);
          else
            ub(nb) = c.solve(b);
            c = etaR(ix,1) - etaR(ix,3);
            c = c.subsVarsPartial([a],[av]);
            lb(nb) = c.solve(b);
          end
          
          %lb(nb) = etaR(ix,2).f
          %ub(nb) = etaR(ix,3).f
          nb = 2;
          c = etaR(jx,1) - etaR(jx,2);
          c = c.subsVarsPartial([a],[av]);
          coef = coeffs(c.f);
          if (coef(end) > 0)
            lb(nb) = c.solve(b);
            c = etaR(jx,1) - etaR(jx,3);
            c = c.subsVarsPartial([a],[av]);
            ub(nb) = c.solve(b);
          else
            ub(nb) = c.solve(b);
            c = etaR(jx,1) - etaR(jx,3);
            c = c.subsVarsPartial([a],[av]);
            lb(nb) = c.solve(b);
          end
         
          lb
          ub
          
          %lb(nb) = etaR(jx,2).f
          %ub(nb) = etaR(jx,3).f

          for j=1:size(etaV,2)
            if (lV(j))
              continue;
            end
            etak = etaV(j)
           % nc = nc+1;
            % make this a function
            c = etah - etak  % <= 0   easier for substitution
            c = c.subsVarsPartial([a],[av]);
            c.print
            disp('here')
            z0 = c.solve(b);
            nb = nb+1;
            if z0(1) < z0(2)
              lb(nb) = z0(1);
              ub(nb) = z0(2);
            else
              ub(nb) = z0(1);
              lb(nb) = z0(2);
            end
          end
          lb
          ub
          %return
          for j=1:size(etaE,1)
            if (lE(j))
              continue;
            end
            for k = 1:size(etaE,2)
              etak = etaE(j,k)
              %nc = nc+1;
              c = etah - etak;  % <= 0   easier for substitution
              c = c(nc).subsVarsPartial([a],[av]);
              c.print
              z0 = c.solve(b);
              nb = nb+1;
              if z0(1) < z0(2)
              lb(nb) = z0(1);
              ub(nb) = z0(2);
            else
              ub(nb) = z0(1);
              lb(nb) = z0(2);
              end
            end
          end
          lb;
          ub;
          mlb = max(lb);
          mub = min(ub);
          if mlb < mub
          psi0=-(alpha0-qh)^2/(4*mh)+alpha0*x;
          psi1=-(alpha1+mh)*(alpha0-qh)/(2*mh) + y -qh + alpha1*x;
          psi2=(alpha1+mh)^2/(4*mh);
          %for i = 1: nb
              f0 = simplify(psi1^2/psi2+psi0);
              %fix this
              %r0 = simplify(psi1<ub(i)*psi2-psi1);
              envfs = [envfs, f0];
              % check this
              envds = [envds, region(f0)];
              f0 = simplify(-mlb^2*psi2 + 2*mlb*psi1 + psi0);
              r0 = simplify(psi1-mlb*psi2);
              envfs = [envfs, f0];
              envds = [envds, region(r0)];
              f0 = simplify(-mub^2*psi2 + 2*mub*psi1 + psi0);
              r0 = simplify(mub*psi2-psi1);
              envfs = [envfs, f0];
              envds = [envds, region(r0)];
              
          %end
          end
        end


        function [envfs, envds] = solve (obj, ix,jx,vix, vjx, ixd, jxd, etaV, etaE, etaR, a, b, x, y)
            envfs = [];

            envds = [];
            disp("solve")
            
            for i=1:size(ix,2)
              i
              for j = 1:size(etaV,2)
                lV(j) = false;
              end
              for j = 1:size(etaE,1)
                lE(j) = false;
              end
            
                
              %nc = 0;
              if (vix(i)==1)
                  %c(1)=etaR(1,1);
                  [etah, mh, qh] = edgeInfoInSolve(obj, etaE, etaR, ix(i), ixd(i)); %, nc, c);
                  lE(ix(i)) = true;
              else
                  etah = etaV(ix(i));
                  lV(ix(i)) = true;
              end
              if (vjx(i)==1)
                  [etaw, mw, qw] = edgeInfoInSolve(obj, etaE, etaR, jx(i), jxd(i)); %, nc, c);
                  lE(jx(i)) = true;
              else
                  etaw = etaV(jx(i));
                  lV(jx(i)) = true;
              end
              etah
              etaw
              degreeh = polynomialDegree(etah.f);
              degreew = polynomialDegree(etaw.f);
              %continue
                if (degreeh==1 & degreew==1)
                    disp("lin-lin")
                    %continue
                    obj0 = etah + functionF(a*x+b*y);
                    if ixd(i) == 0
                        etaRi=functionF;
                    else
                        etaRi = etaR(ix(i),:);
                    end
                    if jxd(i) == 0
                        etaRj=functionF;
                    else
                        etaRj = etaR(jx(i),:);
                    end
                    %etaRi.printL
                    %etaRj.printL
                    %ixd(i)
                    %jxd(i)
                    %disp('Out')
                    %envfs
                    [envfs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRi, ixd(i), etaRj, jxd(i), x, y, a, b, envfs, envds) ;
                    % fix if ix or jx doesnt exist
                    
                    if (~lSol)
                        disp("lSol")
                        [envfs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRi, ixd(i), etaRj, jxd(i), x, y, b, a, envfs, envds) ;
                    end

                    %disp("envfs0")
                    %for i = 1: size(envfs,2)
                    %    envfs(i).print
                    %end

                    
                end 
                

                if (degreeh==2 & degreew==1)
                    disp("quad-lin")
                  %  continue;
                    [envfs, envds] = solveQuadLinear1 (obj, mh, a, b, x, y, etah, etaw, etaR, etaV, lV, etaE, lE, ix(i), envfs, envds);
            %        lSolve
            %        if ~lSolve
            %            disp("b a flipped")
            %            [envfs, envds, lSolve] = solveQuadLinear1 (obj, mh, a, b, b, a, x, y, etah, etaw, etaR, etaV, lV, etaE, lE, ix(i), envfs, envds);
            %        end
            % flipping a,b gives same answer
                end
                
                if (degreeh==1 & degreew==2)
                    disp("lin-quad")
                    
                %   continue
                    [envfs, envds] = solveQuadLinear1 (obj, mw, a, b, x, y, etaw, etah, etaR, etaV, lV, etaE, lE, jx(i), envfs, envds)
                end
                            
                if (degreeh==2 & degreew==2)
                    disp("quad-quad")
                    
                   % continue
                    obj0 = etah + functionF(a*x+b*y);
                    disp("h1")
                    if (mh == mw)
                       %av = (2*mh*b+qh+qw)/2;
                       alpha1 = mh;
                       alpha0 = (qh+qw)/2
                       [envfs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0, alpha1, mh, qh, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envds)
                    else
                       disp("mh /= mw")
                       av = (qw*mh-qh*mw+sqrt(mh*mw)*((mh-mw)*b+qh-qw))/(mh-mw)
                       
                       alpha = coeffs(av,b)
                       size(alpha)
                       if (size(alpha) == 2)
                           alpha1 = alpha(2)
                           alpha0 = alpha(1)
                       else
                           alpha1 = alpha(1)
                           alpha0 = 0
                       end
                       [envfs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0,  alpha1, mh, qh, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envds)
                       disp("second")
                       av = (qw*mh-qh*mw-sqrt(mh*mw)*(mh-mw)*b+qh-qw)/(mh-mw)
                       alpha = coeffs(av,b)
                       size(alpha)
                       if (size(alpha) == 2)
                           alpha1 = alpha(2)
                           alpha0 = alpha(1)
                       else
                           alpha1 = alpha(1)
                           alpha0 = 0
                       end

                       
                       [envfs, envds] = solveQuadQuad1(obj, etah,  x, y, a, b, alpha0, alpha1, mh, qh, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envds)
                    end
                end
                if (i == 5) 
                  %  return
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
             %       i
              %      j
               %     s
                %    vertex
                 %   etaR(i,2)
                 %   etaR(i,3)
                    
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
                
                
                %obj.f
                f1 = obj.f.subsVarsPartial([y],[edgey]);
                %f1.f
                df1 = f1.dfdx(x);
                %df1
                
                etaR(i,1) = functionF(a+obj.d.mE(i)*b);
                %df1 = functionF(a+obj.d.mE(i)*b)

                b1 = df1.subsVarsPartial ([x,y],[xv1,yv1])
                b2 = df1.subsVarsPartial ([x,y],[xv2,yv2])
                if (double(b1.f) < double(b2.f))
                  etaR(i,2) = b1;  %subs(df1,[x,y],[xv1,yv1]);
                  etaR(i,3) = b2; %subs(df1,[x,y],[xv2,yv2]);
                else
                  etaR(i,3) = b1;  %subs(df1,[x,y],[xv1,yv1]);
                  etaR(i,2) = b2; %subs(df1,[x,y],[xv2,yv2]);
                  %t = etaE(i,1);
                  %etaE(i,1) = etaE(i,3);
                  %etaE(i,3) = t;
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