classdef plq_1piece
    properties
        f;
        d;
        envf=functionF.empty();
        envExpr = convexExpr.empty();
        envd = region.empty();
        conjf=functionF.empty();
        conjd = region.empty();
    end

    methods % creation & print
         function obj = plq_1piece(d,f)
            % put checks for type of f and d
            if nargin > 0
              obj.f = f;
              obj.d = d;
            end 
         end

         function print(obj)
         
           disp("Domain")
           obj.d.print
           fprintf("\n")
           disp("Function")
           obj.f.print
           fprintf("\n\n\n")
              
           disp("Convex Envelope")
           size(obj.envf)
           for j=1:size(obj.envf,2) 
             disp('Function')  
             obj.envf(j).print
             disp("Expr")
             obj.envExpr(j).print
             disp('Domain')
             obj.envd(j).print
           end
         end

    end

    methods % utility
         % li returns region no if entire region else 0
        function li = entireRegion (obj)
            li = 0;
            for i=1:size(obj.envd,2)
              if (obj.envd(i) == obj.d.polygon)
                  li = i;
                  return
              end
            end
        end

        function obj = unique(obj)
            rem = [];
             for i = 1:size(obj.envd,2)
                for j = i+1:size(obj.envd,2)
                    if (obj.envd(i) == obj.envd(j))
                        if (obj.envf(i) == obj.envf(j))
                            rem = [rem,j]
                        end
                    end
                end
             end
             obj.envd(rem) = [];
             obj.envf(rem) = [];
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


    end

    methods % convex
        function obj = convexEnvelope(obj)
            disp("in convexEnvelope")
            %for i = 1:size(obj.f,2)
              vars = obj.f.getVars;
              if (size(vars,2)==2)
                  x = vars(1);
                  y = vars(2);
              else
                  disp("not bivariate in 'plq.m'")
                  return
              end
              obj = convexEnvelope1 (obj,x,y);
              
              obj = obj.maxEnvelope([x,y]);
              %return
              li = obj.entireRegion ();
              if li > 0
                obj = obj.removeNMax (li,[x,y]);
              end
              obj = obj.maxEnvelopeIntersect([x,y]);
              obj = obj.unique();
              disp("b4 vertices")
              size(obj.envd,2)
              for j=1:size(obj.envd,2)
                obj.envd(j) = obj.envd(j).getVertices(vars);
              end
            %end
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
            %return
            %ix=[2]
            %jx=[1]
            %vix=[1]
            %vjx=[0]
            %ixd=[2]
            %jxd=[0]
            %return
            %disp("solve")
            [envfs, envxs, envds] = solve (obj, ix,jx,vix, vjx,ixd, jxd, etaV, etaE, etaR,a, b, x, y);
            %disp("solve done")
            
            for i = 1:size(envfs,2)

                %r = obj.d.polygon + envds(i).removeDenominator;
                %r = unique(r);
                r  = obj.d.polygon + envds(i);
                r = r.removeDenominator;
                r = unique(r);
                if (r.isFeasible)
                  obj.envf = [obj.envf, envfs(i)];
                  obj.envExpr = [obj.envExpr, envxs(i)];
                  obj.envd = [obj.envd, r.removeDenominator];
                  
                end
            end
            %obj.envd(1).print
            %obj.envd(2).print
            %d = intersection(obj.envd(1),obj.envd(2))
            %d.print
            % put code for max 
        end 

         % returns etah, slope and y intercept of edge
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

        %change name
        % gives bound of linear inequality c <= 0
        function [nb,lb, ub] = getBound1 (obj,linfeasible, c,b,nb,lb,ub)
            
          c1 =  c.solve(b);
          linfeasible = isempty(c1);
          if (linfeasible) 
                      return;
          end

          s = coeffs(c.f,b);
          if (size(c1,1) > 1)
              if isreal(c1(1))
              if(s(end) > 0)
              nb = nb+1;
              lb(nb)=c1(1);
              ub(nb)=c1(2);
              end
              end
          else
              if (isreal(c1(1)))
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
          
        end

        % get bounds depending on region
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
    
        

        function [envfs, envxs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRix, ixd, etaRjx, jxd, x, y, a, b, envfs, envxs, envds) 
            %if (etah == etaw) wont work
            f0 = etah-etaw;
            % if a gets eliminated we exit this routine and try again
            % exchanging a and b
            linfeasible = false;
            av = f0.solve(a);
            % Get a in terms of b 
            if (isempty(av))
                lSol = false;
                return;
            end
            lSol = true;
            
            %bcoeffs = coeffs(av,b); 
            %alpha0 = bcoeffs(1);
            %alpha1 = bcoeffs(2);

            %substitute a in objective and get coefficients of b
            obj0 = obj0.subsVarsPartial([a],[av]);
            objfacts = coeffs(obj0.f,b);
            eta0 = objfacts(1);
            eta1 = objfacts(2);
            %return
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
            
            nb0 = nb;
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
            
            
            nb1 = nb0; 
            for i = nb0+1:nb
                wB = false;
                for j = 1:nb0
                    if (lb(i) < lb(j))
                        wB = true;
                        break;
                    end
                    if (ub(i) > ub(j))
                        wB = true;
                        break;
                    end
                end
                if wB
                    continue;
                end
                nb1 = nb1+1;
                lb(nb1) = lb(i);
                ub(nb1) = ub(i);
            end
            nb = nb1;
            %for i = 1:nb
            %  mlb = lb(i);
            %  mub = ub(i);
            mlb = max(lb);
            mub = min(ub);
            %size(envfs)
            if (mlb > mub)
                disp('infeasible')
            else
              if (mub == inf)
               %   disp('case1')
                envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                envds = [envds, region(objfacts(2)) ];
              elseif (mlb == -inf)
              %    disp('case2')
             %     obj0
             %     obj0.subsVarsPartial([b],[mub])
                envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                envxs = [envxs, convexExpr(3,eta0,eta1,mub)];

                envds = [envds, region(-objfacts(2)) ];
              %elseif (mlb == mub)
              %  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
              %  envds = [envds, region(objfacts(2)) ];
                
              elseif (mlb ~= mub)
              %    disp('case3')

              % take max here and only 1 entry
                envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];

                envds = [envds, region(objfacts(2)) ];
                envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                envxs = [envxs, convexExpr(3,eta0,eta1,mub)];

                envds = [envds, region(-objfacts(2)) ];
              end
              
            end
            %end
            %size(envfs)        
            return
            
        end

        
        
        function [envfs, envxs, envds] = solveQuadLinear1 (obj, m, a, b, x, y, etah, etaw, etaRix, ixd, etaRjx, jxd, etaV, lV, etaE, lE, ix, envfs, envxs, envds)
          %disp("one quad")
          %etah 
          %etaw
          z = sym('z');
          %av = z - m*b;
          %z = a+m*b;
          %lSolve = true;
          av = solve(z-a-m*b,a);
          eq = etah-etaw;
          eq = eq.subsVarsPartial([a],[av]);
          bv = eq.solve(b);
          linfeasible = false;

          if (isempty(bv))
              return
              %zv = solve(eq,z);
              
              %obj0 = etah + functionF(a*x+b*y);
              %obj0 = obj0.subsVarsPartial([a],[av]);
              %obj0 = obj0.subsVarsPartial([z],[zv]);
          
              %disp ("Fix this")
           %   lSolve = false;
            %  return
          end
          nb = 0 ;
          lb = [];
          ub = [];
            if ixd == 1  
              nb = nb+1;
              lb(nb) = -inf;
              ub(nb) = etaRix(2).f;
            elseif ixd == 2
              nb = nb+1;
              lb(nb) = etaRix(2).f;
              ub(nb) = etaRix(3).f;
            elseif ixd == 3
              nb = nb+1;
              lb(nb) = etaRix(3).f;
              ub(nb) = inf;

            end
          if jxd == 1  
              nb = nb+1;
              lb(nb) = -inf;
              ub(nb) = etaRjx(2).f;
            elseif jxd == 2
              nb = nb+1;
              lb(nb) = etaRjx(2).f;
              ub(nb) = etaRjx(3).f;
            elseif jxd == 3
              nb = nb+1;
              lb(nb) = etaRjx(3).f;
              ub(nb) = inf;

            end  
          nb0 = nb;
          %lb
          %ub

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
              %disp ('negative sqr')
              continue
            end
            z0 = c.solve(z);
            % check this part
            % as z is always satisfied
            if isreal(z0(1))
              nb = nb+1;
              lb(nb) = min(z0);
              if (lb(nb)  < lb(1))
                  return
              end
              ub(nb) = max(z0);
              if (ub(nb)  > ub(1))
                  return
              end
              
            end
          end
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
              c.print

              % -x^2 always <= 0 hence skipping
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
          nb1 = nb0; 
            for i = nb0+1:nb
                wB = false;
                for j = 1:nb0
                    if (lb(i) < lb(j))
                        wB = true;
                        break;
                    end
                    if (ub(i) > ub(j))
                        wB = true;
                        break;
                    end
                end
                if wB
                    continue;
                end
                nb1 = nb1+1;
                lb(nb1) = lb(i);
                ub(nb1) = ub(i);
            end
            nb = nb1;
          %lb
          %ub
           
          %lb
          %ub
          obj0 = etah + functionF(a*x+b*y);
          obj0 = obj0.subsVarsPartial([a],[av]);
          obj0 = obj0.subsVarsPartial([b],[bv]);
          deg = polynomialDegree(obj0.f,z);
          if (deg ~= 2) 
              disp("check degree of polynomial in z in quad-lin")
              return
          end
          objfacts = coeffs(obj0.f,z);
          
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

          %psi0
          %psi1
          %psi2
          %lb
          %ub

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
          
          mlb = min(lb);
          mub = max(ub);
          if mlb == -inf
              return
          elseif mub == inf
              return
          end
            %for i = 1:size(lb,2)
             %   mlb = lb(i);
             %   mub = ub(i);
            if (mlb >= mub)
                disp('infeasible')
            else
          

            %for i = 1:nb
              %f0 = simplify(psi1^2/psi2+psi0);
              f0 = functionF(psi1^2,psi2) + functionF(psi0);
              r0 = -simplify(mub*psi2-psi1);
              r1 = simplify(mlb*psi2-psi1);
              % put in r1
              %envfs = [envfs, functionF(f0)];
              envfs = [envfs, f0];
              envxs = [envxs, convexExpr(1,psi0,psi1,psi2)];
              envds = [envds, region([r0,r1])];
              %envfs = [envfs, functionF(f0)];
              %envds = [envds, region(r1)];
              
              
              f0 = -mub^2*psi2 +2*mub*psi1+psi0 ;
              r0 = simplify(mub*psi2-psi1);
              envfs = [envfs, functionF(f0)];
              envxs = [envxs, convexExpr(2,psi0,psi1,psi2)];
              
              envds = [envds, region(r0)];


              f0 = -mlb^2*psi2 +2*mlb*psi1+psi0;
              r0 = simplify(-mlb*psi2+psi1);
              envfs = [envfs, functionF(f0)];
              envxs = [envxs, convexExpr(2,psi0,psi1,psi2)];
              
              envds = [envds, region(r0)];
              %size(envds)
            %end
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
          c = c.subsVarsPartial([a],[av]);
          coef = coeffs(c.f);
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
         
          %lb
          %ub
          
          %lb(nb) = etaR(jx,2).f
          %ub(nb) = etaR(jx,3).f

          for j=1:size(etaV,2)
            if (lV(j))
              continue;
            end
            etak = etaV(j)
           % nc = nc+1;
            % make this a function
            c = etah - etak;  % <= 0   easier for substitution
            c = c.subsVarsPartial([a],[av]);
            %c.print
            %disp('here')
            z0 = c.solve(b);
            if isreal(z0(1))
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
          %lb
          %ub
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
           %   c.print
              z0 = c.solve(b);
              if isreal(z0(1))
            
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
          end
          lb;
          ub;
          mlb = max(lb);
          mub = min(ub);
          %if mlb < mub
          psi0=-(alpha0-qh)^2/(4*mh)+alpha0*x;
          psi1=-(alpha1+mh)*(alpha0-qh)/(2*mh) + y -qh + alpha1*x;
          psi2=(alpha1+mh)^2/(4*mh);
          %for i = 1: nb
          %    mlb = lb(i);
           %   mub = ub(i);
              if mlb > mub
                  disp('infeasible')
              else
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
              end
          %end
          %end
        end


        % solves all subproblems
        function [envfs, envxs, envds] = solve (obj, ix,jx,vix, vjx, ixd, jxd, etaV, etaE, etaR, a, b, x, y)
          envfs = [];
          envxs = [];
          envds = [];
          for i=1:size(ix,2)
             i
             %if i > 8
                % return
             %end
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
              if (etah == etaw)
                  continue;
              end
              degreeh = polynomialDegree(etah.f);
              degreew = polynomialDegree(etaw.f);
              if (degreeh==1 & degreew==1)
                    disp("lin-lin")
                    %continue
                    %objective function set here as we can exchange a and b
                    %if required

                    %etah
                    %etaw
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
                    %size(envfs,2)
                    [envfs, envxs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRi, ixd(i), etaRj, jxd(i), x, y, a, b, envfs, envxs, envds) ;
                    %size(envfs,2)
                    % fix if ix or jx doesnt exist
                    
                    if (~lSol)
                        disp("lSol")
                        [envfs, envxs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRi, ixd(i), etaRj, jxd(i), x, y, b, a, envfs, envxs, envds) ;
                    end

                    %for i = 1:size(envfs,2)
                    %    disp('envfs')
                    %    envfs(i).print
                    %    disp('envds')
                        
                    %    envds(i).print
                    %end
                    
                end 
                
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
                if (degreeh==2 & degreew==1)
                    disp("quad-lin")
                    %continue;
                    
                    [envfs, envxs, envds] = solveQuadLinear1 (obj, mh, a, b, x, y, etah, etaw, etaRi, ixd(i), etaRj, jxd(i), etaV, lV, etaE, lE, ix(i), envfs, envxs,envds);
            % flipping a,b gives same answer
                end
                
                if (degreeh==1 & degreew==2)
                    disp("lin-quad")
                    %continue;
                    etah
                    etaw


                    % check if lV lE need to be updated -don't think so tho
                    
                    [envfs, envxs, envds] = solveQuadLinear1 (obj, mw, a, b, x, y, etaw, etah, etaRj, jxd(i), etaRi, ixd(i), etaV, lV, etaE, lE, jx(i), envfs, envxs,envds);
                end
                            
                if (degreeh==2 & degreew==2)
                    disp("quad-quad")
                    continue;
                    obj0 = etah + functionF(a*x+b*y);
                    disp("h1")
                    if (mh == mw)
                       %av = (2*mh*b+qh+qw)/2;
                       alpha1 = mh;
                       alpha0 = (qh+qw)/2
                       [envfs, envxs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0, alpha1, mh, qh, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envxs, envds);
                    else
                       disp("mh /= mw")
                       c0 = (qw*mh-qh*mw+sqrt(mh*mw))/(mh-mw);
                       %av = (qw*mh-qh*mw+sqrt(mh*mw)*((mh-mw)*b+qh-qw))/(mh-mw)
                       av = c0 * ((mh-mw)*b+qh-qw);
                       alpha = coeffs(av,b);
                       if (size(alpha) == 2)
                           alpha1 = alpha(2);
                           alpha0 = alpha(1);
                       else
                           alpha1 = alpha(1);
                           alpha0 = 0;
                       end
                       [envfs, envxs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0,  alpha1, mh, qh, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envxs, envds);
                       disp("second")
                       av = (qw*mh-qh*mw-sqrt(mh*mw)*(mh-mw)*b+qh-qw)/(mh-mw);
                       alpha = coeffs(av,b);
                       if (size(alpha) == 2)
                           alpha1 = alpha(2);
                           alpha0 = alpha(1);
                       else
                           alpha1 = alpha(1);
                           alpha0 = 0;
                       end

                       
                       [envfs, envxs, envds] = solveQuadQuad1(obj, etah,  x, y, a, b, alpha0, alpha1, mh, qh, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envxs, envds);
                    end
                end
                    
            end 
           
        end

        

        function l=positivePsi (obj,p,x,y)
            l = 0;
            for i = 1: obj.d.polygon.nv
               if (subs(p,[x,y],[obj.d.polygon.vx(i),obj.d.polygon.vy(i)]) < 0)
                   return
               end
            end
            l = 1;
        end


                %currently edited to take all pairs - to be fixed
        %feasible pairs
        % (ix,jx) feasible pair
        % (vix,vjx) 0 : from V, 1 : from E
        % ixd,jxd : if from E - gives region no - 1,2,3
        function [ix,jx,vix, vjx, ixd, jxd] = feasiblePairs (obj, etaR,a,b)
            n = 0;
            x1 = sym('x1');
            y1 = sym('y1');
            for i = 1:size(etaR,1)
                vj1 = obj.d.E(i,1);  
                vertex = [obj.d.getVertex(vj1)];
                vjx1 = vertex(1);
                vj2 = obj.d.E(i,2);  
                vertex = [obj.d.getVertex(vj2)];
                vjx2 = vertex(1);
                if (vertex(1) < vjx1)
                    vjx2 = vjx1;
                    vjx1 = vertex(1);
                end
                eq1 = y1 - obj.d.mE(i)*x1 - obj.d.cE(i);
                for j = 1:size(obj.d.V,2)
                    vertex = [obj.d.getVertex(obj.d.V(j))];

                    c = vertex(2) + vertex(1) * (1/obj.d.mE(i));
                    eq2 = y1 + (1/obj.d.mE(i))*x1 - c;
                    proj = solve(eq1==0,eq2==0);
                    s = proj.x1;
                    subs(eq1,proj);
                    subs(eq2,proj);
                    %s = etaR(i,1).subsVarsPartial ([a,b],vertex);
                    % s = df(x) = 2mx + q
                    %s = functionF(2*obj.d.mE(i)*vertex(1) + obj.d.cE(i));
         
                    %s = vertex(1)
                    %obj.d.mE(i)
                    %obj.d.cE(i)
                    %i
                    %j
                    %s
                %    vertex
                    %etaR(i,2)
                    %etaR(i,3)
                    
                    %if ((etaR(i,2)<=s)&(s<=etaR(i,3)))
                    %if ((vjx1<=s)&(s<=vjx2))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 2;
                      jxd(n) = 0;
                    %end  
                    %if (s <= vjx2)
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 1;
                      jxd(n) = 0;
                    %end  
                    %if (s >= vjx1)
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 3;
                      jxd(n) = 0;
                    %end
                    
                end
            end
            for i = 1:size(etaR,1)
                for j = i+1:size(etaR,1)
                      for i1 = 1:3
                          for j1=1:3
                              %if (i1==1 & j1==1)
                              %    continue
                              %end
                              %if (i1==3 & j1==3)
                              %    continue
                              %end
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
            return
            % vertex vertex cases
            size(obj.d.V,2)
            for i = 1:size(obj.d.V,2)
              for j = i+1:size(obj.d.V,2)
                 n = n + 1;
                 ix(n) = i;
                 jx(n) = j;
                 vix(n) = 0;
                 vjx(n) = 0;
                 ixd(n) = 0;
                 jxd(n) = 0;
                    
              end
            end
            

        end 
        
        function [etaV, etaE, etaR] = getEtaFunctions (obj,x,y,a,b)
            eta = obj.f - functionF(a*x+b*y);
            
            % Eta for Edges
            for i = 1:obj.d.nE
                xv1 = obj.d.polygon.vx(obj.d.E(i,1));
                yv1 = obj.d.polygon.vy(obj.d.E(i,1));
                etaE(i,1) = eta.subsVarsPartial([x,y],[xv1,yv1]);
            
                xv2 = obj.d.polygon.vx(obj.d.E(i,2));
                yv2 = obj.d.polygon.vy(obj.d.E(i,2));
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

                b1 = df1.subsVarsPartial ([x,y],[xv1,yv1]);
                b2 = df1.subsVarsPartial ([x,y],[xv2,yv2]);
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
                xv = obj.d.polygon.vx(obj.d.V(i));
                yv = obj.d.polygon.vy(obj.d.V(i));
                etaV(i) = eta.subsVarsPartial([x,y],[xv,yv]);
            end
            
        end

    end

    methods % conjugate
        function obj = conjugate (obj)
          for i=1:size(obj.envf,1)
              conjf = obj.conjugateFunction(i);
              %conjd = obj.envd(i).conjugate;
          end
        end

        function conjf = conjugateFunction (obj,i)
            vars = obj.f.getVars
            if obj.envExpr(i).type == 1

              t = sym('t');
              s1 = sym('s1');
              s2 = sym('s2');
              
              cpsi2 = obj.envExpr(i).vpsi2.getLinearCoeffs (vars);
             
              cpsi1 = obj.envExpr(i).vpsi1.getLinearCoeffs (vars);
              
              cpsi0 = obj.envExpr(i).vpsi0.getLinearCoeffs (vars);
              
              vs1 = s1 - (2*cpsi1(1)*t - cpsi2(1)*t^2 + cpsi0(1));
              vs2 = s2 - (2*cpsi1(2)*t - cpsi2(2)*t^2 + cpsi0(2));

              vt = solve (cpsi2(2)*vs1 - cpsi2(1)*vs2, t );
              crs = subs(vs1,t, vt);    


              
                
              



              
            end
            
            NCV = obj.getNormalConeVertex(i, s1, s2)
            NCE = obj.getNormalConeEdge(i, s1, s2)

            dualVars = [s1,s2];
        
            subdV = obj.getSubdiffVertex (i, NCV, dualVars)
            
        end

        function subdV = getSubdiffVertex (obj, i, NCV, dualVars)
            subdV = sym(zeros(obj.envd(i).nv,2));
            vars = obj.f.getVars
            drx1 = obj.envf(i).dfdx(vars(1));
            drx2 = obj.envf(i).dfdx(vars(2)) ;

            [ldrx1,limdrx1] = obj.limits (i, drx1, vars)
            [ldrx2,limdrx2] = obj.limits (i, drx2, vars)
            for j = 1:obj.envd(i).nv
              if ~ldrx1(j)
                  continue;
              end
              if ~ldrx2(j)
                  continue;
              end
              
              f = functionF(NCV(j,1));
              coef = f.getLinearCoeffs (dualVars);
              if (coef(2) == 0)
                subdV(j,1) = NCV(j,1)-limdrx1(j);
              elseif (coef(2) < 0)
                m = double(diff(NCV(j,1),dualVars(1)));
                c = yIntercept(m, [limdrx1(j),limdrx2(j)]);
                subdV(j,1) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -double(diff(NCV(j,1),dualVars(1)));
                c = yIntercept(m, [limdrx1(j),limdrx2(j)]);
                subdV(j,1) = dualVars(2) - m*dualVars(1) - c;
              end 

              f = functionF(NCV(j,1));
              coef = f.getLinearCoeffs (dualVars);
              if (coef(2) == 0)
                subdV(j,2) = NCV(j,2)-limdrx1(j);
              elseif (coef(2) < 0)
                m = double(diff(NCV(j,2),dualVars(1)));
                c = yIntercept(m, [limdrx1(j),limdrx2(j)]);
                subdV(j,2) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -double(diff(NCV(j,2),dualVars(1)));
                c = yIntercept(m, [limdrx1(j),limdrx2(j)]);
                subdV(j,2) = dualVars(2) - m*dualVars(1) - c;
              end 

              

              
             
            end

        end

        
        function [ldrx1,limdrx1] = limits (obj, i, drx1, vars)
            vars2 = [vars(2),vars(1)];
            for j = 1: obj.envd(i).nv
                l1 = drx1.limit(vars,[obj.envd(i).vx(j),obj.envd(i).vy(j)]);
                l2 = drx1.limit(vars,[obj.envd(i).vy(j),obj.envd(i).vx(j)]);
                if (l1 == l2)
                    ldrx1(j) = true;
                    limdrx1(j) = double(l1.f);
                else
                    ldrx1(j) = false;
                end
            end
        end 
            
        
        function NC = getNormalConeVertex(obj, i, s1, s2)
            
            NC = sym(zeros(obj.envd(i).nv,2));
            vars = obj.f.getVars;
            meanx = sum(obj.envd(i).vx)/obj.envd(i).nv;
            meany = sum(obj.envd(i).vy)/obj.envd(i).nv;
            
             for j = 1: obj.envd(i).nv-1
                slope = obj.envd(i).slope(j,j+1);
                pslope = -1/slope;
                if pslope ~= inf
                    q = obj.envd(i).yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.envd(i).vx(j);
                end
                if subs(eq,[s1,s2],[meanx,meany]) < 0
                    eq = -eq;
                end
                NC(j,1) = eq;
                if pslope ~= inf
                    q = obj.envd(i).yIntercept (j+1,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.envd(i).vx(j+1);
                end
                if subs(eq,[s1,s2],[meanx,meany]) < 0
                    eq = -eq;
                end
                
                NC(j+1,2) = eq;
                 
             end
             j = obj.envd(i).nv;
             slope = obj.envd(i).slope(j,1);
             pslope = -1/slope;
             if pslope ~= inf
               q = obj.envd(i).yIntercept (j,pslope);
               eq = s2 - pslope*s1 - q;
             else
               eq = s1 - obj.envd(i).vx(j);
             end
             if subs(eq,[s1,s2],[meanx,meany]) < 0
                    eq = -eq;
                end
                
             NC(j,1) = eq;
             if pslope ~= inf
               q = obj.envd(i).yIntercept (1,pslope);
               eq = s2 - pslope*s1 - q;
             else
               eq = s1 - obj.envd(i).vx(1);
             end
             if subs(eq,[s1,s2],[meanx,meany]) < 0
                    eq = -eq;
                end
                
             NC(1,2) = eq;
                
        end

    
    function NC = getNormalConeEdge(obj, i, s1, s2)
            
            NC = sym(zeros(obj.envd(i).nv,2));
            vars = obj.f.getVars;
            meanx = sum(obj.envd(i).vx)/obj.envd(i).nv;
            meany = sum(obj.envd(i).vy)/obj.envd(i).nv;
            
             for j = 1: obj.envd(i).nv-1
                slope = obj.envd(i).slope(j,j+1);
                q = obj.envd(i).yIntercept (j,slope);
                % get edge no
                edge = vars(2)-slope*vars(1)-q;
                if subs(edge,vars,[meanx,meany]) > 0
                    edge = -edge;
                end
                
                for k = 1: obj.envd(i).nv
                    if (obj.envd(i).ineqs(k).f == edge)
                        break;
                    end
                end
                
                
                %%
                pslope = -1/slope;
                if pslope ~= inf
                    q = obj.envd(i).yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.envd(i).vx(j);
                end
                if subs(eq,[s1,s2],[meanx,meany]) > 0
                    eq = -eq;
                end
                NC(k,1) = eq;
                if pslope ~= inf
                    q = obj.envd(i).yIntercept (j+1,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.envd(i).vx(j+1);
                end
                if subs(eq,[s1,s2],[meanx,meany]) > 0
                    eq = -eq;
                end
                
                NC(k,2) = eq;
                 
             end
             j = obj.envd(i).nv;
             
             slope = obj.envd(i).slope(j,1);
             edge = vars(2)-slope*vars(1)-q;
                if subs(edge,vars,[meanx,meany]) > 0
                    edge = -edge;
                end
                
                for k = 1: obj.envd(i).nv
                    if (obj.envd(i).ineqs(k).f == edge)
                        break;
                    end
                end
                
                   
             pslope = -1/slope;
             if pslope ~= inf
               q = obj.envd(i).yIntercept (j,pslope);
               eq = s2 - pslope*s1 - q;
             else
               eq = s1 - obj.envd(i).vx(j);
             end
             if subs(eq,[s1,s2],[meanx,meany]) < 0
                    eq = -eq;
                end
                
             NC(k,1) = eq;
             if pslope ~= inf
               q = obj.envd(i).yIntercept (1,pslope);
               eq = s2 - pslope*s1 - q;
             else
               eq = s1 - obj.envd(i).vx(1);
             end
             if subs(eq,[s1,s2],[meanx,meany]) < 0
                    eq = -eq;
                end
                
             NC(k,2) = eq;
                
        end

    end
    methods % max


        function obj = maxOfConjugate (obj)
          
        end

        % fix for multiple intersections
        % remove outside polytope ineqs
        function obj = maxEnvelope (obj, vars)
            envdT = [];
            envfT = [];
            enveT = [];
            for i = 1:size(obj.envd,2)
                l(i)=1;
            end
            for i = 1:size(obj.envd,2)
                for j = i+1:size(obj.envd,2)
                    if (obj.envd(i) == obj.envd(j))
                      %obj.envd(i).print  
                      envdT = [envdT,obj.envd(i)];
                      envfT = [envfT, pointwise_max(obj.envf(i), obj.envf(j), obj.d.polygon.vx, obj.d.polygon.vy, obj.envd(j).ineqs, vars)]     ;
                      enveT = [enveT,obj.envExpr(i)];
                      l(i) = 0;
                      l(j) = 0;
                    end
                end
            end
            for i = 1:size(obj.envd,2)
              if (l(i) == 1)
                    envdT = [envdT,obj.envd(i)];
                    enveT = [enveT,obj.envExpr(i)];
                    envfT = [envfT, obj.envf(i)];     
              end
            end
            %envfT
            obj.envf = envfT;
            obj.envd = envdT;
            obj.envExpr = enveT;
            return

        end

        function obj = removeNMax (obj, li, vars)
            d = [];
            for i = 1:size(obj.envd,2)
                if (i == li) 
                    continue;
                end
                p = pointwise_max(obj.envf(i), obj.envf(li), obj.d.polygon.vx, obj.d.polygon.vy, obj.envd(i).ineqs, vars);
                if (p == obj.envf(li))
                    d = [d,i];
                end
            end
            if size(d) == 0
                d = [d,li];
            end
            obj.envf(d) = [];
            obj.envd(d) = [];
        end

        % diff to be added
        function obj = maxEnvelopeIntersect (obj, vars)
            envdT = [];
            envfT = [];
            for i = 1:size(obj.envd,2)
                l(i)=1;
            end
            %return
            for i = 1:size(obj.envd,2)
                for j = i+1:size(obj.envd,2)
                    %if (obj.envd(i) == obj.envd(j))
                    d = intersection(obj.envd(i),obj.envd(j));
                    if isempty(d)
                        continue;
                    end
                    
                    %d0 = d.simplify (vars,obj.envd(i));
                    d0 = d.simplify (vars,obj.d.polygon);
                    disp("intersection")
                    obj.d.polygon.print
                    d.print
                    d0.print
                    if isempty(d0)
                        continue;
                    end
                    %obj.envd(i).print
                    %obj.envd(j).print
                    %d.print
                    envdT = [envdT,d0];
                    envfT = [envfT, pointwise_max(obj.envf(i), obj.envf(j), obj.d.polygon.vx, obj.d.polygon.vy, d0.ineqs, vars)]     ;
                    l(i) = 0;
                    l(j) = 0;
                    

                    %setDifference = simplify(obj.envd(i) & ~d)
                    obj.envd(i).print
                    d0.print
                    d1 = obj.envd(i) - d0
                    d1.print
                    d1 = d1.simplify (vars,obj.d.polygon);
                    d1.print
                    %continue
                    envdT = [envdT,d1];
                    envfT = [envfT, obj.envf(i)]     ;
                    
                    
                    obj.envd(j).print
                    d0.print
                    d1 = obj.envd(j) - d0
                    d1.print
                    d1 = d1.simplify (vars,obj.d.polygon);
                    d1.print
                    %d1 = simplify(obj.envd(j) - d0)
                    envdT = [envdT,d1];
                    envfT = [envfT, obj.envf(j)]     ;
                    

                    %end
                end
              
              
            end
            for i = 1:size(obj.envd,2)
              if (l(i) == 1)
                    envdT = [envdT,obj.envd(i)];
                    envfT = [envfT, obj.envf(i)];     
              end
            end
            %return
            
            %envfT
            obj.envf = envfT;
            obj.envd = envdT;
            return

        end




        function f = max(obj,i,j,x,y)
            for i = 1:obj.d.polygon.nv
                if (subs(obj.envd(i).f,[x,y],[obj.d.polygon.vx(i),obj.d.polygon.vy(i)]))    % change this to vertexindomain
                    f0 = obj.envf(i)-obj.envf(j)
                    
                end
            end
            f = functionF
        end

    
    end

end