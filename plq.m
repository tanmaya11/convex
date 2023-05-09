classdef plq
    properties
        f=functionF.empty();
        d=domain.empty();
        envf=functionF.empty();
        envd=functionF.empty();
    end

    methods
        function obj = plq(d,f)
            if nargin > 0
              obj.f = [obj.f,f];
              obj.d = [obj.d,d];
            end 
        end

        function obj = maxEnvelope (obj,x,y)
            disp('chcek')
            for i = 1:size(obj.envd,2)
              obj.envd(i) = removeDenominator (obj.envd(i),x,y);
            end 
            for i = 1:size(obj.envd,2)
                for j = i+1:size(obj.envd,2)
                    if (obj.envd(i) == obj.envd(j))
                        i
                        j
                    end
                end
            end
        end
        function obj = convexEnvelope1 (obj,i1,x,y)
            a=sym('a');
            b=sym('b');
            obj.d(i1).V    
            
            [etaV, etaE, etaR] =  getEtaFunctions (obj,i1,x,y,a,b)
            % put a check that eta are only polynomials 

             % (vix, vjx) is the pair :  1 for edge else 0 
            [ix,jx,vix, vjx, ixd] = feasiblePairs (obj, i1,etaR, a,b)
            [envfs, envds] = solve (obj,i1, ix,jx,vix, vjx,ixd, etaV, etaE, etaR,a, b, x, y);
            for i = 1:size(envfs,2)
                obj.envf = [obj.envf, functionF(envfs(i))];
                obj.envd = [obj.envd, functionF(envds(i))];
            end
            
        end 

        function [envfs, envds] = solve (obj, i1, ix,jx,vix, vjx, ixd, etaV, etaE, etaR, a, b, x, y)
            envfs = [];
            envds = [];

            for i=1:size(ix,2)
                %disp("i")
                %disp(i)
                nc = 0;
                    
                if (vix(i)==1)
                    etah = etaE(ix(i),ixd(i));
                    degreeh = polynomialDegree(etah);
                    % should be mE(ix(i))
                    if size(obj.d(i1).mE,2) == 1
                      mh = obj.d(i1).mE;
                    else  
                      mh = obj.d(i1).mE(ix(i));
                    end
                    if (ixd(i)==1)
                      nc = nc + 1;
                      c(nc) = etaR(ix(i),1)-etaR(ix(i),3);
                    end
                    if (ixd(i)==3)
                      nc = nc + 1;
                      c(nc) = etaR(ix(i),3)-etaR(ix(i),1);
                    end
                    if (ixd(i)==2)
                        disp('ixd2 to be implemented')
                    end
                else
                    etah = etaV(ix(i));
                    degreeh = polynomialDegree(etah);
                end
                if (vjx(i)==1)
                    disp ('not implemented')
                    %etaw = etaE(jx(i));
                    %degreew = polynomialDegree(etaw);
                    %if size(obj.d(i1).mE,2) == 1
                    %  mh = obj.d(i1).mE;
                    %else  
                    %  mh = obj.d(i1).mE(jx(i));
                    %end  
                    
                    %nc = nc + 1;
                    %c(nc) = etaR(jx(i));
                else
                    etaw = etaV(jx(i));
                    degreew = polynomialDegree(etaw);
                end
               % disp("nc")
               % disp(nc)
               % disp(c(nc))
               % disp("eta")
                etah
                etaw
                degreeh
                degreew
                if (degreeh==1 & degreew==1)
                    %disp("here")
                    obj0 = etah + a*x+b*y;
                    bv = solve(etah==etaw,a);
                    obj0 = subs(obj0,a, bv);
                    %disp("factor")
                    objfacts = coeffs(obj0,b);
                    for k = 1:nc
                      c(k) = subs(c(k),a, bv);
                    end
                    %disp("etaV")
                    %size(etaV,2)
                    %etaV
                    for j=1:size(etaV,2)
                        %j
                        %etaV(j)
                        if ((etaV(j) == etah)  )
                            continue;
                        end
                        if ( (etaV(j) == etaw) )
                            continue;
                        end
                        etak = etaV(j);
                        nc = nc+1;
                        %disp("c")
                        c(nc) = etah - etak;  % <= 0   easier for substitution
                        c(nc) = subs(c(nc),a, bv);
                        
                    end
                    for j = 1:nc
                        bv(j) = solve(c(j),b);
                        if (subs(c(j),b, bv(j)+1) < 0)
                            bs(j)=0;
                        else
                            bs(j)=1;
                        end 
                    end
                    % stop for nc /= 2

                    if (bs(1) == 0) & (bs(2) == 1)
                        if (bv(1) <= bv(2))
                            envfs = [envfs, subs(obj0,b,bv(1))];
                            envds = [envds, objfacts(2) ];
                            envfs = [envfs, subs(obj0,b,bv(2))];
                            envds = [envds, -objfacts(2) ];
                            disp('feasible')
                        end
                    elseif (bs(1) == 1) & (bs(2) == 0)
                        if (bv(1) >= bv(2))
                            envfs = [envfs, subs(obj0,b,bv(1))];
                            envds = [envds, -objfacts(2) ];
                            envfs = [envfs, subs(obj0,b,bv(2))];
                            envds = [envds, objfacts(2) ];
                            
                            disp('feasible')
                        end
                    else
                        disp('infeasible')
                    end

                    
                end 
                
                if (degreeh==2 & degreew==1)
                    disp("one quad")
                    z = sym('z');
                    av = z - mh*b
                    eq = etah-etaw
                    eq = subs(eq,a,av)
                    bv = solve(eq,b)
                    nc
                    for j=1:size(etaV,2)
                        %j
                        %etaV(j)
                        if ((etaV(j) == etah)  )
                            continue;
                        end
                        if ( (etaV(j) == etaw) )
                            continue;
                        end
                        etak = etaV(j);
                        nc = nc+1;
                        %disp("c")
                        c(nc) = etah - etak;  % <= 0   easier for substitution
                        c(nc) = subs(c(nc),a, av)
                        c(nc) = subs(c(nc),b, bv)
                        
                    end
                    obj0 = etah + a*x+b*y
                    obj0 = subs(obj0,a,av)
                    obj0 = subs(obj0,b,bv)
                    deg = polynomialDegree(obj0,z)
                    
                    objfacts = coeffs(obj0,z)
                    size(objfacts)
                    if (deg == size(objfacts,2)+1)
                        psi0 = objfacts(1)
                        psi1 = objfacts(2)/2
                        psi2 = -objfacts(3)
                    elseif (deg == size(objfacts,2))
                        psi0=0
                        psi1 = objfacts(1)/2
                        psi2 = -objfacts(2)
                    else
                        psi0=0
                        psi1=0
                        psi2 = -objfacts(1)
                    end
                    positivePsi (obj,-psi2,x,y)
                    if positivePsi (obj,-psi2,x,y)==1
                       disp('+ to be implemented')
                    else
                       f0 = simplify(psi1^2/psi2+psi0)
                       r0 = -simplify(etaR(ix(i),3)*psi2-psi1)
                       envfs = [envfs, f0];
                       envds = [envds, r0];
                       
                       % multiply by lcm to simplify
                       f0 = -etaR(ix(i),3)^2*psi2 +2*etaR(ix(i),3)*psi1+psi0
                       r0 = simplify(etaR(ix(i),3)*psi2-psi1)
                       envfs = [envfs, f0];
                       envds = [envds, r0];
                            
                    end
                    cons = subs(c(1),a,av)
                    cons = subs(cons,b,bv)
                    %domain(cons)
                end
                        
            end     
        end

        function l=positivePsi (obj,p,x,y)
            p
            l = 0
            for i = 1, obj.d.nVertices
               if (subs(p,[x,y],[obj.d.vx(i),obj.d.vy(i)]) < 0)
                   return
               end
            end
            l = 1
        end

        
        function [ix,jx,vix, vjx, ixd] = feasiblePairs (obj, i1, etaR,a,b)
            n = 0;
            for i = 1:size(etaR,1)
                for j = 1:size(obj.d(i1).V,2)
                    vertex = [obj.d(i1).getVertex(obj.d(i1).V(j))];
                    s = subs(etaR(i,1),[a,b],vertex);
                    
                    if (s <= etaR(i,2))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 1;
                    end
                    if (s >= etaR(i,3))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 3;
                    end
                    if ((etaR(i,2)<=s)&(s<=etaR(i,3)))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 2;
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
                etaT1 = etaT;
                
                %etaE(3*(i-1)+1) = etaT;

                


                xv2 = obj.d(i1).vx(obj.d(i1).E(i,2));
                yv2 = obj.d(i1).vy(obj.d(i1).E(i,2));
                
                fT = obj.f(i1).subsF(x,y,xv2,yv2);
                etaT = subs(eta,[x,y,f],[xv2,yv2,fT]);
                etaT3 = etaT;
                %etaE(3*(i-1)+3) = etaT;



                m = obj.d(i1).mE(i);
                c = obj.d(i1).cE(i);
                yT = subs(edgey);
                
                etaT = subs(eta,y,yT);
                                %dfdx = diff (etaT,x)
                xp = solve(diff (etaT,x),x);

                
                etaM = simplify(subs(etaT,x,xp));
                etaT2 = etaM;

                %etaE(3*(i-1)+2) = etaM;
                
                
                f1 = obj.f(i1).subsF(x,y,x,yT);
                df1 = diff(f1,x);
                % store only <= in order to simplify later  (>= better but currently no multiplications involved)
                %etaR(3*(i-1)+1) = a + m*b - subs(df1,[x,y],[xv1,yv1]); % <= 0
                %etaR(3*(i-1)+2) =  piecewise( -a - m*b  + subs(df1,[x,y],[xv1,yv1])<=0,  (a + m*b - subs(df1,[x,y],[xv2,yv2])<=0));
                %etaR(3*(i-1)+2) =   subs(df1,[x,y],[xv1,yv1]) <= a + m*b <=   subs(df1,[x,y],[xv2,yv2])

                %etaR(3*(i-1)+2,2) = a + m*b -                 
                %etaR(3*(i-1)+3) = -a - m*b + subs(df1,[x,y],[xv2,yv2]); % <=0
                %etaE(i) = piecewise(etaT1,a + m*b <subs(df1,[x,y],[xv1,yv1]), etaT2, subs(df1,[x,y],[xv1,yv1]) <= a + m*b <=   subs(df1,[x,y],[xv2,yv2]), etaT3, a + m*b >   subs(df1,[x,y],[xv2,yv2]))
                
                %etaE(i) = piecewise(a + m*b <= subs(df1,[x,y],[xv1,yv1]),etaT1, subs(df1,[x,y],[xv1,yv1])<a+m*b<subs(df1,[x,y],[xv2,yv2]),etaT2, a+m*b>=subs(df1,[x,y],[xv2,yv2]), etaT3)
                
                etaE(i,1) = etaT1;
                etaE(i,2) = etaT2;
                etaE(i,3) = etaT3;
                etaR(i,1) = a+m*b;
                etaR(i,2) = subs(df1,[x,y],[xv1,yv1]);
                etaR(i,3) = subs(df1,[x,y],[xv2,yv2]);
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