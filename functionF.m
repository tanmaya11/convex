classdef functionF
    properties (Access=private)
        % to be re written restricting functions to be in terms of defined
        % variables only
        nv; 
        vars; 
        num = sym('num');
        den = sym('den');
        symf = sym('f')   % for sake of substitution
    end    
    properties  
        f = sym('f') ;
    end

    methods

        
        function obj = functionF(num, den)
            if nargin == 0
              obj.num=0;
              obj.den=1;
            elseif nargin == 1
              obj.num=num;  
              obj.den=1;  
            elseif nargin == 2
              obj.num=num;
              obj.den=den;
            end
            obj.f= obj.num/obj.den;
            obj.vars = symvar(obj.f);
            obj.nv = size(obj.vars,2);
        end
        
        function num = getNum(obj)
            num = obj.num;
        end   
        function den = getDen(obj)
            den = obj.den;
        end   
        
        function print(obj)
            disp(obj.f);
        end
        
        
        function vars = getVars(obj)
            vars = obj.vars;
        end

        function f = subsF (obj,x,y,xv,yv)
            x = xv;
            y = yv;
            f = subs(obj.f);
        end    

           
        
        function f = subsVarsPartial (obj,vars,varVals)
            f0 = simplify(subs(obj.f, vars, varVals));
            f = functionF(f0);
            
        end    

        function d = double (obj)
            c = obj.f.coeffs();
            if (size(c) == 0)
                d = 0;
                return
            end
            if (size(c) ~= 1)
                disp("Error in double in functionF");
                return
            end
            d = c(1);
        end

        function f = dfdx (obj,x)
            f = functionF(diff(obj.f,x));
        end    

        function f = solve (obj,x)
            f = solve(obj.f,x);
        end    
        
        function f = plus(obj1,obj2)
            f = functionF(obj1.f + obj2.f);
        end


        function f = minus(obj1,obj2)
            f = functionF(obj1.f - obj2.f);
        end


        function f = unaryminus(obj)
            f = functionF(-obj.f);
        end

        function res = eq(obj1,obj2)
            res = false;
            if (obj1.f==obj2.f)
                res = true;
            end 
        end
        
        function res = le(obj1,obj2)
            res = false;
            if (obj1.f<=obj2.f)
                res = true;
            end 
        end
        
        function res = ge(obj1,obj2)
            res = false;
            if (obj1.f>=obj2.f)
                res = true;
            end 
        end
        
        %obj2 double
        function res = lt(obj1,obj2)
            res = false;
            if (obj1.f<obj2)
                res = true;
            end 
        end
        
        function res = gt(obj1,obj2)
            res = false;
            if (obj1.f>obj2)
                res = true;
            end 
        end
        
        
        
        function f = removeDenominator (obj,x,y)
            %f = obj.num;
            %num = f.num
            cy = [];
            cx = coeffs(obj.num,x);
            for i = 1:size(cx,2)
                g  = cx(i);
                cy = [cy,coeffs(cx(i),y)];
                
            end
            cz=[];
            for i = 1:size(cy,2)
              cz = [cz,1/cy(i)];
            end
            if (size(cz)>0)
                mult= lcm(cz);
            end
            f = obj.f* abs(mult);
        end
        function printL (l)
            disp("Printing list")
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).print;
                end
            end
        end

        % temp code - needs to be fixed
        function max_func = pointwise_max(objf, objg, vx, vy, ineqs, ineqs2)
          % POINTWISE_MAX computes the pointwise maximum of two convex functions f and g
          % and returns a function handle to the resulting maximum function.


          
          % Define a function handle to the maximum function
          f = objf.f;
          g = objg.f;
          
          vars = objf.vars;
          max_func = @(vars) max(f(vars), g(vars));
          minx = min(vx);
          maxx = max(vx);
          miny = min(vy);
          maxy = max(vy);
          n = 0;
          step = 10;
          Z = [];
          tol = 1.0d-6;
          for i = step*minx:step*maxx
              for j = step*miny:step*maxy
                  xi = i/step;
                  xj = j/step;
                  l = true;
                  for k = 1:size(ineqs,2)
                      %ineqs(k)
                      l = l&(double(subs(ineqs(k).f,vars,[xi,xj]))<=0);
                  end
                  for k = 1:size(ineqs2,2)
                      %ineqs2(k)
                      l = l&(double(subs(ineqs2(k).f,vars,[xi,xj]))<0);
                  end
                  if (l)
                    n = n+1;
                    lf(n) = double(subs(f,vars,[xi,xj])) >= double(subs(g,vars,[xi,xj]));
                    le(n) = abs(double(subs(f,vars,[xi,xj])) - double(subs(g,vars,[xi,xj]))) <= tol;
                    %Z(n) = max(subs(f,vars,[xi,xj]),subs(g,vars,[xi,xj]));
                    %D(n,1)=xi;
                    %D(n,2)=xj;
                  end 
              end
          end
          if (all(lf) == 0)
              max_func = objg;
          else
              max_func = objf;
          end
          
        end
        function l = isZero(obj)
            l = (obj.f==0);
        end

        function l = isNegativeSqr(obj,z)
            l = false;
            z0 = simplify(obj.f);
            c = coeffs(obj.f);
            if c(end) > 0
                return
            end
            z1 = obj.solve(z);
            s = z1(1);
            for i = 2:size(z1,2)
                if (s ~= z1(i)) 
                    return;
                end
               
            end
             l = true;
        end
             
% not working 
            function l = isSubset (obj1, obj2)
                obj1.f <= 0
                obj2.f <= 0
                (obj1.f<=0) <= (obj2.f<=0)
              l = isAlways((obj1.f<=0) <= (obj2.f<=0))
            end

            
        

        
    end

end