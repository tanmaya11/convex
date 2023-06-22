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

    
    methods  % init + display
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
            
            fprintf(char(simplify(obj.f)));
            fprintf("\n")
        end
        
        function printIneq(obj)
            fprintf(char(obj.f));
            fprintf(" <= 0 \n")
        end

        function printL (l, first, last)

            if nargin == 1
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).print;
                end
            end
            else
                for i = 1: size(l,1)
                for j = first: last
                    l(i,j).print;
                end
            end
            
            end
        end

        function printLIneq (l)
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).printIneq;
                end
            end
        end

    
    end
    
    methods % operations
        function f = plus(obj1,obj2)
            f = functionF(obj1.num*obj2.den + obj2.num*obj1.den, obj1.den*obj2.den);
        end
        function f = minus(obj1,obj2)
            f = functionF(obj1.num*obj2.den - obj2.num*obj1.den, obj1.den*obj2.den);
        end
        function f = unaryminus(obj)
            f = functionF(-obj.num,obj.den);
        end
    end

    methods % derivatives
         function f = dfdx (obj,x)
            
            f = functionF(simplify(diff(obj.f,x)));
         end 

         function f = limit (obj, vars, pt)
             f0 = obj.f;
             for i=1:size(vars,2)
               f0 = limit(f0,vars(i),pt(i));
             end
             f = functionF(f0);
         end

    end
    
    methods

        
        
        function vars = getVars(obj)
            vars = obj.vars;
        end

        % [x, y, const]
        function c = getLinearCoeffs (obj,vars)
           cvars = obj.getVars;
           if (isempty(cvars))
               c(1) = 0;
               c(2) = 0;
               c(3) = 0;
               return;
           end
           if size(cvars,2) == size(vars,2)
              ct = coeffs(obj.f,vars(1));
              c(1) = ct(2) ;
              ct = coeffs(ct(1),vars(2));
              if (size(ct,2) == 2)
                c(3) = ct(1);
              else
                  c(3) = 0;
              end
              ct = coeffs(obj.f,vars(2));
              c(2) = ct(2);
           elseif size(cvars,2) == 1
              if (cvars(1) == vars(1))
                ct = coeffs(obj.f,vars(1));
                if (size(ct,2) == 2)
                  c(2) = ct(1);
                  c(3) = ct(2);
                  c(1) = 0;
                else
                  c(2) = 0;
                  c(3) = 0;
                  c(1) = ct(1);
                end
               else
                ct = coeffs(obj.f,vars(2));
                if (size(ct,2) == 2)
                  c(2) = ct(1);
                  c(3) = 0;
                  c(1) = ct(2);
                else
                  c(2) = ct(1);
                  c(3) = 0;
                  c(1) = 0;
                end
               
              end
           end

        end

        function f = subsF (obj,vars,vals)
            %x = xv;
            %y = yv;
            f = functionF(subs(obj.f, vars, vals));
        end    

           
        function l = isZero(obj)
            l = isAlways(obj.f==0);
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

       
        function f = solve (obj,x)
            f = solve(obj.f,x);
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
        
        
        
        function f = removeDenominator2 (obj)
            %f = obj.num;
            %num = f.num
            %obj.vars
            cx = coeffs(obj.num,obj.vars);
            cz=[];
            for i = 1:size(cx,2)
              cz = [cz,1/cx(i)];
            end
            mult = 1;
            if (size(cz)>0)
                mult= lcm(cz);
            end
            f = obj.f* abs(mult);
            return

            x = obj.vars(1);
            y = obj.vars(2);
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
        % temp code - needs to be fixed
        %function max_func = pointwise_max(objf, objg, vx, vy, ineqs)
        function [index, max_func] = pointwise_max(objf, objg, vx, vy, ineqs, vars)
          % POINTWISE_MAX computes the pointwise maximum of two convex functions f and g
          % and returns a function handle to the resulting maximum function.


          
          % Define a function handle to the maximum function
          f = objf.f;
          g = objg.f;
          
          %varsf = objf.vars;  % Error when ony one variable
          %varsg = objg.vars;
          %if (size(varsf,2) < size(varsg,2))
          %    vars = varsg;
          %else
          %    vars = varsf;
          %end
          max_func = @(vars) max(f(vars), g(vars));
          minx = min(vx);
          maxx = max(vx);
          miny = min(vy);
          maxy = max(vy);
          
          
          n = 0;
          step = 5;
          Z = [];
          tol = 1.0d-6;
          for i = step*minx:step*maxx
              for j = step*miny:step*maxy
                  xi = i/step;
                  xj = j/step;
                  l = true;
                  for k = 1:size(ineqs,2)
                      %ineqs(k).f
                      %vars
                      %subs(ineqs(k).f,vars,[xi,xj])
                      l = l&(double(subs(ineqs(k).f,vars,[xi,xj]))<0);
                  end
                  %for k = 1:size(ineqs2,2)
                  %    %ineqs2(k)
                  %    l = l&(double(subs(ineqs2(k).f,vars,[xi,xj]))<0);
                  %end
                  if (l)
                    n = n+1;
                    lf(n) = double(subs(f,vars,[xi,xj])) >= double(subs(g,vars,[xi,xj]));
                    %if ~lf(n)
                    %    [xi,xj]
                    %    subs(f,vars,[xi,xj])
                    %    subs(g,vars,[xi,xj])
                    %end
                  end 
              end
          end
                  
          % put code for further division
          %lf
          %all(lf)
          if (all(lf) == 0)
              max_func = objg;
              index=2;
          else
              max_func = objf;
              index=1;
          end
          
        end
        %function l = isZero(obj)
        %    l = (obj.f==0);
        %end

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
        

        function l = isPolynomial(obj)
            obj.vars
            l = isPolynomial(obj.f, obj.vars);
        end
% not working 
            function l = isSubset (obj1, obj2)
                obj1.f <= 0
                obj2.f <= 0
                (obj1.f<=0) <= (obj2.f<=0)
              l = isAlways((obj1.f<=0) <= (obj2.f<=0))
            end

            % fix this
            function d = degreeNum(obj)
                
                d = polynomialDegree(obj.num);
            end

            

    end 

    methods  % conjugate
        function obj = conjugateRational(obj)
            
        end

        function obj = conjugate(obj)
            if polynomialDegree(obj.getDen) > 0
                obj = conjugateRational(obj)
            else
                %obj = conjugateR(obj)
            end
            
        end
    end

end