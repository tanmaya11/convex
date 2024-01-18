classdef functionF
    properties (Access=private)
        % to be re written restricting functions to be in terms of defined
        % variables only
        nv; 
        vars; 
        num = sym('num');
        den = sym('den');
    end    
    properties  
        f = sym('f') ;
    end

   
    methods  % init 
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
            if nargin ~= 0
            
            obj.f= obj.num / obj.den;
            obj.vars = symvar(obj.f);
            obj.nv = size(obj.vars,2);
            end
        end
        
        function f = getF(obj)
            f = obj.f;
        end  
        function num = getNum(obj)
            num = obj.num;
        end   
        function den = getDen(obj)
            den = obj.den;
        end   
        
    end
    methods % display
        function print(obj)
          %fprintf(char(simplify(obj.f))); 
          %fprintf("\n")
          %obj.f
          if obj.isPolynomial
          [coef,terms] = coeffs(obj.f);
         
          for i=1:length(terms)
              if double(coef(i)) ~= 1
                fprintf(num2str(double(coef(i))));
              end
              if terms(i) ~= 1
                fprintf(char(terms(i)));
              end
              if i == length(terms)
                  break;
              end
              fprintf(" + ");
          end
          fprintf("\n")
          else
              fprintf(char(simplify(obj.f))); 
          fprintf("\n")
          
          end
        end

        function printLatexWB(obj)
          %fprintf(char(simplify(obj.f))); 
          %fprintf("\n")
          %obj.f
          if obj.isPolynomial
          [coef,terms] = coeffs(obj.f);
          
          for i=1:length(terms)
              if abs(double(coef(i))-1) > 1.0d-8
                  
                [n,d] = numden(coef(i));
                if  double(d) == 1
                  fprintf(num2str(abs(double(coef(i)))));
                else
                  fprintf("\\frac{" + num2str(abs(double(n)))+"}{"+ num2str(double(d))+"}");  
                end
              end
              if terms(i) ~= 1
                fprintf(char(terms(i)));
                
              end
              if i == length(terms)
                  break;
              end
              if double(coef(i+1)) > 0
                fprintf(" + ");
              else  
                fprintf(" - ");  
              end
          end
          else
              [n,d] = numden(obj.f);
              fprintf("\\frac{" + char(n)+"}{"+ char(d)+"}");  
              %fprintf(char(simplify(obj.f))); 
          
          end
        end

        function printLatex(obj)
          %fprintf(char(simplify(obj.f))); 
          %fprintf("\n")
          %obj.f
          %if obj.isPolynomial
          %[coef,terms] = coeffs(obj.f);
          
          fprintf("\\[");
          obj.printLatexWB;
          % for i=1:length(terms)
          %     if double(coef(i)) ~= 1
          % 
          %         [n,d] = numden(coef(i));
          %       if  double(d) == 1
          %         fprintf(num2str(double(coef(i))));
          %       else
          %         fprintf("\\frac{" + num2str(double(n))+"}{"+ num2str(double(d))+"}");  
          %       end
          %     end
          %     if terms(i) ~= 1
          %       fprintf(char(terms(i)));
          % 
          %     end
          %     if i == length(terms)
          %         break;
          %     end
          %     fprintf(" + ");
          % end
          % fprintf("\\]\n")
          % else
          %     [n,d] = numden(obj.f);
          %     fprintf("\\[\\frac{" + char(n)+"}{"+ char(d)+"}");  
          %     %fprintf(char(simplify(obj.f))); 
          fprintf("\\]\n")
          
         % end
        end

        function fprint(obj, uNo)
          fprintf(uNo, char(simplify(obj.f))); 
          fprintf(uNo, "\n")
        end

        
        function plot3d(obj, limits)
            if limits(1) == limits(2)
                limits(1) = limits(1)-30;
                limits(2) = limits(2)+30;
            end
            if limits(3) == limits(4)
                limits(3) = limits(3)-30;
                limits(4) = limits(4)+30;
            end
            ezsurf(obj.f, limits);
        end

        function plot(obj, vars, limits)
            colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k'];
            polyvars = obj.vars;
            if  size(vars,2) == size(polyvars,2)
                vx = linspace(limits(1), limits(2), 100);
                fy = solve(obj.f==0,vars(2));
                vy = subs(fy,polyvars(1),vx);
                %[vx,vy] = getPoints (obj, polyvars, limits)
                plot(vx,vy);
                %fill(vx,vy,colors(1+mod(obj.getGlobalParameter,7)),'FaceAlpha',0.5);
            elseif ismember (vars(1), polyvars)
                fx = solve(obj.f==0,vars(1));
                vy = zeros(1,100);
                vx = double(fx) + vy;
                vy = linspace(limits(1), limits(2), 100);
                plot(vx,vy);
            elseif ismember (vars(2), polyvars)
                fy = solve(obj.f==0,vars(2));
                vx = zeros(1,100);
                vy = double(fy) + vx;
                vx = linspace(limits(1), limits(2), 100);
                plot(vx,vy);
            else
            end
            %obj.setGlobalParameter; 
        end

        function plotIneq (obj, limits)
            ezsurf(obj.f <=0, limits)
        end 

        
        function printIneq(obj)
            fprintf(char(obj.f));
            fprintf(" <= 0 \n")
        end

        function printIneqLatex(obj)
            fprintf("\\[")
            obj.printLatexWB
            %fprintf(char(obj.f));
            fprintf(" \\le 0")
            fprintf("\\] ")
            fprintf("\n")
        end

        function printIneqM(obj)
            fprintf(char(obj.f));
            fprintf(" <= 0 ")
        end

        function printL (l, first, last)

            if nargin == 1
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).f = simplifyFraction(l(i,j).f);
                    l(i,j).print;
                end
            end
            else
                for i = 1: size(l,1)
                for j = first: last
                    l(i,j).f = simplifyFraction(l(i,j).f);
                    l(i,j).print;
                end
            end
            
            end
        end

        function printLLatex (l, first, last)

            if nargin == 1
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).printLatex;
                end
            end
            else
                for i = 1: size(l,1)
                for j = first: last
                    l(i,j).printLatex;
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

        function printLIneqLatex (l)
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).printIneqLatex;
                end
            end
        end

         function printLIneqM (l)
            
             fprintf("{")
                for j = 1: size(l,2)-1
                    l(j).printIneqM;
                    fprintf(",");
                end
                l(size(l,2)).printIneqM;
                    fprintf("}");
            
        end
        function fprintLIneq (l,uNo)
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).fprint(uNo);
                end
            end
        end

        function plotLIneq (l, vars,limit)
            % change this
            %l2 = [];
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    %l(i,j)
                     
                    l(i,j).plot(vars,limit);
                    hold on;
             %       l2 = [l2,l(i,j).f<=0]
                end
            end
            %plot(intersect(l2))
        end

        function plot2 (l)
            figure;
            plot(intersect(l))
        end
    
    end
    
    methods % operations
        function f = plus(obj1,obj2)
            f = functionF(obj1.num*obj2.den + obj2.num*obj1.den, obj1.den*obj2.den);
        end
        function f = minus(obj1,obj2)
            f = functionF(obj1.num*obj2.den - obj2.num*obj1.den, obj1.den*obj2.den);
        end

        % to be disabled
        function f = unaryminus(obj)
            f = functionF(-obj.num,obj.den);
        end

        function f = uminus(obj)
            f = functionF(-obj.num,obj.den);
        end

        
        function f = mtimes(f1,f2)
            num = f1.num * f2.num;
            den = f1.den * f2.den;
            f = functionF(num,den);
            % fix this
            % check def of num den here
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

    methods % inquiry

        function l = isPolynomial(obj)
            l = obj.degreeDen == 0;
        end
        
        function l = isQuad(obj)
         
          if (obj.degreeDen ~= 0)
              l = false;
              return;
          end
          if (obj.degreeNum == 2)
              l = true;
          else
              l = false;
          end
        end

        
        
        function l = isLinear(obj)
         
          if (obj.degreeDen ~= 0)
              l = false;
              return;
          end
          if (obj.degreeNum == 1)
              l = true;
          else
              l = false;
          end
        end

        % fix this
        function l = isConst(obj)

            if size(obj.vars,2) > 0
                l = false;
                return
            end
            cn = coeffs(obj.num,obj.vars);
            cd = coeffs(obj.den,obj.vars);
            l = all(cn(2:end)==0) & all(cd(2:end)==0);
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
               c(3) = coeffs(obj.f,vars(1));
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
                  c(3) = ct(1);
                  c(1) = ct(2);
                  c(2) = 0;
                else
                  c(2) = 0;
                  c(3) = 0;
                  c(1) = ct(1);
                end
               else
                ct = coeffs(obj.f,vars(2));
                if (size(ct,2) == 2)
                  c(3) = ct(1);
                  c(1) = 0;
                  c(2) = ct(2);
                else
                  c(2) = ct(1);
                  c(3) = 0;
                  c(1) = 0;
                end
               
              end
           end

        end

        function [l,c] = quadterm (obj, x)
            
            [qc,qt] = coeffs(obj.f,x);
            l = false;
            
            for i = 1:size(qt,2)
                if isAlways (qt(i) == x^2)
                    l = true;
                    c = qc(i);
                    return
                end
            end
        end
        

        % not for rational functions
        function obj = normalize1 (obj)
            if obj.den ~= 1
                disp('Rational in normalize1')
            end
            %obj.f
            %obj.vars
            
          c = coeffs(obj.f,obj.vars);
          obj = functionF(simplify((1/abs(c(end)))*obj.num), obj.den);
          
        end

        function obj = normalize (obj,vars)
          c = obj.getLinearCoeffs (vars);
             
          if (c(2) == 0)
            % double * f not overloaded
            obj.f = (1/c(1)) * obj.f;
          else
             obj.f = (1/c(2)) * obj.f;
          end
        end

        function f = subsF (obj,vars,vals)
            %x = xv;
            %y = yv;
            if (subs(obj.den, vars, vals) == 0)
                if (subs(obj.num, vars, vals) == 0)
                  f = functionF(sym(nan),1);  
                elseif (subs(obj.num, vars, vals) > 0)    
                  f = functionF(sym(intmax),1);
                else
                  f = functionF(sym(-intmax),1);
                end  
                return;
            end
            
            f = functionF(subs(obj.f, vars, vals));
        end    

           
        function l = isZero(obj)
            %disp("isZero")
            %obj.vars
            %obj.f
            l = false;
            if obj.den ~= 1
                return
            end
            
            c = coeffs(obj.num,obj.vars);

            n = size(c,2);
            for i = 1:n
                if (abs(double(c(i))) > 1.0e-6)
                    return
                end
            end
            l = true;
            
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
        
        function res = lt(obj1,obj2)
            
            if isequal(class(obj2),'double')
                
                res = ltd(obj1,obj2);
            else
            res = false;
            if (obj1.f<obj2.f)
                res = true;
            end 
            end
        end
        
        function res = gt(obj1,obj2)
            if isequal(class(obj2),'double')
                res = gtd(obj1,obj2);
            else
            res = false;
            
            if (obj1.f>obj2.f)
                res = true;
            end 
            end
        end
        
        
        
        %obj2 double
        function res = ltd(obj1,obj2)
            res = false;
            if (obj1.f<obj2)
                res = true;
            end 
        end
        
        function res = gtd(obj1,obj2)
            res = false;
            if (obj1.f>obj2)
                res = true;
            end 
        end
        

        function [vx,vy] = solveF (f2)
            f1x = subs(obj.f, obj.vars,[x,y]);
            f2x = subs(f2.f, obj.vars,[x,y]);
            s = solve ([f1==0,f2==0],[x,y]);
            vx = s.x;
            vy = s.y;
        end
            
        
        function f = removeDenominator2 (obj)
            %f = obj.num;
            %num = f.num
           
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
        
        
% not working 
%             function l = isSubset (obj1, obj2)
%                 obj1.f <= 0
%                 obj2.f <= 0
%                 (obj1.f<=0) <= (obj2.f<=0)
%               l = isAlways((obj1.f<=0) <= (obj2.f<=0))
%             end

            function d = degreeNum(obj)
                if isequal(class(obj.num),'double')
                    d = 0;
                else
                  d = polynomialDegree(obj.num);
                end
            end

            function d = degreeDen(obj)
                if isequal(class(obj.den),'double')
                    d = 0;
                else
                  d = polynomialDegree(obj.den);
                end
            end


            
    %        function fs = filterZero(obj)
    %            fs = [];
    %            for i = 1:size(obj,2)
    %                if isZero(obj.f)
    %                    continue
    %                end
    %                fs = [fs,obj(i)]
    %            end
    %        end

    end 

%     methods  % conjugate
%         function obj = conjugateRational(obj)
%             
%         end
% 
%         function obj = conjugate(obj)
%             if polynomialDegree(obj.getDen) > 0
%                 obj = conjugateRational(obj)
%             else
%                 %obj = conjugateR(obj)
%             end
%             
%         end
%     end

     methods % ineqs

        
%         function [l,obj] = removeParallel(obj, vars, lprint)
%             rm = [];
%             l = false;
%                 
%             for i = 1:size(obj,2)
%               c1 = obj(i).getLinearCoeffs (vars);
%               slope1 = c1(1)/c1(2);
%               for j = i+1:size(obj,2)
%                 c2 = obj(j).getLinearCoeffs (vars);
%                 slope2 = c2(1)/c2(2);
%                 %if nargin == 3
%                     %if lprint
%                     %disp('remove parallel')
%                     %c1
%                     %c2
%                     %end
%                 %end
%                 if (slope1 == slope2 & sign(c1(1)) == sign(c2(1)))
%                     if (c1(2) == 0)
%                         if ((c1(3)/c1(1) >= c2(3)/c2(1)) & c1(1)>0)
%                             rm = [rm,j];
%                         else
%                             rm = [rm,i];
%                         end
%                     else
%                         if (c1(3)/c1(2) >= c2(3)/c2(2))
%                             rm = [rm,j];
%                         else
%                             rm = [rm,i];
%                         end
%                     end
%                     l = false;
%                 elseif (slope1 == slope2 )
%                     if c1(2) ~= 0
%                     %    disp('h1')
%                     if c1(2) > 0
%                         l = c1(3)/c1(2) > - c2(3)/c2(2);
%                     else
%                     %    disp('h2')
%                         l = c2(3)/c2(2) < - c1(3)/c1(2);
%                     end
%                     else 
%                     if c1(1) > 0
%                         l = c1(3)/c1(1) > - c2(3)/c2(1);
%                     else
%                         l = c2(3)/c2(1) > - c1(3)/c1(1);
%                     end
%                     
%                     end
%                 end
%               end
%               
%             end
%             
%             obj(rm) = [];
%         end
% 
%         function l = antiParallelFeasible(obj,vars)
%             for i = 1:size(obj,2)
%               c1 = obj(i).getLinearCoeffs (vars);
%               slope1 = c1(1)/c1(2);
%               if slope1 == -inf
%                   slope1 = inf;
%               end
%               for j = i+1:size(obj,2)
%                 c2 = obj(j).getLinearCoeffs (vars);
%                 slope2 = c2(1)/c2(2);
%                 if slope2 == -inf
%                   slope2 = inf;
%                 end
%               
%                 if (slope1 == slope2 & sign(c1(1)) ~= sign(c2(1)))
%                         if (c1(3)/c1(1) + c2(3)/c2(1) > 0)
%                             l = false;
%                             return
%                         end
%                     
%                 end
%               end
%               
%             end
%             l = true;
%         
%         end
% 
        function [l,obj] = removeSum(obj, lprint)
            rm = [];
            l = false;
            for i = 1:size(obj,2)
              o1 = obj(i);
              for j = i+1:size(obj,2)
                o2 = obj(j) + o1;
                % removing = 0 also although only >0 are invalid
                if lprint
                    disp("o2 b4")
                    o2
                end
                
                o2 = simplify(o2.f < 0);
                if o2 == symfalse
                   l = true;
                  return
                end
                if lprint
                    disp("o2")
                    o2
                end
                for k = 1:size(obj,2)
                  o = simplify(obj(k).f<0);
                  if (o2 == o)
                      rm = [rm,k];
                  end
                end
              
              end
              
            end
            if lprint
            rm 
            end
            obj(rm) = [];
        end

    end

end