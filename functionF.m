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
            for i = 1:size(cx,1)
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
            f = obj.f*mult;
        end
        function printL (l)
            disp("Printing list")
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).print;
                end
            end
        end
        
    end

end