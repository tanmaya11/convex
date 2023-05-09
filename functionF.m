classdef functionF
    properties (Access=private)
        % to be re written restricting functions to be in terms of defined
        % variables only
        nvar; 
        x; 
        num = sym('num');
        den = sym('den');
    end    
    properties  
        f = sym('f') ;
    end

    methods

        function setNum(obj,expr)
            obj.num = expr;
        end     
        function setDen(obj,expr)
            obj.den = expr;
        end     
        
        function obj = functionF0(n)
            if nargin == 0
                disp("No variables")
                return 
            end     
            if isinteger (n)
                obj.setVars(n);
            end
        end


        function obj = functionF(num, den)
            if nargin == 0
              obj.num=1;
              obj.den=1;
            elseif nargin == 1
              obj.num=num;  
              obj.den=1;  
            elseif nargin == 2
              obj.num=num;
              obj.den=den;
            end
            obj.f= obj.num/obj.den;

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
        function f = subsF (obj,x,y,xv,yv)
            x = xv;
            y = yv;
            f = subs(obj.f);
        end    

        function f = dfdx (obj,x)
            f = diff(obj.f,x);
        end    

        function f = minus(obj1,obj2)
            f = obj1.f - obj2.f;
        end

        function res = eq(obj1,obj2)
            res = false;
            if (obj1.f==obj2.f)
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
    end
end