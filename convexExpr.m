classdef convexExpr
    properties
        expr = functionF.empty();
        psi0=sym('psi0');
        psi1=sym('psi1');
        psi2=sym('psi2');
        %zeta=sym('zeta');
        vpsi0 = functionF.empty();
        vpsi1 = functionF.empty();
        vpsi2 = functionF.empty();
    end

    methods
        function obj = convexExpr(type, s0, s1, s2)
            if nargin == 1
              obj.vpsi0=0;  
              obj.vpsi1=0;
              obj.vpsi2=1;
            elseif nargin == 2
              obj.vpsi0=s0;  
              obj.vpsi1=0;
              obj.vpsi2=1;
            elseif nargin == 3
              obj.vpsi0=s0;  
              obj.vpsi1=s1;
              obj.vpsi2=1;
            else
              obj.vpsi0=s0;  
              obj.vpsi1=s1;
              obj.vpsi2=s2;
            end
            d = sym('d');
            if type == 1
                obj.expr= obj.psi1^2 /obj.psi2 + obj.psi0;
            elseif type ==2
                obj.expr= -d^2*obj.psi2 + 2* d* obj.psi1 + obj.psi0;
            else
                obj.expr=obj.psi1 * obj.psi2 + obj.psi0;
            end
        end

        function print(obj)
            
            obj.expr
            fprintf("psi0 = ") 
            obj.vpsi0
            disp("psi1") 
            
            obj.vpsi1
            disp("psi2") 
            
            obj.vpsi2

        end
    end
end