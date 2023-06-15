classdef convexExpr
    properties
        type;
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
            obj.type=type;
            if nargin == 1
                obj.vpsi0=functionF(0);  
              obj.vpsi1=functionF(0);
              obj.vpsi2=functionF(1);
            elseif nargin == 2
              obj.vpsi0=functionF(s0);  
              obj.vpsi1=functionF(0);
              obj.vpsi2=functionF(1);
            elseif nargin == 3
              obj.vpsi0=functionF(s0);  
              obj.vpsi1=functionF(s1);
              obj.vpsi2=functionF(1);
            else
              obj.vpsi0=functionF(s0);  
              obj.vpsi1=functionF(s1);
              obj.vpsi2=functionF(s2);
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