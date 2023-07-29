classdef convexExpr
    properties
        type;
        expr = functionF.empty();
        psi0=sym('psi0');
        psi1=sym('psi1');
        psi2=sym('psi2');
        zeta=sym('zeta');
        vpsi0 = functionF.empty();
        vpsi1 = functionF.empty();
        vpsi2 = functionF.empty();
        vzeta = functionF.empty();
    end

    methods
        function obj = convexExpr(type, s0, s1, s2, dl)
            obj.type=type;
            %if nargin == 1
            %    obj.vpsi0=functionF(0);  
            %  obj.vpsi1=functionF(0);
            %  obj.vpsi2=functionF(1);
            %elseif nargin == 2
            %  obj.vpsi0=functionF(s0);  
            %  obj.vpsi1=functionF(0);
            %  obj.vpsi2=functionF(1);
            %elseif nargin == 3
            %  obj.vpsi0=functionF(s0);  
            %  obj.vpsi1=functionF(s1);
            %  obj.vpsi2=functionF(1);
            %else
            %  disp("Error in convexExpr")
            %end
            d = sym('d');
            
            if type == 1
                obj.expr = obj.psi1^2 /obj.psi2 + obj.psi0;
                obj.vpsi0=functionF(s0);  
                obj.vpsi1=functionF(s1);
                obj.vpsi2=functionF(s2);
            elseif type ==2
                %obj.expr= -d^2*obj.psi2 + 2* d* obj.psi1 + obj.psi0;
                obj.expr = obj.psi1^2 / obj.zeta + obj.psi0;
                obj.vpsi0=functionF(s0);  
                obj.vpsi1=functionF(s1);
                obj.zeta =functionF(s2);
            
            elseif type ==3
                obj.expr= obj.zeta * obj.psi1 + obj.psi0;
                obj.vpsi0=functionF(s0);  
                obj.vpsi1=functionF(s1);
                obj.zeta =functionF(s2);
            elseif type ==4
                % zeta = dl, du
                obj.expr= - obj.zeta^2 * obj.psi2 + 2 * obj.zeta * obj.psi1 + obj.psi0;
                obj.vpsi0=functionF(s0);  
                obj.vpsi1=functionF(s1);
                obj.vpsi2=functionF(s2);
                obj.zeta =functionF(dl);
            else
                disp("Error in convexExpr")
                %obj.expr= obj.psi1^2 /obj.psi2 + obj.psi0;
                % psi2 constant quad-quad case
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

            disp("zeta") 
            
            obj.zeta

        end
    end
end