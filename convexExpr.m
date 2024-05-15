classdef convexExpr
    properties
        type;
        expr = symbolicFunction.empty();
        psi0=sym('psi0');
        psi1=sym('psi1');
        psi2=sym('psi2');
        zeta=sym('zeta');
        vpsi0 = symbolicFunction.empty();
        vpsi1 = symbolicFunction.empty();
        vpsi2 = symbolicFunction.empty();
        vzeta = symbolicFunction.empty();
    end

    methods  % testing
         function l = checkExpr1 (obj)
             x = sym('x');
             y = sym('y');
             l = false;

             if obj.type~= 1
                 return
             end
             f = 0;
             if ~isAlways(obj.vpsi0.f == f)
                 return
             end
             f = y/2;
             if ~isAlways(obj.vpsi1.f == f)
                 return
             end
             f = y/8 -x/8 + 1/4;
             if ~isAlways(obj.vpsi2.f == f)
                 return
             end
             
             l = true;
             return
             
             
         end

         function l = checkExpr2 (obj)
             x = sym('x');
             y = sym('y');
             l = false;

             if obj.type~= 3
                 return
             end
             
             f = x;
             if ~isAlways(obj.vpsi0.f == f)
                 return
             end
             
             f = y - 1;
             if ~isAlways(obj.vpsi1.f == f)
                 return
             end
             
             
             f = 2;
             if ~isAlways(obj.vzeta.f == f)
                 return
             end
             
             
             l = true;
             return
             
             
         end

         function l = checkExpr21 (obj)
             x = sym('x');
             y = sym('y');
             %obj.print
             l = false;

             if obj.type~= 1
                 return
             end
             f = (591*x)/176 - (197*y)/88 + 525/176;
             if ~isAlways(obj.vpsi0.f == f)
                 return
             end
             f = (17*x)/88 + (9*y)/44 - 5/88;
             if ~isAlways(obj.vpsi1.f == f)
                 return
             end
             f = x/44 - y/66 + 25/132;
             if ~isAlways(obj.vpsi2.f == f)
                 return
             end
             
             l = true;
             return
             
             
         end

         function l = checkExpr22 (obj)
             x = sym('x');
             y = sym('y');
            % obj.print
             l = false;

             if obj.type~= 1
                 return
             end
             f = (11*x)/4 - (11*y)/4 + 5/2;
             if ~isAlways(obj.vpsi0.f == f)
                 return
             end
             f = x/4 + y/4;
             if ~isAlways(obj.vpsi1.f == f)
                 return
             end
             f = x/36 - y/36 + 5/18;
             if ~isAlways(obj.vpsi2.f == f)
                 return
             end
             
             l = true;
             return
             
             
         end
     end
%  
    methods
        function obj = convexExpr(type, s0, s1, s2, dl)
            obj.type=type;
            %if nargin == 1
            %    obj.vpsi0=symbolicFunction(0);  
            %  obj.vpsi1=symbolicFunction(0);
            %  obj.vpsi2=symbolicFunction(1);
            %elseif nargin == 2
            %  obj.vpsi0=symbolicFunction(s0);  
            %  obj.vpsi1=symbolicFunction(0);
            %  obj.vpsi2=symbolicFunction(1);
            %elseif nargin == 3
            %  obj.vpsi0=symbolicFunction(s0);  
            %  obj.vpsi1=symbolicFunction(s1);
            %  obj.vpsi2=symbolicFunction(1);
            %else
            %  disp("Error in convexExpr")
            %end
            d = sym('d');
            
            if type == 1
                obj.expr = obj.psi1^2 /obj.psi2 + obj.psi0;
                obj.vpsi0=symbolicFunction(s0);  
                obj.vpsi1=symbolicFunction(s1);
                obj.vpsi2=symbolicFunction(s2);
            elseif type ==2  % disable this
                %obj.expr= -d^2*obj.psi2 + 2* d* obj.psi1 + obj.psi0;
                obj.expr = obj.psi1^2 / obj.zeta + obj.psi0;
                obj.vpsi0=symbolicFunction(s0);  
                obj.vpsi1=symbolicFunction(s1);
                obj.vzeta =symbolicFunction(s2);
            
            elseif type ==3
                obj.expr= obj.zeta * obj.psi1 + obj.psi0;
                obj.vpsi0=symbolicFunction(s0);  
                obj.vpsi1=symbolicFunction(s1);
                obj.vzeta =symbolicFunction(sym(s2));
            elseif type ==4   % disable this
                % zeta = dl, du
                obj.expr= - obj.zeta^2 * obj.psi2 + 2 * obj.zeta * obj.psi1 + obj.psi0;
                obj.vpsi0=symbolicFunction(s0);  
                obj.vpsi1=symbolicFunction(s1);
                obj.vpsi2=symbolicFunction(s2);
                obj.vzeta =symbolicFunction(dl);
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
            
            obj.vzeta

        end

        function printL(objL)
            
            for i = 1: size(objL,2)
                objL(i).print;
            end

        end
    end
end