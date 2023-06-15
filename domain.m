
classdef domain
  properties
    % vertex information  
    %nVertices {mustBeInteger}
    %vx;
    %vy;
    
    % ineq representation
    polygon=region;
    % edge information
    nE=0; % number of convex edges
    E;    % index set 
    mE; % slope
    cE; % y intercept
    
    % remaining vertices
    nV=0;
    V;    % index set

  end 

  methods

      
      function obj = domain(v, x,y)
          if nargin > 0
            %obj.nVertices = size(v,1) ;
            %[obj.vx,obj.vy] = poly2cw(v(:,1),v(:,2));

            obj.polygon.nv=size(v,1) ;
            [obj.polygon.vx,obj.polygon.vy] = poly2cw(v(:,1),v(:,2));
          end
          %obj.nVertices
          %obj.vx
          %obj.vy
          obj = getEdges (obj);
          %obj.E
          %obj.mE
          %obj.cE
          %obj.V
          obj = getAllEdges (obj, x, y);
          %obj.polygon.print;
      end
      
      function vertex = getVertex(obj,i)
          vertex = [obj.polygon.vx(i), obj.polygon.vy(i)];
      end

      function obj = getEdges (obj)
          for i = 1:obj.polygon.nv-1
              l(i) = 0;
          end
          for i = 1:obj.polygon.nv-1
            m = obj.slope (i,i+1);
             if (m > 0 & m < inf)
                obj.nE = obj.nE+1;
                obj.E(obj.nE,1) =  i;
                obj.E(obj.nE,2) =  i+1;
                obj.mE(obj.nE) = m;
                obj.cE(obj.nE) = yIntercept (obj,i,m);
                l(i)=1;
                l(i+1)=1;
            end
          end 
          
          m = obj.slope (obj.polygon.nv,1);
          if (m > 0 & m < inf)
                
            obj.nE = obj.nE+1;
            obj.E(obj.nE,1) =  obj.polygon.nv ;
            obj.E(obj.nE,2) =  1;
            obj.mE(obj.nE) = m;
            obj.cE(obj.nE) = yIntercept (obj,1,m);
            l(1)=1;
            l(obj.polygon.nv) = 1;
          end

          
          for i = 1:obj.polygon.nv
              if (l(i) == 1) 
                  continue
              end    
              obj.nV = obj.nV+1;    
              obj.V(obj.nV) = i;
          end    
          
      end

      function obj = getAllEdges (obj, x, y)
        cx = mean(obj.polygon.vx);
        cy =  mean(obj.polygon.vy);
        for i = 1:obj.polygon.nv
          if (i == obj.polygon.nv) 
            m = obj.slope (i,1);
          else
            m = obj.slope (i,i+1);
          end
          if (m == inf | m == -inf)
            obj.polygon.ineqs(i) = functionF(x  - obj.polygon.vx(i));
          else
            obj.polygon.ineqs(i) = functionF(y - m*x - yIntercept (obj,i,m));
          end
          if obj.polygon.ineqs(i).subsVarsPartial([x,y],[cx,cy]) > 0
              %obj.ineqs(i) = -obj.ineqs(i)
              obj.polygon.ineqs(i) = obj.polygon.ineqs(i).unaryminus;
          end
        end
           
         
      end
      
      function m = slope (obj,i,j)
          m = (obj.polygon.vy(i)-obj.polygon.vy(j))/(obj.polygon.vx(i)-obj.polygon.vx(j));
          % change to
          %m = obj.polygon.slope(i,j);
      end
      
      function c = yIntercept (obj,i,m)
          c = obj.polygon.vy(i)-m*obj.polygon.vx(i);   
          % change to c = obj.yIntercept (i,m)
      end

      function print(obj)
        %disp(["nVertices = ", num2str(obj.polygon.nv)]);
        %fprintf("vx =  ")
        %fprintf("%d  ", obj.polygon.vx);
        %fprintf("\n")
        %fprintf("vy =  ")
        %fprintf("%d  ", obj.polygon.vy);
        %fprintf("\n")
        fprintf("Polygon Ineqs <= 0 \n")
        obj.polygon.print
        disp(["Number of Convex edges = ", num2str(obj.nE)])
        disp("Edges joining vertex numbers")
        disp(obj.E)
        fprintf("Slopes =  ")
        fprintf("%d  ", obj.mE);
        fprintf("\n")
        fprintf("y-intercepts =  ")
        fprintf("%d  ", obj.cE);
        fprintf("\n")
        disp(["Remaining vertices = ", num2str(obj.nV)])
        disp("Vertex number")
        disp(obj.V)
      end
        
      
      
  end
end