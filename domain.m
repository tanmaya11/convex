classdef domain
  properties
    % vertex information  
    nVertices {mustBeInteger}
    vx;
    vy;
    ineqs=functionF;
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
            obj.nVertices = size(v,1) ;
            [obj.vx,obj.vy] = poly2cw(v(:,1),v(:,2));
          end
          %obj.nVertices
          %obj.vx
          %obj.vy
          obj = getEdges (obj);
          %obj.E
          %obj.mE
          %obj.V
          obj = getAllEdges (obj, x, y);
      end
      
      function vertex = getVertex(obj,i)
          vertex = [obj.vx(i), obj.vy(i)];
      end

      function obj = getEdges (obj)
          for i = 1:obj.nVertices-1
              l(i) = 0;
          end
          for i = 1:obj.nVertices-1
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
          
          m = obj.slope (obj.nVertices,1);
          if (m > 0 & m < inf)
            obj.nE = obj.nE+1;
            obj.E(obj.nE,1) =  obj.nVertices;
            obj.E(obj.nE,2) =  1;
            obj.mE(obj.nE) = m;
            obj.cE(obj.nE) = yIntercept (obj,1,m);
            l(1)=1;
            l(obj.nVertices) = 1;
          end

          
          for i = 1:obj.nVertices
              if (l(i) == 1) 
                  continue
              end    
              obj.nV = obj.nV+1;    
              obj.V(obj.nV) = i;
          end    
          
      end

      function obj = getAllEdges (obj, x, y)
        cx = mean(obj.vx);
        cy =  mean(obj.vy);
        for i = 1:obj.nVertices
          if (i == obj.nVertices) 
            m = obj.slope (i,1);
          else
            m = obj.slope (i,i+1);
          end
          if (m == inf)
            obj.ineqs(i) = functionF(x  - obj.vx(i));
          else
            obj.ineqs(i) = functionF(y - m*x - yIntercept (obj,i,m));
          end
          if obj.ineqs(i).subsVarsPartial([x,y],[cx,cy]) > 0
              %obj.ineqs(i) = -obj.ineqs(i)
              obj.ineqs(i) = obj.ineqs(i).unaryminus;
          end
        end
           
         
      end
      
      function m = slope (obj,i,j)
          m = (obj.vy(i)-obj.vy(j))/(obj.vx(i)-obj.vx(j));
      end
      
      function c = yIntercept (obj,i,m)
          c = obj.vy(i)-m*obj.vx(i);   
      end
      
  end
end