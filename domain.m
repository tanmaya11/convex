classdef domain
  properties
    nVertices {mustBeInteger}
    vx;
    vy;
    nE=0;
    E;
    mE;
    cE;
    nV=0;
    V;
    polygon;
    parabola; 
  end 

  methods

      
      function obj = domain(v)
          if nargin > 0
            obj.nVertices = size(v,1) ;
            [obj.vx,obj.vy] = poly2cw(v(:,1),v(:,2));
          end
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
            obj.mE(obj.nE) = m
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

      function m = slope (obj,i,j)
          m = (obj.vy(i)-obj.vy(j))/(obj.vx(i)-obj.vx(j));
      end
      
      function c = yIntercept (obj,i,m)
          c = obj.vy(i)-m*obj.vx(i);   
      end
      
  end
end