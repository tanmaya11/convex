classdef functionNDomain
    properties
      f = functionF.empty();
      d = region.empty();
    end

     methods
         function obj = functionNDomain(f, d)
             obj.f = f;
             obj.d = d;
         end

         function print (obj)
             obj.f.printL
             obj.d.print
         end

         function objL = mtimes (objL1, objL2)
             
             
             
             n = 0;
             for i = 1:size(objL1,2)
               for j = 1:size(objL2,2)
                 rf = objL1(i).d + objL2(j).d;
                 if isempty(rf)
                   continue
                 end
                 rf = rf.simplifyOpenRegion;
                 if isempty(rf)
                   continue
                 end
                 n = n + 1;
                 objL(n) = functionNDomain([objL1(i).f(1), objL2(j).f(1)],rf);
               end
             end
                      
         end
     end
end 