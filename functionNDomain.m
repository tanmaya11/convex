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
             disp("")
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

         function objR = maximumP(objL) %, f, r2)
           n = 0;
           for i = 1:size(objL,2)
             if size(objL(i).f,2) == 1
               n = n + 1;
               objR(n) = objL(i);
               continue;
             end
             [l, fmax, ind, lSing] = objL(i).d.maximum(objL(i).f);
             if lSing
                continue
             end
             if l
               n = n + 1;
               objR(n) = functionNDomain([fmax],objL(i).d);
               continue
             end  
                          
             ineqs = objL(i).d.splitmax3 (objL(i).f(1),objL(i).f(2));
             ineqs1 = sym.empty ;          
             for k = 1: size(objL(i).d.ineqs,2)
               ineqs1(k) = objL(i).d.ineqs(k).f;
             end
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(1).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1 = d1.simplifyOpenRegion;
             n = n + 1;
             objR(n) = functionNDomain([objL(i).f(1)],d1)
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(2).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1 = d1.simplifyOpenRegion;
             n = n + 1;
             objR(n) = functionNDomain([objL(i).f(2)],d1)
           end
           if n == 0
              return
           end
          % add merge here
            
         end 
     end
end 