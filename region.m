classdef region
    % ineqs always stored as <= 0
    properties
        ineqs=functionF;
    end

     methods
         function obj = region(fs)
            % put checks for type of f and d
            if nargin > 0
              for i = 1:size(fs,2)

                  obj.ineqs(i) = functionF(fs(i));
              end
            end 
         end

         function f = plus(obj1,obj2)
            l = []; 
            for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
            end
            %n = size(obj1.ineqs,2);
            for i = 1:size(obj2.ineqs,2)
              l = [l,obj2.ineqs(i).f];
            end 
            f = region(l);
         end
         
         function l = eq(obj1,obj2)
             l = false;
             if (size(obj1.ineqs,2)~=size(obj2.ineqs,2))
                 return;
             end
             for i = 1:size(obj1.ineqs,2)
               if (~(obj1.ineqs(i) == obj2.ineqs(i)))
                   return
               end
            end 
            
            l = true;
         end

         function obj = unique(obj)
             n = 0;
            duplicates = [];
            for i = 1:size(obj.ineqs,2)
                for j = i+1:size(obj.ineqs,2)
                    if (obj.ineqs(i) == obj.ineqs(j))
                        n = n + 1;
                        duplicates(n) = j;
                        lMark(j) = 1;
                    end
                end

            end 
            %duplicates
            obj.ineqs(duplicates) = [];
            %obj.ineqs
         end

         function print(obj)
             obj.ineqs.printL;
         end

         function obj = removeDenominator(obj)
           for i = 1:size(obj.ineqs,2)
              obj.ineqs(i) = obj.ineqs(i).removeDenominator2;
           end 

         end

         function l = isFeasible(obj)
             l = false;
             for i = 1:size(obj.ineqs,2)
               for j = i+1:size(obj.ineqs,2)
                 if ( obj.ineqs(i) == unaryminus(obj.ineqs(j)))
                     %disp('negative')
                     %obj.ineqs(i).f
                     %obj.ineqs(j).f
                    
                     return
                 end
                 if ( obj.ineqs(i) == obj.ineqs(j))
                     continue
                 end
                 
                 s = solve(obj.ineqs(i).f<=0,obj.ineqs(j).f<=0);
                 if isempty(s)
                     
                     %obj.ineqs(i).f
                     %obj.ineqs(j).f
                     
                     return;
                 end
                 %obj.ineqs(i).print
                 %obj.ineqs(j).print
               end 
 
             end 
             
             l = true;
         end

         % stupid way of doing this
     function obj = intersection(obj1, obj2)
         l = [];
           for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
           end 
         for i = 1:size(obj2.ineqs,2)
              l = [l,obj2.ineqs(i).f];
         end 
         obj = region(l);
         if (isFeasible(obj))
             disp('feasible')
             obj = obj.unique;
         else
             obj = region.empty;
         end

         end
    
     end

     
end