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
         
         function print(obj)
             obj.ineqs.printL;
         end

         function obj2 = removeDenominator(obj)
           for i = 1:size(obj.ineqs,2)
              obj.ineqs(i) = obj.ineqs(i).removeDenominator2;
           end 

         end
     end
end