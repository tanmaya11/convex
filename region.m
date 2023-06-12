classdef region
    % ineqs always stored as <= 0
    properties
        ineqs=functionF;
        nv;
        vx;
        vy;
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

         function f = minus(obj1,obj2)
            l = []; 
            for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
            end
            %n = size(obj1.ineqs,2);
            l2 = [];
            for i = 1:size(obj2.ineqs,2)
              ladd = true;  
              for j = 1:size(obj1.ineqs,2)
                if (obj2.ineqs(i) == obj1.ineqs(j))
                    ladd = false;
                    break;
                    
                end

                
              end
              if ladd
                  l2 = [l2,obj2.ineqs(i).f];
              end
            end
            if (size(l2,2) ~= 0)
              mult = (-1) ^ size(l2,2);
              for i = 1: size(l2,2)
                mult = mult * l2(i);
              end
              l = [l,mult];
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

         function set = ithSet (obj, logical_indices)
           set = true;
           if all(logical_indices == false)
             set = false;
           end
                   
           for j = 1:size(logical_indices,2)
              if (logical_indices(j))
                set = set & obj.ineqs(j).f<=0;
              end
           end
            
         end
         
         function obj = removeSuperSet(obj, vars)
           % powerset indices
           %set_size=size(obj.ineqs,2);
           %logical_indices = mat2cell(double(dec2bin(0:2^set_size-1,set_size))==48,ones(2^set_size,1));
           logical_indices = supersetIndices(size(obj.ineqs,2))
           
           for i = 1:size(logical_indices,1)
              lS(i) = false  ;
           end
           disp('superset')
           size(logical_indices,1)
           %for i1 = 1:size(logical_indices,1)
             i1 = 1;  
             set1 = ithSet (obj, logical_indices{i1})

             s1 = solve(set1,vars)
             %if islogical(set1) 
             %    continue
             %end
             %if (~ set1) 
             %    continue
             %end
             %assume(set1);
               
             for i2 = 1:size(logical_indices,1)
                 if i1 == i2
                     continue
                 end
               set2 = ithSet (obj, logical_indices{i2})
               s2 = solve(set2,vars)
               if isequal(s1,s2)
                   disp('found')
               end
               continue;
               assume(set1);
               l2 = isAlways(set2)
               if ( l2)
                 assume(set2);
                 l1 = isAlways(set1)
               end
               
               continue
               
               if isAlways(set2)
                   lS(i2) = true;
                   disp('found')

                   logical_indices{i1}
                   logical_indices{i2}
              %     duplicates = [];
              %     for j = 1:size(logical_indices{i2},2)
              %         logical_indices{i1}(j)
              %         logical_indices{i2}(j)
              %         logical_indices{i1}(j)  ~= logical_indices{i2}(j)
              %         if (logical_indices{i1}(j)  ~= logical_indices{i2}(j) )
              %             duplicates = [duplicates,j]
              %         end
              %     end
              %     obj.ineqs(duplicates) = [];
                   
                   %return
               end
               %return
               
             %end
             %unassume(set1);

           end
           lS
         end

         function print(obj)
             obj.ineqs.printL;
             disp(["Vertices = ", num2str(obj.nv)]);
             disp(obj.vx)
             disp(obj.vy)
         end

         function obj = simplify (obj, vars, obj2)
           rem = [];
           n = 0;
           %intersectingPts = []
           %intersectingEdges = []
           
           for i = 1:size(obj.ineqs,2)
             for j = i+1:size(obj.ineqs,2)
                 s = solve ([obj.ineqs(i).f==0,obj.ineqs(j).f==0],vars);
                 if isempty(s.x)
                     continue;
                 end
                 if obj2.ptFeasible (vars,[s.x,s.y])
                     pointExists = false;
                     for k = 1:n
                         if (intersectingPts(k) == [s.x,s.y])
                             pointExists = true;
                             break;
                         end
                     end
                     if pointExists
                       intersectingEdges{k} = [intersectingEdges{k},[i,j]];
                       
                     else
                       n = n + 1;
                       intersectingPts{n} = [s.x,s.y];
                       intersectingEdges{n} = [i,j];
                     end
                 end
                 
             end 
           end 
           for i = 1:size(intersectingEdges,2)
               lp(i) = false;
           end
           for i = 1:size(obj.ineqs,2)
               ls(i) = false;
               
           end
           for i = 1:size(intersectingEdges,2)
               if obj.ptFeasible (vars,intersectingPts{i})
                 continue;
               end
               lp(i) = true;
               for j = 1:  size(intersectingEdges{i},2) 
                   ls(intersectingEdges{i}(j)) = true;
               end
           end
           for i = 1:size(intersectingEdges,2)
              if (lp(i))
                  continue
              end
              
              if size(intersectingEdges{i},2) == 2 
                 for j = 1:  size(intersectingEdges{i},2) 
                   ls(intersectingEdges{i}(j)) = false;
                 end
              end
           end
           
           for i = 1:size(intersectingEdges,2)
              if (lp(i))
                  continue
              end
              edges = [];
              for j = 1:  size(intersectingEdges{i},2)
                  if mod(j,2)==0
                      continue
                  end
                  if ~(        ls(intersectingEdges{i}(j))|ls(intersectingEdges{i}(j+1))  )  
                    edges = [edges,[intersectingEdges{i}(j),intersectingEdges{i}(j+1)]];
                  end
              end
              intersectingEdges{i} = edges;
           end
           
           % mark singletons and filter 
           for i = 1:size(obj.ineqs,2)
               ls(i) = true;
               
           end
           
           for i = 1:size(intersectingEdges,2)
              if (lp(i))
                  continue
              end
              if size(intersectingEdges{i},2) ~= 2
                  continue
              end
              lp(i)=true;
              for j = 1:  size(intersectingEdges{i},2)
                ls(intersectingEdges{i}(j)) = false  ;
              end
              intersectingEdges{i} = edges;
           end
           
           for i = 1:size(intersectingEdges,2)
              if (lp(i))
                  continue;
              end
              lM = false;
              for j = 1:  size(intersectingEdges{i},2)
                if mod(j,2)==0
                  continue;
                end
                if ls(intersectingEdges{i}(j)) & ls(intersectingEdges{i}(j+1))
                    lM=true;
                    break;
                end
              end
              if lM
                  lP(i)=true;
              end
           end
           
           % put code for lP false        
           
           for i = 1:size(ls,2)
               
              if (ls(i))
                rem = [rem,i];
              end
           end
           obj.ineqs(rem) = [];
         end

         function obj = simplify2 (obj, vars, obj2)
           rem = [];
             for i = 1:size(obj.ineqs,2)
             for j = i+1:size(obj.ineqs,2)
                 s = solve ([obj.ineqs(i).f==0,obj.ineqs(j).f==0],vars);
                 if isempty(s.x)
                     continue;
                 end
                 if obj2.ptFeasible (vars,[s.x,s.y])
                 if obj.ptFeasible (vars,[s.x,s.y])
                     continue;
                 end
                 rem = [rem,i];
                 rem = [rem,j];
                 end
                 
             end 
           end 
           
           obj.ineqs(rem) = [];
         end

         function l = ptFeasible(obj, vars, point)
           l = true;
           for i = 1:size(obj.ineqs,2)
               if subs ([obj.ineqs(i).f],vars,point) > 0
                   l = false;
                   return
               end
           end  
         end

         function obj = removeDenominator(obj)
           for i = 1:size(obj.ineqs,2)
              obj.ineqs(i) = obj.ineqs(i).removeDenominator2;
           end 

         end

         function l = isFeasible(obj)
             l = false;
             for i = 1:size(obj.ineqs,2)
                 %obj.ineqs(i).print
               for j = i+1:size(obj.ineqs,2)
                   %obj.ineqs(j).print
                 if ( obj.ineqs(i) == unaryminus(obj.ineqs(j)))
                 %if ( obj.ineqs(i) == -obj.ineqs(j))
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
    
     

     % wont work for degree > 2
     function obj = getVertices(obj,vars)
       obj.nv=0;
       for i = 1:size(obj.ineqs,2)  
           for j = i+1:size(obj.ineqs,2)  
               s = solve ([obj.ineqs(i).f==0,obj.ineqs(j).f==0],vars);
               if isempty(s)
                   continue;
               end
               if (obj.ptFeasible(vars, [s.x,s.y]))
                   disp('here')
                   obj.nv=obj.nv+1;
                   obj.vx(obj.nv) = s.x;
                   obj.vy(obj.nv) = s.y;
               end
               
           end
           
       end
       [obj.vx,obj.vy] = poly2cw(obj.vx,obj.vy);
       
     end
     end
     
end