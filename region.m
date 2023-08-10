classdef region
    % ineqs always stored as <= 0
    properties
        ineqs=functionF;
        lineqs;
       % not;
        % not impies not of union of ineqs

        % else its intersection of ineqs

        nv;
        vx;
        vy;
        vars;
    end

     methods
         function obj = region(fs, vars) %, not)
            % put checks for type of f and d
            %disp("region")
            %disp(nargin)

            %fs = fs.filterzero()
         
         %   if nargin == 2
         %     obj.not = false;
         %   elseif nargin == 3
         %     obj.not = not;
         %   else
         %     obj.not = false;
         %     return
         %   end

            if nargin == 0
                return
            end
            m = size(fs,1);
            n = size(fs,2);
            nineq = 0;
            for i = 1:m
              for j = 1:n
                f = functionF(fs(i,j));        
                if (~f.isZero)
                  nineq = nineq + 1;  
                  obj.ineqs(nineq) = f;
                end  
              end
            end
            obj.vars = vars;  
             
         end

         % add code for not
         function f = plus(obj1,obj2)
          %   if obj1.not | obj2.not
          %       disp("Implement not in plus")
          %   end
            l = []; 
            for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
            end
            %n = size(obj1.ineqs,2);
            for i = 1:size(obj2.ineqs,2)
              l = [l,obj2.ineqs(i).f];
            end 
            f = region(l,obj1.vars);
         end

         function f = minus(obj1,obj2)
%              if obj1.not | obj2.not
%                  disp("Implement not in minus")
%              end
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
            f = region(l,obj1.vars);
         end
         
         
         function l = eq(obj1,obj2)
             l = false;
             if (size(obj1.ineqs,2)~=size(obj2.ineqs,2))
                 return;
             end
             if obj1.nv ~= obj2.nv
                return
             end
             for i = 1:obj1.nv
                 l1(i) = false;
                 l2(i) = false;
             end

             for i = 1:obj1.nv
               for j = 1:obj2.nv
                   if (obj1.vx(i) == obj2.vx(j)) & (obj1.vy(i) == obj2.vy(j))
                       l1(i) = true;
                       l2(j) = true;
                       break;
                   end
               end
             end

             if (all(l1)==true )
                 l = true;
                 
                 
             else
                 l = false;
             end
             
             
             
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

         function m = slope (obj,i,j)
          m = (obj.vy(i)-obj.vy(j))/(obj.vx(i)-obj.vx(j));
         end

          function c = yIntercept (obj,i,m)
          c = obj.vy(i)-m*obj.vx(i);   
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
             disp("Variables")
             obj.vars
             disp(["nVertices = ", num2str(obj.nv)]);
             fprintf("vx =  ")
             fprintf("%d  ", obj.vx);
             fprintf("\n")
             fprintf("vy =  ")
             fprintf("%d  ", obj.vy);
             fprintf("\n\n")
%              disp("not")
%              obj.not
%              if obj.not
%                disp("Not of union of following ineqs")
%                obj.ineqs.printLIneq;
%              else
               disp("Intersection of following ineqs")
               obj.ineqs.printLIneq;
%              end
         %    disp(["Vertices = ", num2str(obj.nv)]);
         %    disp(obj.vx)
         %    disp(obj.vy)
         end

         function plot (obj)
             
             l1 = min(min(obj.vx),min(obj.vy));
             if l1 < -10
                 l1 = -10;
             end
             l2 = max(max(obj.vx),max(obj.vy));
             if l2 > 10
                 l2 = 10;
             end
             l1 = -6;
             l2 = 6;
           obj.ineqs.plotLIneq (obj.vars, [l1,l2])   ;

         end

         function plotL (obj)
             figure;
             size(obj)
             for i = 1:size(obj,2)
                 i
                 obj(i).plot
                 obj(i).plotRegion;
             end
         end


       
         function value = getGlobalParameter(obj)
 %           persistent icolor;
 %           value = icolor
             persistent icolor;
            if isempty (icolor)
                icolor=1
            else
                icolor=icolor+1
            end
            value = icolor
         end


         function [vx,vy] = plotRegion (obj)
            limitsx = [-6,6];
            limitsy = [-15,20];
            pts = 75;
            colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k'];
            stepx = (limitsx(2)-limitsx(1))/pts;
            stepy = (limitsy(2)-limitsy(1))/pts;
            n = 0;
            vx = [];
            vy = [];
            ci = limitsx(1);
            for i = 1:pts
                cj = limitsy(1);
                for j = 1:pts
                  
                  if obj.ptFeasible (obj.vars,[ci,cj]);
                      n = n+1;
                      vx(n) = ci;
                      vy(n) = cj;
                  end
                  cj = cj+stepy;    
                end
                ci = ci+stepx;    
            end
            c = colors(1+mod(obj.getGlobalParameter,7))
            if n == 0
                disp ('region not displayed')
            end
            fill(vx, vy, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
         end

         function plotByVertex (obj)
             %figure;
             limit = 6;
             for i = 1:obj.nv
                 vertices_ineq1(i, 1) = obj.vx(i);
                 if (vertices_ineq1(i, 1) > limit) 
                     vertices_ineq1(i, 1) = limit
                 end
                 if (vertices_ineq1(i, 1) < -limit) 
                     vertices_ineq1(i, 1) = -limit
                 end

                 vertices_ineq1(i, 2) = obj.vy(i);

                 if (vertices_ineq1(i, 2) > limit) 
                     vertices_ineq1(i, 2) = limit
                 end
                 if (vertices_ineq1(i, 2) < -limit) 
                     vertices_ineq1(i, 2) = -limit
                 end
                 
             end
             fill(vertices_ineq1(:, 1), vertices_ineq1(:, 2), 'b', 'FaceAlpha', 0.5);
         end
        
         % max over a region
         function [l, fmax, index] = maxArray (obj, f1, f2) 
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          for i = 1:size(fv1,2)
              sv1(i) = fv1(i).f;
              sv2(i) = fv2(i).f;
              if (isnan(sv1(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              elseif (isnan(sv2(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              end
          end
          l = true;
          if all(sv1 <= sv2)
              fmax = f2;
              index=2;
          elseif all(sv2 <= sv1)
              fmax = f1;
              index=1;
          else
              l = false;
              fmax = 0;
              index=0;
          end
         
        end
        

         % removes redundant ineq by finding points of intersection
         % check this code again
         function obj = simplify (obj, vars, obj2)
           rem = [];
           n = 0;
           x = sym('x');
           y = sym('y');
           %intersectingPts = []
           %intersectingEdges = []
           
           for i = 1:size(obj.ineqs,2)
             f1 = subs(obj.ineqs(i).f, vars,[x,y]);
             for j = i+1:size(obj.ineqs,2)
                 f2 = subs(obj.ineqs(j).f, vars,[x,y])  ;
                 s = solve ([f1==0,f2==0],[x,y]);
                 %disp('s.x')
                 %s.x
                 if isempty(s.x)
                     continue;
                 end
                 l =true;
                 if nargin > 2
                     l = obj2.ptFeasible (vars,[s.x,s.y]);
                 end
                 %if obj2.ptFeasible (vars,[s.x,s.y])
                 if l
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

           %n
           %intersectingPts{n}
           %intersectingEdges{n}

           % if intersecting point is not in feasible region mark it true
           for i = 1:size(intersectingEdges,2)
               %disp('intersecting point')
               %intersectingPts{i}
               %obj.ptFeasible (vars,intersectingPts{i})
               if obj.ptFeasible (vars,intersectingPts{i})
                 continue;
               end
               lp(i) = true;
               for j = 1:  size(intersectingEdges{i},2) 
                   ls(intersectingEdges{i}(j)) = true;
               end
           end
           %lp
           %ls

           % if valid point created by intersection of only 2 edges - mark
           % false
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
           %disp('lp2')
           %lp
           %ls
           
           % remove edges where point is to be removed
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
           %disp('lp3')
           %lp
           %ls
           
           for i = 1:size(intersectingEdges,2)
              if (lp(i))
                  continue
              end
              %if size(intersectingEdges{i},2) == 2
                lp(i)=true;
                for j = 1:  size(intersectingEdges{i},2)
                  ls(intersectingEdges{i}(j)) = false  ;
                end
              %else
              %  lp(i)=true;
              %  for j = 1:  2
              %    ls(intersectingEdges{i}(j)) = false  ;
              %  end  
              %end 
              intersectingEdges{i} = edges;
           end
           %disp('lp4')
           %lp
           %ls
           
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
                  lp(i)=true;
              end
           end
           %disp('lp5')
           %lp
           %ls
           
           % put code for lP false        
           %all(lp)
           %lp
           if  all(lp)==true
           for i = 1:size(ls,2)
               
              if (ls(i))
                rem = [rem,i];
              end
           end
           end
           %rem
           obj.ineqs(rem) = [];
         end

         function obj = simplify2 (obj, vars)
           rem = [];
             for i = 1:size(obj.ineqs,2)
             for j = i+1:size(obj.ineqs,2)
                 s = solve ([obj.ineqs(i).f==0,obj.ineqs(j).f==0],vars)
                 if isempty(s.x)
                     continue;
                 end
                 if obj2.ptFeasible (vars,[s.x,s.y])
                 
                 %obj.ptFeasible (vars,[s.x,s.y])
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
               %subs ([obj.ineqs(i).f],vars,point)
               for j = 1:size(point,1)
                %obj.ineqs(i).f
                %point(j,:)
                %subs ([obj.ineqs(i).f],vars,double(point(j,:)))
               if subs ([obj.ineqs(i).f],vars,double(point(j,:))) > 0
                   l = false;
                   return
               end
               end
           end  
         end

         function obj = removeDenominator(obj)
           for i = 1:size(obj.ineqs,2)
              
              obj.ineqs(i) = obj.ineqs(i).removeDenominator2;
           end 

         end

         function l = isFeasibleWBPts(obj)
             l = false;
             for i = 1:obj.nv
               if ~ obj.ptFeasible(obj.vars, [obj.vx(i),obj.vy(i)])
                   return
               end
             end
             l = true;
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
%          if obj1.not | obj2.not
%                  disp("Implement not in intersection")
%              end
         l = [];
           for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
           end 
         for i = 1:size(obj2.ineqs,2)
              l = [l,obj2.ineqs(i).f];
         end 
         obj = region(l, obj1.vars);
         if (isFeasible(obj))
             %disp('feasible')
             obj = obj.unique;
         else
             obj = region.empty;
         end

     end

     function [linter, objR] = intersection2(obj1, obj2, lprint)
         % split into linear and quad first
         %disp("Intersection2")

         
         linter = true;
         
         objR = obj1;
         obj = [obj1,obj2];
         nl = 0;
         nq = 0;
         vars = obj1.vars;
         
         for j = 1:size(obj,2)
           for i = 1:size(obj(j).ineqs,2)  
             if obj(j).ineqs(i).isLinear
                nl = nl+1; 
                l(nl) = obj(j).ineqs(i);
             else
                nq = nq+1; 
                q(nq) = obj(j).ineqs(i);
             end
           end
            
            
         end
      if lprint
          disp('nl nq')
          nl
          nq
      end

         if nl > 0
             
             %[notF,l] = l(1:nl).removeParallel(vars, lprint);
             %if notF
             %    linter = false  ;
                 %disp("Not feasible 1")
           
             %    return
             %end
             
             %if lprint
             %disp('removeParallel')
             %    l.printL
             %end
             [notF,l] = l.removeSum(lprint);
             if notF
                 linter = false  ;
                 %disp("Not feasible 2")
           
                 return
             end
             if lprint
             disp('removeSum')
             l.printL
             end
             l2 =[];
             for i = 1:size(l,2)
                 l2 = [l2,l(i).f];
             end
             lR = region(l2, obj1.vars);
             if lprint
                 disp('vertex simplify')
             lR = lR.getVertices();
             lR.isFeasibleWBPts
             lR.print;
             lR.plot;
             end
             %size(lR.ineqs,2)
             lR = lR.simplify (vars);
             %disp('aft simplify')
             %lR.print;
             
             if size(lR.ineqs,2) == 0
                 linter = false  ;
                % disp("Infeasible 3")
                 return
             end
             %disp("Linear")
             %lR.print
             
         end
         %nl
         %nq
         if nq > 0
         %disp("Quadratic")
         
          %   q
             q2 =[];
             for i = 1:size(q,2)
                 q2 = [q2,q(i).f];
             end
             
             qR = region(q2,obj1.vars);
           %  disp("lR")
           %  lR.print
            % disp("qR")
           %  qR.print;
           %  lR.print;
             objR = lR + qR;
         else    
             objR = lR;
         end
         %disp("Not in intersection2")
         %objR.not
         objR = objR.getVertices();
         if lprint
             disp("objR b4 simplify")
             
             objR.print
             disp("objR")
             
             objR = objR.simplify (objR.vars)
             objR.print
         end
         
         if objR.nv == 0
           linter=false  ;
         else
           linter = isFeasible(objR);
         end
     end
     
     function [linter,obj] = intersection3(obj1, obj2, lprint)
         % split into linear and quad first
         %disp("Intersection3")
         obj = region.empty();
       %  obj1.ineqs.printL
       %  obj2.ineqs.printL
         linter = false;
         n = 0;
         [linter, inter] = intersection2(obj1, obj2,lprint);
         if linter
             n = n + 1  ;
             obj(n) = inter;
         else
             obj(1) = obj1;     
         end  
     end
     
     function obj = normalizeEdge(obj)
         for i = 1:size(obj.ineqs,2)  
             c = getLinearCoeffs (obj.ineqs(i),obj.vars)
             
             if (c(2) == 0)
                 % double * f not overloaded
                 obj.ineqs(i).f = (1/c(1)) * obj.ineqs(i).f
             else
                 obj.ineqs(i).f = (1/c(2)) * obj.ineqs(i).f
             end
         end
     end

     % wont work for degree > 2
     function obj = getVertices(obj, lprint)
         
       obj.nv=0;
       t1 = sym('t1');
       t2 = sym('t2');
       varsTemp = [t1,t2];
       for i = 1:size(obj.ineqs,2)  
           f1 = obj.ineqs(i);
           f1 = f1.subsF (obj.vars,varsTemp);
           for j = i+1:size(obj.ineqs,2)  
               f2 = obj.ineqs(j);
               f2 = f2.subsF (obj.vars,varsTemp);
               s = solve ([f1.f==0,f2.f==0],varsTemp);
               
               if isempty(s)
                   continue;
               elseif isempty(s.t1)
                   continue;
               elseif isempty(s.t2)
                   continue;
               end
               if (obj.ptFeasible(obj.vars, [s.t1,s.t2]))
                   
                   obj.nv=obj.nv+1;
                   obj.vx(obj.nv) = s.t1;
                   obj.vy(obj.nv) = s.t2;
               end
               
           end
           
       end
       %[obj.vx,obj.vy] = poly2cw(obj.vx,obj.vy);

       % putting intmax for inf to avoid Nans 
       % intmax + intmax = intmax
       %disp('ptFeasile')
       %obj.print
       n = obj.nv;
       if obj.ptFeasible(obj.vars, [intmax,intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = intmax;
          obj.vy(obj.nv) = intmax;
          %disp("++")
       end
       if obj.ptFeasible(obj.vars, [intmax,-intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = intmax;
          obj.vy(obj.nv) = -intmax;
          %disp("+-")
       end
       if obj.ptFeasible(obj.vars, [-intmax,intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = -intmax;
          obj.vy(obj.nv) = intmax;
          %disp("-+")
       end
       if obj.ptFeasible(obj.vars, [-intmax,-intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = -intmax;
          obj.vy(obj.nv) = -intmax;
          %disp("--")
       end
       for i = 1:n
           if nargin == 2
               intmax,obj.vy(i)
               obj.ptFeasible(obj.vars, [intmax,obj.vy(i)])
           end
           if obj.ptFeasible(obj.vars, [obj.vx(i),intmax])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = obj.vx(i);
             obj.vy(obj.nv) = intmax;
           end
           if obj.ptFeasible(obj.vars, [obj.vx(i),-intmax])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = obj.vx(i);
             obj.vy(obj.nv) = -intmax;
           end
           if obj.ptFeasible(obj.vars, [intmax,obj.vy(i)])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = intmax;
             obj.vy(obj.nv) = obj.vy(i);
           end
           if obj.ptFeasible(obj.vars, [-intmax,obj.vy(i)])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = -intmax;
             obj.vy(obj.nv) = obj.vy(i);
           end
           
       end
       if obj.nv == 0
           return;
       end
       for i = 1:obj.nv
           V(i,1) = obj.vx(i);
           V(i,2) = obj.vy(i);
       end 
       %disp('unique vertices')
       %obj.vx,obj.vy
       V = unique(V,"rows");
       if obj.nv ~= size(V,1)
         obj.nv = size(V,1);
        % obj.vx
        % obj.vy
         obj.vx = V(:,1)';
         obj.vy = V(:,2)';
       end

     end

     function obj = conjugate(obj)
     end

     % If function is linear - check with lin ineqs 
                             % if coef are same try to resolve  
     function [l,less, fIneq, nf] = funcIneq (obj, f)
         % fix sign
         less = true;
       if ~f.isLinear
         fIneq = 0;
         l = false;
         nf = true
         return
       end
       c0 = f.getLinearCoeffs (obj.vars)
       for i =1:size(obj.ineqs,2)
           if ~obj.ineqs(i).isLinear
               continue;
           end
           c1 = obj.ineqs(i).getLinearCoeffs (obj.vars)
           if sign(c1(1)) ~= sign(c0(1))
               less = false;
               c1 = -c1;
           end
           if (c1(1) ~= c0(1))
               c1 = (c0(1)/c1(1)) * c1;
           end
           if (c1(2) == 0 & c0(2) == 0)
               fIneq = -c1(3);
               l = true;
               nf = false
               return
           elseif (c1(1)/c0(1) == c1(2)/c0(2))
               fIneq = -double(c1(3)/c1(2))+double(c0(3)/c0(2));
               l = true;
               nf = false
               return
           end
       end
       fIneq = 0;
       l = true;
       nf = true
         
     end
     
     % return function values at vertices
     function fv = funcVertices (obj, f)
         for i =1:obj.nv
             fv(i) = f.subsF(obj.vars,[obj.vx(i),obj.vy(i)]);
         end
     end


     function [l, maxf, index] = maximum(obj, f)
          %obj.print
          %f(1).print
          %f(2).print
          
          [l,  maxf, index] = obj.maxArray (f(1), f(2)) ;
     end

     
     % wont work cause of intersection vs union
     function [l,obj] = merge (obj, obj2)
         %obj.print;
         %obj2.print;
         l = false;
         n = 0;
         obj = obj + obj2;
         obj = obj.unique;
         %disp("unique")
         %obj.print;
         mark = [];
         for i =1:size(obj.ineqs,2)
           for j =1:size(obj.ineqs,2)
             if obj.ineqs(i) == obj.ineqs(j).unaryminus
                 l = true;
                 n = n + 1;
                 mark(n) = i;
                 n = n + 1;
                 mark(n) = j;
             end
           end
         
         end
         obj.ineqs(mark) = []; 
     end

     
     end

     
end