classdef region
    % ineqs always stored as <= 0
    properties
        ineqs=symbolicFunction;
        nv;
        vx=sym.empty();
        vy=sym.empty();
        vars;
    end

    
%  57 methods
     methods

         % check when multiple regions are created by ineqs
         function obj = region(fs, vars) %, not)
        
            if nargin == 0
                return
            end
            m = size(fs,1);
            n = size(fs,2);
           
            nineq = 0;
            for i = 1:m
              for j = 1:n
                  f = symbolicFunction(fs(i,j));
                if (~f.isZero)
                  nineq = nineq + 1;  
                  obj.ineqs(nineq) = f;
                end  
              end
            end
            
            obj.vars = vars;  
            obj = obj.normalize1;
            obj = obj.unique;
            %obj.ineqs.printL
            obj = obj.getVertices  ;
            if obj.nv == 0
                obj = region.empty;
            end
            % put simplify in here
         end

         function [l, pos] = isVertexIrrational (obj)
             l = false;
             n = 0;
             pos = [];
             for i = 1:obj.nv
                 if symFunType(obj.vx(i)) == 'plus' | symFunType(obj.vy(i)) == 'plus'
                     l = true;
                     n = n + 1;
                     pos(n) = i;
                 end
             end
         end
         function obj = regionWPts(obj, vx, vy, x, y)
             n = 0;
             for i = 1: size(vx,2)-1
                 if vx(i) == vx(i+1)
                   n = n + 1;
                   ineq(n) = x-vx(i)  ;
                 else
                   m = (vy(i+1)-vy(i))/(vx(i+1)-vx(i));
                   c = vy(i)-m*vx(i);
                   n = n + 1;
                   ineq(n) = y-m*x-c
                 end 
             end
             if vx(1) == vx(end)
               n = n + 1;
               ineq(n) = x-vx(1)  
             else
               m = (vy(end)-vy(1))/(vx(end)-vx(1));
               c = vy(1)-m*vx(1);
               n = n + 1;
               ineq(n) = y-m*x-c
             end
             meanx = sum(vx)/3;
             meany = sum(vy)/3;
             for i = 1:n
                 if double(subs(ineq(i),[x,y],[meanx,meany])) > 0
                     ineq(i) = -ineq(i);
                 end

             end
             obj = region(ineq,[x,y]);
         end

         function f = plus(obj1,obj2)
            l = []; 
            for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
            end
            for i = 1:size(obj2.ineqs,2)
               lunique = true;
               for j = 1:size(obj1.ineqs,2)
                   %obj1.ineqs(j)
                   %obj2.ineqs(i)
                 if (obj1.ineqs(j) == obj2.ineqs(i))
                     lunique = false;
                     break;
                 end
                 if (obj1.ineqs(j).f == -obj2.ineqs(i).f)
                     f = region.empty;
                     return;
                 end
               end
               if lunique
                   l = [l,obj2.ineqs(i).f];
               end
            end 
            %l
            f = region(l,obj1.vars);
         end

         function obj = poly2orderUnbounded(obj)
           rad = cart2pol(obj.vx(1:obj.nv), obj.vy(1:obj.nv));
           radWrapped = mod(rad,2*pi);
           radWrapped(radWrapped==0 & rad>0) = 2*pi; 
           [~, sortIdx] = sort(radWrapped, 'descend'); 
           obj.vx(sortIdx) = obj.vx(1:obj.nv);
           obj.vy(sortIdx) = obj.vy(1:obj.nv);
           for i = 1:obj.nv
             edges = obj.getEdges(obj.vx(i),obj.vy(i));
             if size(edges,2) == 1
                 break;
             end
           end
           vx = obj.vx;
           vy = obj.vy;
           obj.vx(1:obj.nv-i+1) = vx(i:obj.nv); 
           obj.vy(1:obj.nv-i+1) = vy(i:obj.nv);
           obj.vx(obj.nv-i+2:obj.nv) = vx(1:i-1); 
           obj.vy(obj.nv-i+2:obj.nv) = vy(1:i-1); 
         end

         function obj = poly2order(obj)
            % obj.print
             vx(1) = obj.vx(1);
             vy(1) = obj.vy(1);
             for i = 1:obj.nv
               lineqs(i) = false ;
             end
             for i = 1:obj.nv-1
                 ineqs = obj.ineqThroughVertex (i);
                 for j = 1:2
                     if lineqs(j)
                         continue
                     end
                     break;
                 end
                 %ineqs(j)
                 lineqs(ineqs(j)) = true;
                 v = obj.getEndpoints (ineqs(j));
                 if (v(1,1) == vx(i) & v(1,2) == vy(i))
                     vx(i+1) = v(2,1);
                     vy(i+1) = v(2,2);
                 else
                     vx(i+1) = v(1,1);
                     vy(i+1) = v(1,2);
                 end


             end
             obj.vx = vx;
             obj.vy = vy;
             cx = sum(obj.vx) / obj.nv;
             cy = sum(obj.vy) / obj.nv;
%obj.print
             theta1 = atan((vy(1)-cy)/(vx(1)-cx));
             if theta1 < 0
                 theta1 = theta1 + pi;
             end
             theta2 = atan((vy(2)-cy)/(vx(2)-cx));
             if theta2 < 0
                 theta2 = theta2 + pi;
             end
             if theta1 > theta2
                 return;
             end
             for i = 1:obj.nv
                 obj.vx(i) = vx(obj.nv-i+1);
                 obj.vy(i) = vy(obj.nv-i+1);
             end
         end

         function ineqs = ineqThroughVertex (obj,j)
           ineqs = [];
           for i = 1:size(obj.ineqs,2)
               if abs(double(subs(obj.ineqs(i).f,obj.vars,[obj.vx(j),obj.vy(j)]))) <= 1.0e-6
                   ineqs = [ineqs,i];
               end
           end
         end

         function v = getEndpoints (obj, i)
             n = 0;
             for j = 1:obj.nv
               if abs(double(subs(obj.ineqs(i).f,obj.vars,[obj.vx(j),obj.vy(j)]))) <= 1.0e-6
                   n = n + 1;
                   v(n,1) = obj.vx(j);
                   v(n,2) = obj.vy(j);
               end
             end
            
             v = sortrows(v);
         end

         function n = getN(obj)
             n = 0;
             for i = 1: obj.nv
                 if abs(obj.vx(i)) == intmax
                     continue
                 end
                 if abs(obj.vy(i)) == intmax
                     continue
                 end
                 n = n + 1;
             end
         end

         function l = in(obj1,obj2)
             l = false;
             for i = 1:obj2.nv
                 if ~obj1.ptFeasible (obj1.vars,[obj2.vx(i),obj2.vy(i)])
                     return
                 end
             end
             l = true;
         end

         function f = minus(obj1,obj2)
             
             if (obj1 == obj2)
                 f = [region.empty];
                 return
             end     
            l = [];
            for i = 1:size(obj1.ineqs,2)
                v1(i,:,:) = obj1.getEndpoints(i);
                l = [l,obj1.ineqs(i).f];
            end
            l2 = [];
            rm = [];
            for i = 1:size(obj2.ineqs,2)
              v2(i,:,:) = obj2.getEndpoints(i);
              ladd = true;
              lsub = false;
              for j = 1:size(obj1.ineqs,2)
                %  i,j
                  %obj1.ineqs(j).f
                if (obj2.ineqs(i) == obj1.ineqs(j)) 
                    lv = true;
                    for i1 = 1:size(v1,2)
                        for i2 = 1:size(v1,3)
                            lv = lv & (v1(j,i1,i2) == v2(i,i1,i2));
                        end
                    end
                    if lv
                      rm = [rm,j];
                      ladd = false;
                    else 
                      ladd = false;  
                      lsub = true;  
                    end  
                %f = [f0.simplify];
                    
                    break;
                end
              end
              %ladd
              %obj2.ineqs(i) == obj1.ineqs(j)
              %obj2.ineqs(i).f
              if ladd
                  l2 = [l2,obj2.ineqs(i).f];
                 
              end
              % if lsub
              %     l2 = [l2,-obj2.ineqs(i).f];
              % 
              % end
            end
            
            %l2
            l(rm) = [];
            % l
            % l2
            for i = 1: size(l2,2)
                l = [l,-l2(i)];
            end
            % l
            % size(l)
            if size(l,2) <= 2
                f = [region.empty];
                return
            end
            % disp('here')
            f0 = region(l,obj1.vars);
            % f0.print
            %f0.getN
            if ~isempty (f0) & f0.getN > 2 
                f = [f0];
                return
            end
            % disp('here2')
            rL = region.empty();
            
            for i = 1:size(l,2)
                
                for j = i+1:size(l,2)
                    for k = j+1:size(l,2)
                        l0 = [l(i),l(j),l(k)];
                        f0 = region(l0,obj1.vars);
                        % i,j,k
                         % f0.print
                         
                        if isempty(f0)
                            continue;
                        end
                        if f0.getN < 3
                             continue;
                        end
                        if ~ obj1.in(f0)
                            continue
                        end
                        f0 = f0.simplify;
                        if f0.nv ~= size(f0.ineqs,2)
                            continue
                        end
                        rL = [rL,f0]; 
                    end

                end

            end
            % disp('triple')
            % for i = 1:size(rL,2)
            %     rL(i).print
            % end
            for i = 1:size(l,2)
                for j = i+1:size(l,2)
                    for k = j+1:size(l,2)
                       for il = k+1:size(l,2)
                        l0 = [l(i),l(j),l(k),l(il)];
                        f0 = region(l0,obj1.vars);
                        if isempty(f0)
                            continue;
                        end
                                   if f0.getN < 3
                             continue;
                        end
                        if ~ obj1.in(f0)
                            continue
                        end
                        f0 = f0.simplify;
                        if f0.nv ~= size(f0.ineqs,2)
                            continue
                        end
                        rL = [rL,f0]; 
                       end 
                    end

                end

            end

            for i = 1:size(l,2)
                for j = i+1:size(l,2)
                    for k = j+1:size(l,2)
                       for il = k+1:size(l,2)
                         for im = il+1:size(l,2)
               
                           l0 = [l(i),l(j),l(k),l(il),l(im)];
                           f0 = region(l0,obj1.vars);
                           if isempty(f0)
                             continue;
                           end
                           if f0.getN < 3
                             continue;
                           end
                           if ~ obj1.in(f0)
                            continue
                           end
                        
                           f0 = f0.simplify;
                        if f0.nv ~= size(f0.ineqs,2)
                            continue
                        end
                        rL = [rL,f0]; 
                         end
                       end 
                    end

                end

            end
            if (size(l,2) > 6)
                disp('minus to implement ')
                size(l,2)
                l
                obj1.print
                obj2.print
            end
            %rL

            if size(rL,2) > 0
                f = rL.uniqueL;
            else
                f = [];
                return
            end
            rm = [];
            for i = 1:size(f,2)
             if f(i) == obj1
                 rm = [rm,i];
             end
            end
            
            
            f(rm) = [];
             % disp('finak')
             % for i = 1:size(f,2)
             % f(i).print
             % end
            return

            f0 = region(l,obj1.vars);
            disp('in minus')
            l
            f0.print
            if isempty(f0)
                f = [f0];
                return
            end
            if f0.nv >= 3
              f = [f0.simplify];
              return
            else  
              f1 = f0.divideRegions(obj1);
              f = [];
              for i = 1: size(f1,2)
                  if f1(i).nv <= 2
                      continue
                  end
                  
                  f0 = f1(i).simplify;
                  if f0 == obj1
                      continue
                  end
                  f = [f,f0];
              end
            end

            
            
            
         end
         
%          function f = minus0(obj1,obj2)
% %              if obj1.not | obj2.not
% %                  disp("Implement not in minus")
% %              end
%             l = []; 
%             for i = 1:size(obj1.ineqs,2)
%               l = [l,obj1.ineqs(i).f];
%             end
%             %n = size(obj1.ineqs,2);
%             l2 = [];
%             for i = 1:size(obj2.ineqs,2)
%               ladd = true;  
%               for j = 1:size(obj1.ineqs,2)
%                 if (obj2.ineqs(i) == obj1.ineqs(j))
%                     ladd = false;
%                     break;
% 
%                 end
% 
% 
%               end
%               if ladd
%                   l2 = [l2,obj2.ineqs(i).f];
%               end
%             end
%             if (size(l2,2) ~= 0)
%               mult = (-1) ^ size(l2,2);
%               for i = 1: size(l2,2)
%                 mult = mult * l2(i);
%               end
%               l = [l,mult];
%             end
%             f = region(l,obj1.vars);
%          end
         
         
         function l = eqVertices(obj1,obj2)
             l = false;
             %obj1.nv
             if obj1.nv ~= obj2.nv
                return
             end
             for i = 1:obj1.nv
                 l1(i) = false;
                 l2(i) = false;
             end

             for i = 1:obj1.nv
               for j = 1:obj2.nv
                   if l2(j)
                       continue
                   end
                   if (abs(obj1.vx(i) - obj2.vx(j))<1.0d-15) & (abs(obj1.vy(i) - obj2.vy(j)) < 1.0d-15)
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

         function l = eq(obj1,obj2)
             l = false;
             
             if (size(obj1.ineqs,2)~=size(obj2.ineqs,2))
                 return;
             end
             
             if ~ obj1.eqVertices(obj2)
                 return
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
           
            obj.ineqs(duplicates) = [];
         end


         function obj = uniqueL(obj)
             n = 0;
             
            duplicates = [];
            for i = 1:size(obj,2)
                for j = i+1:size(obj,2)
                    if (obj(i) == obj(j))
                        n = n + 1;
                        duplicates(n) = j;
                        
                    end
                end

            end
           
            obj(duplicates) = [];
         end
         function m = slope (obj,i,j)
          m = (obj.vy(i)-obj.vy(j))/(obj.vx(i)-obj.vx(j));
         end

         function m = slopeIneq(obj,i,pt)
             if obj.ineqs(i).isLinear
                 c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                 
                 if c(2) == 0
                     m = inf;
                 else
                     m = -c(1)/c(2);
                 end
             else
                  vars =  obj.vars;
                  drx1 = obj.ineqs(i).dfdx(vars(1));
                  drx2 = obj.ineqs(i).dfdx(vars(2));
             %     subs(drx2.f,vars,pt)
                  if abs(subs(drx2.f,vars,pt)) < 1.0d-6
                      
                      m = intmax;
                  else
                      m0 = -drx1.f/drx2.f;
                      
                      m = subs(m0,vars,pt);
                  end
                 % m = -intmax;
                 % disp("to be implemented in slopeIneq")

                 
             end
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
         
       

         function print(obj)
             if isempty(obj)
                 disp("Empty region")
                 return
             end
             disp("Variables")
             obj.vars
             disp(["nVertices = ", num2str(obj.nv)]);
             fprintf("vx =  ")
             fprintf("%f  ", obj.vx);
             fprintf("\n")
             fprintf("vy =  ")
             fprintf("%f  ", obj.vy);
             fprintf("\n\n")
               disp("Intersection of following ineqs")
               obj.ineqs.printLIneq;
         end

         function printLatex(obj)
             disp("Variables")
             fprintf("\\[")
             for i = 1: size(obj.vars,2)
               fprintf(char(obj.vars(i)) );
               if i == size(obj.vars,2)
                 break;
               end
               fprintf(",");
             end
             n = 0;
             for i = 1: obj.nv
                 if (abs(obj.vx(i)) == intmax)
                     continue
                 end
                 if (abs(obj.vy(i)) == intmax)
                     continue
                 end
                 n = n + 1;
             end
             fprintf("\\]\n\\[\\text{Number of Vertices = }")
             fprintf(num2str(n));
             
             fprintf("\\]\n\\[v =  ")
             for i = 1: obj.nv
                 if (abs(obj.vx(i)) == intmax)
                     continue
                 end
                 if (abs(obj.vy(i)) == intmax)
                     continue
                 end
             fprintf("(" + num2str(obj.vx(i)) + ","+num2str(obj.vy(i))+")");
             if i == obj.nv
                 break;
             end
             fprintf(",");
             end
             fprintf("\\]\n")
             
               disp("Intersection of following ineqs")
               obj.ineqs.printLIneqLatex;
         end

         
         function printMaple(obj)
             obj = obj.subsF;
             
             obj.ineqs.printLIneqM;
             
         end

         function fprint(obj, uNo)
             fprintf(uNo, num2str(obj.nv)+"\n");
             for i = 1:obj.nv
               fprintf(uNo, num2str(obj.vx(i)) + "  " + num2str(obj.vy(i)) + "\n")  
             end
             obj.ineqs.fprintLIneq(uNo);
         end

          function plot (obj)
             
             l1 = min(min(obj.vx),min(obj.vy));
             if l1 < -25
                 l1 = -25;
             end
             l2 = max(max(obj.vx),max(obj.vy));
             if l2 > 20
                 l2 = 20;
             end
             l1 = -25;
             l2 = 20;
           obj.ineqs.plotLIneq (obj.vars, [l1,l2])   ;

         end


         function [vx,vy] = plotRegionC (obj, textR, c)
            limitsx = [-25,17];
            limitsy = [-15,20];
            
            pts = 75;
            colors = ['b', 'r', 'g', 'm', 'c', 'y'];
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
            %c = colors(1+mod(obj.getGlobalParameter,6))
            if n == 0
                disp ('region not displayed')
            end
            avx = sum(vx)/n;
            avy = sum(vy)/n;
            m = 0;
            for i = 1:n
                if (abs(vx(i) -avx)<0.1 & abs(vy(i) -avy)<0.1)
                    continue
                end    
                m = m+1;
                vx1(m) = vx(i);
                vy1(m) = vy(i);
            end
            %text = "R";
            text(avx,avy,textR,'FontSize',12, 'FontWeight', 'bold','Color', 'k')
            if m==0
                fill(vx, vy, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            else
                fill(vx1, vy1, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            end
            
         end

         function [vx,vy] = plotRegion (obj, textR)
            limitsx = [-6,6];
            limitsy = [-15,20];
            limitsx = [-25,17];
            limitsy = [-6,11];
            
            pts = 75;
            colors = ['b', 'r', 'g', 'm', 'c', 'y'];
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
            c = colors(1+mod(obj.getGlobalParameter,6))
            if n == 0
                obj.print 
                obj.ptFeasible (obj.vars,[2,-9])
                disp ('region not displayed')
            end
            avx = sum(vx)/n;
            avy = sum(vy)/n;
            m = 0;
            for i = 1:n
                if (abs(vx(i) -avx)<0.1 & abs(vy(i) -avy)<0.1)
                    continue
                end    
                m = m+1;
                vx1(m) = vx(i);
                vy1(m) = vy(i);
            end
            %text = "R";
            text(avx,avy,textR,'FontSize',12, 'FontWeight', 'bold','Color', 'k')
            if m==0
                fill(vx, vy, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            else
                fill(vx1, vy1, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            end
            
         end

         % function plotByVertex (obj)
         %     %figure;
         %     limit = 6;
         %     for i = 1:obj.nv
         %         vertices_ineq1(i, 1) = obj.vx(i);
         %         if (vertices_ineq1(i, 1) > limit) 
         %             vertices_ineq1(i, 1) = limit
         %         end
         %         if (vertices_ineq1(i, 1) < -limit) 
         %             vertices_ineq1(i, 1) = -limit
         %         end
         % 
         %         vertices_ineq1(i, 2) = obj.vy(i);
         % 
         %         if (vertices_ineq1(i, 2) > limit) 
         %             vertices_ineq1(i, 2) = limit
         %         end
         %         if (vertices_ineq1(i, 2) < -limit) 
         %             vertices_ineq1(i, 2) = -limit
         %         end
         % 
         %     end
         %     fill(vertices_ineq1(:, 1), vertices_ineq1(:, 2), 'b', 'FaceAlpha', 0.5);
         % end
        
         function obj = subsF(obj)
             x= sym('x');
             y= sym('y');
             for i = 1:size(obj.ineqs,2)
                 obj.ineqs(i) = obj.ineqs(i).subsF(obj.vars,[x,y]);
             end
         end
         function m = slopes (obj)
           n = 0 ;  
           for i = 1:size(obj.ineqs,2)
              if obj.ineqs(i).isLinear
                  c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                  n = n + 1;
                  m(n) = -c(1)/c(2);
              end
           end
         end

         
         function m = slopeAtVertex (obj, pi, pt)
            % pt
           n = 0 ;  
           %pi
           %pt
           %obj.ineqs.printL
           for j = 1:size(pi,2)
              i = pi(j); 
              if obj.ineqs(i).isLinear
                  c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                  n = n + 1;
                  m(n) = -c(1)/c(2);
              else
                  %check if point lies on curve
                  %
            %      obj.ineqs(i).print
                  vars =  obj.vars;
                  drx1 = obj.ineqs(i).dfdx(vars(1));
                  drx2 = obj.ineqs(i).dfdx(vars(2));
             %     subs(drx2.f,vars,pt)
                  if abs(subs(drx2.f,vars,pt)) < 1.0d-6
                      n = n + 1;
                      m(n) = intmax;
                  else
                      m0 = -drx1.f/drx2.f;
                      n = n + 1;
                      m(n) = subs(m0,vars,pt);
                  end
              end
           end
         end

         
         function m = slopes2 (obj)
           n = 0 ;  
           % l = false;
           for i = 1:size(obj.ineqs,2)
              if obj.ineqs(i).isLinear
                  c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                  n = n + 1;
                  m(n) = -c(1)/c(2);
              else
                  vars =  obj.vars;
                  
                  for j = 1:obj.nv
                      if abs(obj.ineqs(i).subsF(vars,[obj.vx(j),obj.vy(j)]).f) < 1.0d-8
                          pt = [obj.vx(j),obj.vy(j)];
                          break;
                      end
                  end
                  obj.ineqs(i)
                  drx1 = obj.ineqs(i).dfdx(vars(1));
                  drx2 = obj.ineqs(i).dfdx(vars(2));
                  %drx1
                  %drx2
                  % check sign
                  m0 = drx1.f/drx2.f;
                  n = n + 1;
                  %disp("Slope of tangent")
                  %pt
                  % l = true;
                  if abs(subs(drx2.f,vars,pt)) < 1.0d-6
                      m(n) = intmax;
                  else
                      m(n) = subs(m0,vars,pt);
                  end
                  

              end
              % if l
              %      disp("Slope of tangent")
              %      pt
              %      m
              % end

           end
         end

         % max over a region
         function [l, fmax, index, lsing] = maxArray (obj, f1, f2) 
          lsing = false;  
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          m = obj.slopes2();
          for i = 1:size(fv1,2)
              sv1(i) = fv1(i).f;
              sv2(i) = fv2(i).f;
          end
          %sv1
          %sv2
          l = true;
          if all(abs(double(sv1 - sv2))< 1.0d-14)
              if size(sv1) == 1
                  fmax = f1;
                  index=1;
              end
              nx = 0;
              for iv = 1:obj.nv
                  if (abs(obj.vx(iv))==intmax|abs(obj.vy(iv))==intmax)
                      continue
                  end
                  nx = nx + 1;
                  vx0(nx) = obj.vx(iv);
                  vy0(nx) = obj.vy(iv);
              end
              
              if nx > 1
                sx = mean(vx0);
                sy = mean(vy0);
                [l,fmax,index] = maxFromPt(obj, [sx,sy], [f1,f2]);
                if l 
                    return
                end
              end
              % use slope mid pt to get directions
              for i = 1:size(m,2)
                
                for j = i+1: size(m,2)
                  if (abs(m(i))~= inf) & (abs(m(j))~= inf)
                    d =  (m(i)+m(j) )/2;
                  else
                     % disp('infinity')
                  if (abs(m(i))==inf)
                        d = tan((pi/2 + atan(m(j)))/2);
                    else
                        d = tan((pi/2 + atan(m(i)))/2);
                  end
                  end
                  c = vy0(1) - d * vx0(1);
                  px = vx0(1) + 0.1;
                  py = d*px+c;
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                  px = vx0(1) - 0.1;
                  py = d*px+c;
                  
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                  
                  %disp("perpen")
                  
                  if d == 0
                      px = vx0(1); 
                      py = vy0(1)+ 0.1;
                  
                      [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                      px = vx0(1); 
                      py = vy0(1)- 0.1;
                  
                      [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                  
                  else
                      d = -1/d;
                  
                  c = vy0(1) - d * vx0(1);
                  px = vx0(1) + 0.1;
                  py = d*px+c;
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                  px = vx0(1) - 0.1;
                  py = d*px+c;
                  
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                  end
                end
              end 
              if size(m,2) == 1
                  d = -1/m(1)
                  c = vy0(1) - d * vx0(1);
                  px = vx0(1) + 0.1;
                  py = d*px+c;
                  
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                  px = vx0(1) - 0.1;
                  py = d*px+c;
                  
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  if l 
                    return
                  end
                  
        
              end
              %disp("SINGLETON REGION")
              %obj.print
              
              lsing = true;
              
              end
              

          if all(double(sv1) <= double(sv2))
              fmax = f2;
              index=2;
          elseif all(double(sv2) <= double(sv1))
              fmax = f1;
              index=1;
          else
              l = false;
              fmax = 0;
              index=0;
          end
         
        end
        
        function [l,fmax,index] = maxFromPt(obj, s, f)
          l = false;
          fmax = f(1);
          index=1;
          if obj.ptFeasible(obj.vars,s)
            fv01 = f(1).subsF(obj.vars,s).f;
            fv02 = f(2).subsF(obj.vars,s).f;
            
          else
              return;
          end
          l = true;
          if abs(double(fv01 - fv02))> 1.0d-14
            if double(fv01) < double(fv02)
              fmax = f(2);
              index=1;
            else
              fmax = f(1);
              index=2;
            end
            return
          end
          l = false;
          fmax = f(1);
          index=1;
        end

        function [r] = splitmax2 (obj, f1, f2) 
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          f = f1-f2;
          vars = f.getVars;
          x = sym('x');
          y = sym('y');
          f = subs(f.f, vars ,[x,y]);
          
          fx = [];
          fy = [];
          n = 0;
          fxy = [];
          for i = 1:size(obj.ineqs,2)
            g = subs(obj.ineqs(i).f, vars,[x,y]);
            s = solve ([f==0,g==0],[x,y]);
            if isempty(s)
              continue;
            end
            
            if isempty(s.x)
              continue;
            end
            sx = double(s.x);
            sy = double(s.y);
            if ~obj.ptFeasible (vars,[sx,sy])
                continue
            end
            fx = [fx,sx];
            fy = [fy,sy];
            n = n+1;
            fxy(n,1)= sx;
            fxy(n,2)= sy;
          end
          fxy
          fxy = unique(fxy,"rows");
          fx = fxy(:,1);
          fy = fxy(:,2);
          if size(fx,1) == 2
              m = (fy(2)-fy(1))/(fx(2)-fx(1));
              c = fy(1) - m*fx(1);
              ineq = vars(2)-m*vars(1)-c;
          end
          for i = 1:size(fv1,2)
              if (double(fv1(i).f) > double(fv2(i).f))
                  if subs(ineq,vars,[obj.vx(i),obj.vy(i)]) <= 0
                    r = [ineq,-ineq];
                    return
                  else
                    r = [-ineq,ineq];
                    return  
                  end
              end

          end
          
        end

        function [r] = splitmax3 (obj, f1, f2) 

          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          f = f1-f2;
           vars = f.getVars;
          ineq = f;
          for i = 1:size(fv1,2)
              if (abs(double(fv1(i).f)) >= abs(double(fv2(i).f)))
                  %ineq.f
                  %subs(ineq.f,vars,[obj.vx(i),obj.vy(i)])
                  if subs(ineq.f,vars,[obj.vx(i),obj.vy(i)]) < 0
                    r = [ineq,-ineq];
                    return
                  else
                    r = [-ineq,ineq];
                    return  
                  end
              end

          end
          
        end

        
        function [nl, v1, v2] = splitmaxArray (obj, f1, f2) 
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          v1 = [];
          v2 = [];
          for i = 1:size(fv1,2)
              sv1(i) = fv1(i).f
              sv2(i) = fv2(i).f
              % replace with limit
              if (isnan(sv1(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              elseif (isnan(sv2(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              end
              if sv1(i) == sv2(i)
                  v2 = [v2,i];
                  v1 = [v1,i];
              elseif sv1(i) > sv2(i)
                  v1 = [v1,i];
              else    
                  v2 = [v2,i];
              end

          end
          
          if size(v1,2) == 3
              nl(1) = 1;
          else
              nl(1) = 0;
          end
          if size(v2) == 3
              nl(2) = 1;
          else
              nl(2) = 0;
          end
          
        end
        
        function [lelim, obj] = deleteIneq (obj, vars)
          lelim = false;

          for i = 1:size(obj.ineqs,2)
              l(i) = obj.ineqs(i).f;
              mark(i) = false;
              nPts = 0;

              for j = 1:obj.nv
                  if obj.ineqs(i).subsF(obj.vars,[obj.vx(j),obj.vy(j)]).isZero
                      nPts = nPts+1;
                  end
              end
              if nPts <= 1
                  mark(i)=true;
              end
          end
          for i = 1:size(obj.ineqs,2)
              if ~mark(i)
                  continue
              end
              if obj.ineqs(i).isQuad
                  continue;
              end
             l1 = l;
             l1(i) = [];
             r = region (l1,vars);
             if isempty(r)
                 return
             end
             if obj.eqVertices(r)
                 lelim = true;
                 obj = r;
                 return;
             end
          end


        end

        function obj = simplifyOpenRegion1 (obj, nP, px, py)
            % remove ineqs that dont go through a vertex
            n = 0;
            mark = [];
            for i = 1:size(obj.ineqs,2)
                l = false;
                for j = 1:nP
                    if (abs(double(obj.ineqs(i).subsF(obj.vars,[px(j),py(j)]).f)) < 1.0d-8)
                        l = true;
                        break
                    end
                end
                if l
                    continue;
                end
                n = n + 1;
                mark(n) = i;
            end
            obj.ineqs(mark) = [];
           
        end

        function obj = removeTangent0 (obj, nP, px, py)
            n = 0;
            mark = [];
           
            for i = 1:size(obj.ineqs,2)
                vars = obj.ineqs(i).getVars();
                %obj.ineqs(i).subsF(vars,[px(i),py(i)])
              
              if ~ obj.ineqs(i).isQuad 
                continue
              end
              % get feasible point
              for j = 1:nP
                if ~obj.ineqs(i).subsF(vars,[px(j),py(j)]).isZero  
                  continue
                end
              sx = px(j)-0.1;
              sy = py(j);

              if ~ isAlways(obj.ineqs(i).subsF(vars,[sx,sy])>0)
                sx = px(j)+0.1;
                sy = py(j);
                if ~ isAlways(obj.ineqs(i).subsF(vars,[sx,sy])>0)
                  sx = px(j);
                  sy = py(j)+0.1;
                  if ~ isAlways(obj.ineqs(i).subsF(vars,[sx,sy])>0)
                    sx = px(j);
                    sy = py(j)-0.1;
                    if ~ isAlways(obj.ineqs(i).subsF(vars,[sx,sy])>0)
                        disp('feasible pt not found')
                        return
                    end
                  end
                end
              end
              % disp('sx')
              % double(sx), double(sy)
              
                tangent = obj.ineqs(i).tangent(px(j),py(j));
                tangent = tangent.normalize1;
                if ~ isAlways(tangent.subsF(vars,[sx,sy])> 0)
                    tangent = -tangent
                end
                disp('tangent')
                %px(j),py(j)
                % obj.ineqs(i).f
                % tangent.f
                % d1 = sqrt((sx-px(j))^2+(sy-py(j))^2)
                % lc = tangent.getLinearCoeffs (vars)
                % d2 = abs(double(tangent.subsF(vars,[sx,sy]).f/sqrt(lc(1)^2+lc(2)^2)))
                % pt = isAlways(d1 <= d2)
                %continue
                for k = 1:size(obj.ineqs,2)
                  l1 = isAlways(obj.ineqs(k).f==tangent.f);
                  %l2 = isAlways(obj.ineqs(k).f==-tangent.f);
                  % if l2
                  %     obj = region.empty();
                  %     return
                  % end
                  if ~l1 %& ~l2
                    continue
                  end

                  % fix pt
                      n = n + 1;
                      mark(n) = k;
                  break;
                  % check whether to remove tangent or parabola
                  % x = px(j) + 0.1;
                  % y = py(j);
                  % l = isAlways (obj.ineqs(i).subsF(x,y).f<= 0);
                  % if l1
                  %     %isAlways (tangent.subsF(x,y).f <= 0)
                  %   if (l & isAlways (tangent.subsF(x,y).f <= 0))
                  %     n = n + 1;
                  %     mark(n) = i;
                  %     break;
                  %   else
                  %     n = n + 1;
                  %     mark(n) = k;
                  %     break;
                  %   end
                  % end
                  % if l2
                  %     %isAlways (-tangent.subsF(x,y).f <= 0)
                  %   if (l & isAlways (-tangent.subsF(x,y).f <= 0))
                  %     n = n + 1;
                  %     mark(n) = i;
                  %     break;
                  %   else
                  %     n = n + 1;
                  %     mark(n) = k;
                  %     break;
                  %   end
                  % end
                end 
              end
            end
            obj.ineqs(mark) = [];
        end

        function obj = removeTangent (obj, nP, px, py)
            n = 0;
            mark = [];
            vars = obj.vars;
            for i = 1:size(obj.ineqs,2)
              if ~ obj.ineqs(i).isQuad 
                continue
              end
              % get feasible point
              for j = 1:nP
                if ~obj.ineqs(i).subsF(vars,[px(j),py(j)]).isZero  
                  continue
                end
                tangent = obj.ineqs(i).tangent(px(j),py(j));
                tangent = tangent.normalize1;
                %disp('tan')
                %tangent.print
                % Get point in parabola
                sx = px(j) + 0.1;
                p = obj.ineqs(i).subsF([vars(1)],[sx]);
                s0 = solve(p.f,vars(2));
                sy = s0(1);
                mx = (px(j) + sx)/2;
                my = (py(j) + sy)/2;
                lin = isAlways(obj.ineqs(i).subsF(vars,[mx,my]).f <= 0);
                for k = 1:size(obj.ineqs,2)
                  l1 = isAlways(obj.ineqs(k).f==tangent.f);
                  if isAlways(obj.ineqs(k).f==-tangent.f);
                      tangent = -tangent;
                      l1 = true;
                  end
                  %obj.ineqs(k).f
                  if ~l1
                      continue
                  end
                   [nvi, vxi, vyi] = obj.vertexOfEdge(k);
                  if nvi > 1
                      continue
                  end
                  tin = isAlways(tangent.subsF(vars,[mx,my]).f<= 0);
                  if lin & tin
                    n = n + 1;
                    mark(n) = k;
                  end  
                  if lin & ~tin
                    obj = region.empty;
                    return
                  end  
                  if ~lin & tin
                    % do nothing
                  end  
                  if ~lin & ~tin
                    n = n + 1;
                    mark(n) = i;
                  end  
                  %mark
                  break;
                end 
              end
            end
            %mark
            obj.ineqs(mark) = [];
            % check
            nv = obj.nv;
            obj = obj.getVertices;
            if nv ~= obj.nv
                disp('Error in removeTangent')
                disp('This is okay when 3 vertices are on parabola - its removing middle one')
            end
        end
        
        function obj = simplifyOpenRegion (obj)
            nP = 0;
            for j = 1:obj.nv
              if abs(obj.vx(j)) == intmax
                continue
              end
              if abs(obj.vy(j)) == intmax
                continue
              end
              nP = nP+1;
              px(nP) = obj.vx(j);
              py(nP) = obj.vy(j);
            end
            obj = obj.simplifyOpenRegion1 (nP, px, py);
            for j = 1:nP
                nPoint(j) = 0;
            end
            keep = [];
            for i = 1:size(obj.ineqs,2)
                keep(i) = false;
            end
            for i = 1:size(obj.ineqs,2)
                for j = 1:nP
                    
                    if (abs(obj.ineqs(i).subsF(obj.vars,[px(j),py(j)]).f) > 1.0d-8)
                        continue;
                    end
                    nPoint(j)=nPoint(j)+1;
                    point(j,nPoint(j)) = i;
                    if obj.ineqs(i).isQuad
                        continue
                    end
                    
                end
            end
            if all(nPoint == 2)
                return;
            end
            for i = 1:size(obj.ineqs,2)
               if obj.ineqs(i).isQuad
                   return
               end
            end
            m0 = obj.slopes;
            n0 = 0;
            n1 = 0;
            markF0 = [];
            for i = 1:size(obj.ineqs,2)
                markF0(i) = i;
               if obj.ineqs(i).isQuad
                n0 = n0+1;
                % replace by slope of tangent and keep it
                m(n0) = intmax;
               else
                n0 = n0+1;
                n1 = n1+1;
                m(n0) = m0(n1);
               end
            end
           
            
            markF = [];
            for ip = 1:nP
                sx = px(ip);
                sy = py(ip);
                pi0 = point(ip,1:nPoint(ip));
                
                mark = [];
                for j = 1:nPoint(ip)
                    if obj.ineqs(pi0(j)).isQuad
                        mark = [mark,j];
                    end
                end
                
                pi0(mark) = [];
                mp = m(pi0);
                [sorted_m, indices] = sort(mp);
              for i = 1:size(indices,2)
                m1 = sorted_m(i);
                if i == size(indices,2)
                    m2 = sorted_m(1);
                else
                    m2 = sorted_m(i+1);
                end
                if (abs(m1)~= inf) & (abs(m2)~= inf)
                  d =  (m1+m2 )/2;
                else
                    % put temp on 14 th May
                    return
                  if (abs(m1)==inf)
                    d = tan((pi/2 + atan(m2))/2);
                  else
                      
                    d = tan((pi/2 + atan(m1))/2);
                  end
                end
                if i == size(indices,2)
                    d = -d;
                end
                c = sy - d * sx;
                tx = sx + 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                   continue
                end
                tx = sx - 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                  continue
                end


                if i == 1
                    m2 = sorted_m(size(indices,2));
                else
                    m2 = sorted_m(i-1);
                end
                if (abs(m1)~= inf) & (abs(m2)~= inf)
                  d =  (m1+m2 )/2;
                else
                  if (abs(m1)==inf)
                    d = tan((pi/2 + atan(m2))/2);
                  else
                      
                    d = tan((pi/2 + atan(m1))/2);
                  end
                end
                if i == 1
                    d = -d;
                end
                c = sy - d * sx;
                tx = sx + 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                   continue
                end
                tx = sx - 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                  continue
                end
                markF = [markF,pi0(indices(i))];
                
              end
              for i = 1:size(obj.ineqs,2)
                  keep0(i) = false;
              end
              for i = 1:size(pi0,2)
                  keep0(pi0(i)) = true;
              end
              for i = 1:size(markF,2)
                  keep0(markF(i)) = false;
              end
              keep = keep | keep0; 
            end
            markF0 = [];
            for i = 1:size(pi0,2)
                if keep(pi0(i))
                    continue
                end
                markF0 = [markF0,pi0(i)];
            end
            
            if size(markF0,2) == size(obj.ineqs,2)
                obj = region.empty;
                disp("Singleton")
                return
                
            end
            obj.ineqs(markF0) = [];
            % check vertices


        end
            
        function obj = simplifyUnboundedRegion (obj)
            
            if isempty(obj)
                return
            end
            % Removing inf from vertices - to be removed later
            nP = 0;
            for j = 1:obj.nv
              if abs(obj.vx(j)) == intmax
                continue
              end
              if abs(obj.vy(j)) == intmax
                continue
              end
              nP = nP+1;
              px(nP) = obj.vx(j);
              py(nP) = obj.vy(j);
            end

            for i = 1:size(obj.ineqs,2)
              for j = i+1:size(obj.ineqs,2)
                  if isAlways(obj.ineqs(i) == -obj.ineqs(j))
                      obj = region.empty();
                      return
                  end
              end
            
            end
            obj = obj.simplifyOpenRegion1 (nP, px, py);
            % disp('b4 remove tangent')
            % obj.print
            %obj = obj.removeTangent (nP, px, py);
            if isempty(obj)
                return
            end
            % disp('aft remove tangent')
            % obj.print
            
            for j = 1:nP
                nPoint(j) = 0;
            end
            keep = [];
            for i = 1:size(obj.ineqs,2)
                keep(i) = false;
            end

            % get Ineqs going through a vertex
            for i = 1:size(obj.ineqs,2)
                for j = 1:nP
                    if symFunType(obj.ineqs(i).subsF(obj.vars,[px(j),py(j)]).f) == 'plus'
                        continue;
                    end
                    if (abs(obj.ineqs(i).subsF(obj.vars,[px(j),py(j)]).f) > 1.0d-8)
                        continue;
                    end
                    nPoint(j)=nPoint(j)+1;
                    point(j,nPoint(j)) = i;
                    
                end
            end
            %disp('point')
            %nPoint
            %point
            if all(nPoint == 0)
                obj.print
                nP
                px
                py
            end
            
             
            
            % make a function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nP == 1 & nPoint(1)==2
                mp = obj.slopeAtVertex(point(1,1:2), [px,py]);
                theta = atan2(mp,1);
                t1 = theta(1);
                t2 = theta(2);
                ind = [1,2];
                d = tan((t1+t2)/2);
                
                if isnan(d) | isinf(d)
                    tx = px(1);
                    ty = py(1)+0.1;
                else
                    c = py(1) - d * px(1);
                    
                    tx = px(1) + 0.1;
                    ty = d*tx+c;
                end
                
                lm = obj.ptFeasible (obj.vars,[tx,ty]);
                if isnan(d)| isinf(d)
                    tx = px(1);
                    ty = py(1)-0.1;
                else
                
                    tx = px(1) - 0.1;
                    ty = d*tx+c   ;
                end
                
                lm = lm | obj.ptFeasible (obj.vars,[tx,ty]);
                d = -1/d;
                
                if isnan(d) | isinf(d)
                    tx = px(1);
                    ty = py(1)+0.1;
                else
                    c = py(1) - d * px(1);
                    
                    tx = px(1) + 0.1;
                    ty = d*tx+c;
                end
                
                lm = lm | obj.ptFeasible (obj.vars,[tx,ty]);
                if isnan(d)| isinf(d)
                    tx = px(1);
                    ty = py(1)-0.1;
                else
                
                    tx = px(1) - 0.1;
                    ty = d*tx+c   ;
                end
                
                lm = lm | obj.ptFeasible (obj.vars,[tx,ty]);
                 


                
               if ~lm 
                   obj = region.empty()
              
               end
               
               return
              
            end   


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if all(nPoint == 2)
                return;
            end
            
            %point
            for ip = 1:nP
                sx = px(ip);
                sy = py(ip);
                pi0 = point(ip,1:nPoint(ip));
                markF = [];
                if size(pi0,2) == 2
                    keep(pi0) = true;
                    continue
                end
                if nPoint(ip) == 0
                    continue
                end
                

                mp = obj.slopeAtVertex(pi0, [sx,sy]);
                
                for i = 1:size(mp,2)
                    if mp(i) == -inf
                        mp(i) = inf;
                    end
                end

                

                [sorted_mp, indices] = sort(mp);
                theta = atan2(mp,1);
                sorted_theta = atan2(sorted_mp,1);
                % theta
                % [temp_theta, temp_indices] = sort(theta)
                %return
              for i = 1:size(indices,2)
                 
                t1 = sorted_theta(i);
                 
                
                if i == size(indices,2)
                    t2 = sorted_theta(1);
                    ind = [pi0(indices(i)),pi0(indices(1))];
                    l90 = true;
                else
                    t2 = sorted_theta(i+1);
                    ind = [pi0(indices(i)),pi0(indices(i+1))];
                    l90 = false;
                end
                if t1 == t2
                    return
                    %break
                end
                
                if l90
                    d = tan(pi/2+ (t1+t2)/2);
                else
                    d = tan((t1+t2)/2);
                end
                %d = tan((t1+t2)/2);
%                m1, m2, double(d)
                if isnan(d) | isinf(d)
                    tx = sx;
                    ty = sy+0.1;
                else
                    c = sy - d * sx;
                    tx = sx + 0.1;
                    ty = d*tx+c;
                end
                %double(tx),double(ty) 
                %l = obj.ptFeasible (obj.vars,[tx,ty]);

                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                  % disp('c1')     
                   continue
                end
                if isnan(d)| isinf(d)
                    tx = sx;
                    ty = sy-0.1;
                else
                
                    tx = sx - 0.1;
                    ty = d*tx+c   ;
                end
                %double(tx),double(ty) 
                %l = obj.ptFeasible (obj.vars,[tx,ty])
                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                 % disp('c2')     
                  continue
                end
                


                
                if i == 1
                    %m2 = sorted_m(size(indices,2));
                    t2 = sorted_theta(size(indices,2));
                    ind = [pi0(indices(i)),pi0(indices(size(indices,2)))];
                    l90 = true;
                else
                    %m2 = sorted_m(i-1);
                    t2 = sorted_theta(i-1);
                    ind = [pi0(indices(i)),pi0(indices(i-1))];
                    l90 = false;
                end
                if t1 == t2
                    return
                    %break
                end
                
                
                if l90
                    d = tan(pi/2+ (t1+t2)/2);
                else
                    d = tan((t1+t2)/2);
                end
 %               m1, m2, double(d)
         %       atan(m1)
         %       m1,m2,d
                % if i == 1
                %     d = -d;
                % end
               % d
                %isnan(d)
                if isnan(d)| isinf(d)
                    tx = sx;
                    ty = sy+0.1;
                else
                    c = sy - d * sx;
                    tx = sx + 0.1;
                    ty = d*tx+c;
                end
                
                %d,tx,ty
                %double(tx),double(ty) 
                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                   %disp('c3')     
                   continue
                end
                if isnan(d)| isinf(d)
                    tx = sx;
                    ty = sy-0.1;
                else
                
                    tx = sx - 0.1;
                    ty = d*tx+c   ;    
                end
                %double(tx),double(ty) 
                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                 % disp('c4')     
                  continue
                end
                %disp('end');
                %disp('marked')
                %i
             %   indices(i)
                %pi0(indices(i))
                markF = [markF,pi0(indices(i))];
                
              end
              %markF
              for i = 1:size(obj.ineqs,2)
                  keep0(i) = false;
              end
              for i = 1:size(pi0,2)
                  keep0(pi0(i)) = true;
              end
              for i = 1:size(markF,2)
                  keep0(markF(i)) = false;
              end
              %keep0
              %keep
              keep = keep | keep0 ;
              %keep
              
            end
            %keep
            %markF
            %if ~isempty (markF)
                if size(markF,2) == size(obj.ineqs,2)
                disp("Singleton1")
                %obj.print
                
                    obj = region.empty;
                    disp("Singleton")
                    return
                    
                end
            %end
            %obj.ineqs(markF) = [];
            %obj.ineqs = obj.ineqs(markF);
            % check vertices
            markF0 = [];
            for i = 1:size(obj.ineqs,2)
                if keep(i)
                    continue
                end
                markF0 = [markF0,i];
            end
            %markF0
            if size(markF0,2) == size(obj.ineqs,2)
                disp("Singleton0")
                obj.print
                obj = region.empty;
                disp("Singleton")
                return
                
            end
            objTemp = obj;
            obj.ineqs(markF0) = [];
            nv = obj.nv;
            obj = obj.getVertices;
            %obj = obj.removeTangent (nP, px, py);
            if nv ~= obj.nv
                disp ("Error in simplify");
                %sorted_m
                objTemp.print
                obj.print
            end


        end
        


        function obj = simplify (obj) 
          for i = 1:size(obj.ineqs,2)
            [lelim, obj] = obj.deleteIneq (obj.vars);
            if ~lelim
                return
            end
          end

        end


        %function simplifyClosedRegion
        %end
         function l = ptFeasible(obj, vars, point)
           l = true;
           
           % disp('point')
           % point
           % double(point)
           % size(point,1)
           
           for i = 1:size(obj.ineqs,2)
               %obj.ineqs(i).f
               for j = 1:size(point,1)
                %  disp('subs')
                 %  double(point(j,:))
                  
                  %i
                   % isAlways(subs ([obj.ineqs(i).f],vars,point(j,:))>0)
                   % double(subs ([obj.ineqs(i).f],vars,point(j,:)))
                %if    (subs ([obj.ineqs(i).f],vars,point(j,:)) > 0)
                % double used here as sqrt giving issues
               % if double(subs ([obj.ineqs(i).f],vars,double(point(j,:)))) > 1.0d-12
               % point(j,:)
               if (isAlways(subs ([obj.ineqs(i).f],vars,point(j,:))>0))
                   l = false;
                   return
               end
               end
           end  
         end

         function l = ptFeasibleSubset(obj, vars, point, index)
           l = true;
          % point
           
           
           for i = 1:size(index,2)
               %obj.ineqs(index(i)).f    
               for j = 1:size(point,1)
               %    subs ([obj.ineqs(index(i)).f],vars,point(j,:))
                if    (subs ([obj.ineqs(index(i)).f],vars,point(j,:)) > 0)
                %if double(subs ([obj.ineqs(i).f],vars,double(point(j,:)))) > 1.0e-12
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
               for j = i+1:size(obj.ineqs,2)
                 if ( obj.ineqs(i) == unaryminus(obj.ineqs(j)))
                     return
                 end
                 if ( obj.ineqs(i) == obj.ineqs(j))
                     continue
                 end
                 
                 s = solve(obj.ineqs(i).f<=0,obj.ineqs(j).f<=0);
                 if isempty(s)
                     return;
                 end
               end 
 
             end 
             
             l = true;
         end

         % stupid way of doing this
         function obj = intersection (obj1, obj2)
         l = [];
           for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
           end 
         for i = 1:size(obj2.ineqs,2)
             lf = true;
             for j = 1:size(obj1.ineqs,2)
               if obj1.ineqs(j)==obj2.ineqs(i)
                   lf = false;
                   break;
               end
             end
             if lf
               l = [l,obj2.ineqs(i).f];
             end
         end 
         obj = region(l, obj1.vars);
         if isempty(obj)
             return
         end
         %obj.print
         if (isFeasible(obj))
             %disp('feasible')
             obj = obj.unique;
             if obj.nv <= 2
                 %disp('degenerate ');
                 obj = region.empty;
             else
                 if obj.nv == 0
           disp('in intersection nv 0')
         end

             end
         else
             obj = region.empty;
         end

     end
     
     function obj = normalizeEdge(obj)
         for i = 1:size(obj.ineqs,2)  
             c = getLinearCoeffs (obj.ineqs(i),obj.vars);
             
             if (c(2) == 0)
                 % double * f not overloaded
                 obj.ineqs(i).f = (1/c(1)) * obj.ineqs(i).f;
             else
                 obj.ineqs(i).f = (1/c(2)) * obj.ineqs(i).f;
             end
         end
     end

     function obj = normalize1(obj)
         for i = 1:size(obj.ineqs,2)
             obj.ineqs(i) = obj.ineqs(i).normalize1;
         end
     end
     
     function V = getIntersectingFeasiblePts(obj, f0)
         
       nv=0;
       V = sym.empty();
       %obj.vx=[];
       %obj.vy=[];
       t1 = sym('t1');
       t2 = sym('t2');
       varsTemp = [t1,t2];
       f2 = f0.subsF (obj.vars,varsTemp);
       for i = 1:size(obj.ineqs,2)  
           f1 = obj.ineqs(i);
           f1 = f1.subsF (obj.vars,varsTemp);
               
               s = solve ([f1.f==0,f2.f==0],varsTemp);
               
               if isempty(s)
                   continue;
               elseif isempty(s.t1)
                   continue;
               elseif isempty(s.t2)
                   continue;
               end
               
               
               for k = 1:size(s.t1,1)
                    if ~ isreal(double(s.t1(k)))  % weird error - not detecting complex in symbolic z^4
                        continue
                    end
                    if ~ isreal(double(s.t2(k)))
                        continue
                    end
                    
                    if obj.ptFeasible(obj.vars, [s.t1(k),s.t2(k)])
                       
                          nv=nv+1;
                          V(nv,1) = s.t1(k);
                          V(nv,2) = s.t2(k);
                        
                    end
               end 
           
       end
     end
     
     % wont work for degree > 2
     % removed complex vertices
     
     function obj = linear3pt(obj)
       mark = [];  
      
       for i = 1:size(obj.ineqs,2) 
         if obj.ineqs(i).isQuad  
             continue
         end
         [nvi, vxi, vyi] = obj.vertexOfEdge(i);
         if nvi <= 2
             continue
         end
         % assuming nvi = 3
         [sortvx, ind] = sort(double(vxi));  % rational sorting goofed up
         x = (vxi(ind(1)) + vxi(ind(2)))/2;
         y = (vyi(ind(1)) + vyi(ind(2)))/2;
         %obj.ptFeasible (x,y) 
         if ~obj.ptFeasible (obj.vars,[x,y])
             mark(1) = ind(1);
         else
             mark(1) = ind(3);
         end
         break;

       end  
       
       if isempty(mark)
           return
       end
       %obj.print
       %mark
       for i = 1:obj.nv
           if (obj.vx(i) ==  vxi(mark(1))) & (obj.vy(i) ==  vyi(mark(1)))
               mark(1) = i;
           end
       end
       obj.vx(mark(1)) = [];
       obj.vy(mark(1)) = [];
       obj.nv = obj.nv-1;
     end 

     function obj = removeInfV (obj)
         n = 0;
         for i = 1:obj.nv
             
             if isAlways (abs(obj.vx(i)) == intmax)
                 continue
             end
             if isAlways (abs(obj.vy(i)) == intmax)
                 continue
             end
             n = n + 1;
             obj.vx(n) = obj.vx(i);
             obj.vy(n) = obj.vy(i);
             
         end
         obj.nv = n;
     end 

     function obj = getVertices(obj)
         
       obj.nv=0;
       vx=sym.empty();
       vy=sym.empty();
       %obj.vx=[];
       %obj.vy=[];
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
               
               
               for k = 1:size(s.t1,1)
                    if ~ isreal(double(s.t1(k)))  % weird error - not detecting complex in symbolic z^4
                        continue
                    end
                    if ~ isreal(double(s.t2(k)))
                        continue
                    end
                    
                    if obj.ptFeasible(obj.vars, [s.t1(k),s.t2(k)])
                       %i,j,obj.nv
                       %isAlways(subs ([obj.ineqs(i).f],obj.vars,[s.t1(k),s.t2(k)])<=0)
                       %double(subs ([obj.ineqs(i).f],obj.vars,[s.t1(k),s.t2(k)]))
                          obj.nv=obj.nv+1;
                          obj.vx(obj.nv) = s.t1(k);
                          obj.vy(obj.nv) = s.t2(k);
                        
                    end
               end 
           end
               
       
       end
       if obj.nv ~= 0
           %disp('b4 sort')
           %obj.nv
           %obj.vx
           %double(obj.vx)
       for i = 1:obj.nv
           V(i,1) = obj.vx(i);
           V(i,2) = obj.vy(i);
       end 
       
       V = unique(V,"rows");
       if obj.nv ~= size(V,1)
         obj.nv = size(V,1);
         obj.vx = V(:,1)';
         obj.vy = V(:,2)';
       end
           obj = obj.linear3pt; 
       end  
       %[obj.vx,obj.vy] = poly2cw(obj.vx,obj.vy);

       % putting intmax for inf to avoid Nans 
       % intmax + intmax = intmax
       %disp('ptFeasile')
       %obj.print
       n = obj.nv;
       % if n== 0
       %     return
       % end
       %return
       % to remove inf - fix + in convex envelope solve
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
               %intmax,obj.vy(i)
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
       
     end

     function [edgeNo] = getEdgeNosInf(obj, vars)
         edgeNo = zeros(size(obj.ineqs,2),1);
         
         for i = 1:size(obj.ineqs,2)
            
           [nv, vx, vy] = obj.vertexOfEdge(i);
           if nv == 1
               if obj.vx(1)==vx(1) & obj.vy(1)==vy(1)
                   edgeNo(i) = 1;
               else
                   edgeNo(i) = obj.nv+1;
               end
               continue
           end
            for j = 1:obj.nv
                  if obj.vx(j)==vx(1) & obj.vy(j)==vy(1)
                      break;
                  end    
            end
            edgeNo(i) = j+1;
           
         end
         
     end

     function [edgeNo] = getEdgeNos(obj, vars)
       edgeNo = zeros(obj.nv,1);
       for j = 1: obj.nv-1
         slope = obj.slope(j,j+1);
         if slope == -inf
           slope = inf;
         end
         if slope == inf
           edge = vars(1) -obj.vx(j) ; 
         else
           q = obj.yIntercept (j,slope);
           edge = vars(2)-slope*vars(1)-q;
         end
         for k = 1: size(obj.ineqs,2)
           e0 = obj.ineqs(k);
           e0 = e0.normalize (vars);
           if (e0.f == edge)
             break;
           end
         end
         edgeNo(j)=k;
       end
       j = obj.nv;
       slope = obj.slope(j,1);
       if slope == inf
         edge = vars(1) -obj.vx(j)  ;
       else
         q = obj.yIntercept (j,slope);
         edge = vars(2)-slope*vars(1)-q;
       end
       for k = 1: size(obj.ineqs,2)
         e0 = obj.ineqs(k);
         e0 = e0.normalize (vars);
         if (e0.f == edge)
           break;
         end
       end
       edgeNo(j)=k;
     end

      function [nv, vx, vy] = vertexOfEdge(obj,ind)
              nv = 0;
              vx = sym.empty();
              vy = sym.empty();
              %disp('vertex of edge')
              %obj.nv
              %obj.ineqs(ind)
               for i = 1:obj.nv
                 f1 = obj.ineqs(ind).subsF (obj.vars,[obj.vx(i),obj.vy(i)]);
                 if f1.isZero
                     nv = nv + 1;
                     vx(nv) = obj.vx(i);
                     vy(nv) = obj.vy(i);
                 end
               end
      end

      function edges = getEdges(obj,vx,vy)
          n = 0;
          for i = 1:size(obj.ineqs,2)
              if isZero(obj.ineqs(i).subsF(obj.vars,[vx,vy]))
                  n = n + 1;
                  edges(n) = i;
              end
          end
          
      end

      function f = divideRegions(obj,obj2)

       obj.nv=0;
         t1 = sym('t1');
         t2 = sym('t2');
         varsTemp = [t1,t2];
         ir = 0;
         % in this loop we are getting all points of intersection of pair
         % of ineqs.
         % remove duplicate points 
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
               lpt = false;
               for k  = 1:ir
                   if all(points(k,:) == [s.t1,s.t2])
                       lpt = true;
                       break;
                   end
               end
               if lpt
                   continue
               end
               if ~obj2.ptFeasible (obj2.vars,[s.t1,s.t2])
                 continue
               end

               ir = ir + 1;
               index(ir) = 0;
               points(ir,1:2) = [s.t1,s.t2]; 
               for k = 1:size(obj.ineqs,2)  
                 if double(subs ([obj.ineqs(k).f],obj.vars,[s.t1,s.t2])) < 1.0e-12
                     index(ir) = index(ir)+1;
                     r(ir,index(ir)) = obj.ineqs(k).f;
                 end
               end

           end

         end
 
          is = 0;
         s = r;
         for i=1:ir
                 leq = false;
                 for j = 1:is
                     if indexS(j) ~= index(i)
                         continue;
                     end
                     leq = true;
                     for k = 1:index(i)
                         if ~ isequal(r(i,k),s(j,k))
                             leq = false;
                             break;
                         end
                     end
                     if leq
                       break;
                     end
                 end
                 if leq
                    continue;
                 end
                 is = is+1;

                 indexS(is) = index(i);
                 for j = 1:indexS(is)
                     s(is,j) = r(i,j);
                 end

         end

         f = [];
         for i = 1:is
           f = [f, region(s(i,:),obj.vars)]; 
         end

       end

  
 
     function fv = funcVertices (obj, f)
         n = 0;
         for i =1:obj.nv
             if (abs(obj.vx(i))==intmax)
                 continue
             end
             if (abs(obj.vy(i))==intmax)
                 continue
             end
             n = n + 1;
             fv(n) = f.subsF(obj.vars,[obj.vx(i),obj.vy(i)]);
             if (isnan(fv(n).f))
                 fv(n) = f.limit(obj.vars,[obj.vx(i),obj.vy(i)]);
             end
         end
     end


     function [l, maxf, index, lSing] = maximum(obj, f)
          
          if f(1)==f(2)
              l = true;
              lSing = false;
              maxf = f(1);
              index = 1;
              return
          end
          [l,  maxf, index, lSing] = obj.maxArray (f(1), f(2)) ;
     end

     function [f2, fe, d] = splitMax(obj, f, expr)
       [nl, v1, v2] = obj.splitmaxArray (f(1), f(2)) 
       f2 = [];
       fe = [];
       d = [];
       v{1}=v1;
       v{2}=v2;
       for i = 1:2
           if nl(i) == 0
               continue
           end
           f2 = [f2,f(i)];
           fe = [fe,expr(i)];
           d0 = obj.regionWPts(obj.vx(v{i}), obj.vy(v{i}), obj.vars(1), obj.vars(2));
           d = [d,d0];
           % disp('Regiopn')
           % d0.print
       end
       
     end
     
     function e = getOtherEdgeAtVertex (obj, ind, vertex)
         e = 0;
         %disp('vertex')
         %vertex
         for i = 1:size(obj.ineqs,2)
            if i == ind
                continue
            end
            f1 = obj.ineqs(i).subsF (obj.vars,vertex);
            %double(obj.ineqs(i).subsF (obj.vars,vertex).f)
            if f1.isZero
                
                e = i;
                break
            end

         end
     end


     function l = isconvex (obj, obj2, i, j, vx, vy)
       l = false;
      
       edgeiNo = obj.getOtherEdgeAtVertex (i,[vx,vy]);
       if edgeiNo == 0
          return
       end
       % if ~obj.ineqs(edgeiNo).isLinear 
       %   return
       % end
       edgejNo = obj2.getOtherEdgeAtVertex (j,[vx,vy]);
       if edgejNo == 0
         return
       end
       % if ~obj2.ineqs(edgejNo).isLinear 
       %   return
       % end
       
       if isAlways (obj.ineqs(edgeiNo) == obj2.ineqs(edgejNo))
           l = true;
           return
       end

       if obj.slopeIneq(edgeiNo,[vx,vy]) == inf
           x1 = vx;
           y1 = vy + 0.1;
           if ~ obj.ptFeasible (obj.vars,[x1,y1])
             y1 = vy - 0.1;
         
           end
       
       else

       x1 = vx+0.1;
       f1 = obj.ineqs(edgeiNo).subsF([obj.vars(1)],[x1]);
       y0 = solve(f1, obj.vars(2));
       %if (size(y0,1) == 1)
           y1 = y0(1);
       %end
       
       if ~ obj.ptFeasible (obj.vars,[x1,y1])
         x1 = vx-0.1;
         f1 = obj.ineqs(edgeiNo).subsF([obj.vars(1)],[x1]);
         y0 = solve(f1, obj.vars(2));
       %if (size(y0,1) == 1)
           y1 = y0(1);
       %end
       
       end
       end
       %x1,y1
       if obj2.slopeIneq(edgejNo,[vx,vy]) == inf
           x2 = vx;
           y2 = vy + 0.1;
           if ~ obj2.ptFeasible (obj2.vars,[x2,y2])
             y2 = vy - 0.1;
         
           end
       
       else

       
       x2 = vx+0.1;
       f1 = obj2.ineqs(edgejNo).subsF([obj2.vars(1)],[x2]);
       y0 = solve(f1, obj2.vars(2));
       %if (size(y0,1) == 1)
           y2 = y0(1);
       %end
       if ~ obj2.ptFeasible (obj2.vars,[x2,y2])
         x2 = vx-0.1;
         f1 = obj2.ineqs(edgejNo).subsF([obj2.vars(1)],[x2]);
         y0 = solve(f1, obj2.vars(2));
         %if (size(y0,1) == 1)
           y2 = y0(1);
         %end
       
       end
       end
       %x2,y2
       xm = (x1 + x2)/2;
       ym = (y1 + y2)/2;
       if obj.ptFeasible (obj.vars,[xm,ym]) | obj2.ptFeasible (obj2.vars,[xm,ym])
          l = true;
       end
       %disp('is convex')
       %l
      end


%      % wont work cause of intersection vs union
     function [l,obj] = merge (obj, obj2)
         l = false;
         n = 0; 
%          lQuad = false;
%          for i =1:size(obj.ineqs,2)
%               if obj.ineqs(i).isQuad
%                   lQuad = true;
%                   %return
%               end
%          end
%          for i =1:size(obj2.ineqs,2)
%               if obj2.ineqs(i).isQuad
%                   lQuad = true;
%                   %return
%               end
%          end

         % to be tested and added
         lQuad1 = false;
         nmq1 = 0;
         %size(obj.ineqs,2)
         for i =1:size(obj.ineqs,2)
             if obj.ineqs(i).isQuad
                 lQuad1 = true;
                 nmq1 = nmq1 + 1;
                 mq1(nmq1) = i;
             end
         end
         lQuad2 = false;
         nmq2 = 0;

         for i =1:size(obj2.ineqs,2)
             if obj2.ineqs(i).isQuad
                 lQuad2 = true;
                 nmq2 = nmq2 + 1;
                 mq2(nmq2) = i;

             end
         end
         %lQuad1
         %lQuad2
         if (lQuad1 & lQuad2)
             %disp("Quad merge")
             %obj.print
             %obj2.print
             marki = [];
             markj = [];
             for i = 1:nmq1
               for j = 1:nmq2
                   if obj.ineqs(mq1(i)) == -obj2.ineqs(mq2(j))
                       n = n + 1;
                       marki(n) = mq1(i);
                       markj(n) = mq2(j);
                   end
               end
             end
             if n > 0
               l = true;
               obj3 = obj;
               obj.ineqs(marki) = []; 
               obj2.ineqs(markj) = [];
               obj = obj+obj2;
               %disp("returning from quad")
               if isempty(obj)
                   %disp('reverting')
                   l = false;
                   obj=obj3;
               end
               return
             end
%          elseif lQuad1
%              return
%          elseif lQuad2
%              return
          end
         lQuad = lQuad1 | lQuad2;
         if lQuad1 & lQuad2
             for i = 1:nmq1
               for j = 1:nmq2
                   if obj.ineqs(mq1(i)).f ~= obj2.ineqs(mq2(j)).f
                       return;
                   end
               end
             end
             
             lQuad= false;
         end
         if lQuad
    %         lQuad1, lQuad2
             if lQuad1
                 
                 for i = 1:nmq1
                     V = obj2.getIntersectingFeasiblePts(obj.ineqs(mq1(i)));

                     if ~ isempty(V)
                         return
                     end
                 end
             end
            % V
             
         end
         % cant do + as it returns empty due to contadictory ineqs
         %obj = obj + obj2;
         % hence remove those then do +
         %obj.print;
       %  disp('in merge')
       %  obj.print
       %  obj2.print
         %obj = obj.unique;
         %disp("unique")
         %obj.print;

         %disp('after q')
         marki = [];
         markj = [];
         for i =1:size(obj.ineqs,2)

           for j =1:size(obj2.ineqs,2)
               %j
             if obj.ineqs(i) == -obj2.ineqs(j)
                 [nvi, vxi, vyi] = obj.vertexOfEdge(i); % move to outer loop
                 [nvj, vxj, vyj] = obj2.vertexOfEdge(j);
                 % obj.ineqs(i)
                 % nvi,nvj
                 % nvi 3 : crazy case - line with 3 vertices - example in
                 % testmerge5
                 if nvi ~= nvj
                     continue
                 end
%                  if nvi ~= 2
%                      continue
%                  end
                 %nvi
                 if nvi == 1 & (~lQuad)
                     
                     if obj.isconvex (obj2, i, j, vxi(1), vyi(1))

                     l = true;
                     n = n + 1;
                     marki(n) = i;
                     markj(n) = j;
                     end
                     break;

                 end
                 if nvi == 2

                 [tvxi, sortedIndices] = sort(abs(vxi));
                 vxi = vxi(sortedIndices);
                 vyi = vyi(sortedIndices);

                 [tvxj, sortedIndices] = sort(abs(vxj));
                 vxj = vxj(sortedIndices);
                 vyj = vyj(sortedIndices);

                 % check same vertices and convex angles 
                 % implement this part
                 % [l, pos] = obj.isVertexIrrational 
                 % [l, pos] = obj2.isVertexIrrational 
                 if all(vxi == vxj) & all(vyi == vyj)
                  %   disp('here')
                  %   obj.isconvex (obj2, i, j, vxi(1), vyi(1))
                     if ~obj.isconvex (obj2, i, j, vxi(1), vyi(1))
                       continue
                     end
                     
                     % edgeiNo = obj.getOtherEdgeAtVertex (i,[vxi(1),vyi(1)]);
                     % if edgeiNo == 0
                     %     continue
                     % end
                     % if ~obj.ineqs(edgeiNo).isLinear 
                     %     continue;
                     % end
                     % edgejNo = obj2.getOtherEdgeAtVertex (j,[vxi(1),vyi(1)]);
                     % if edgejNo == 0
                     %     continue
                     % end
                     % 
                     % if ~obj2.ineqs(edgejNo).isLinear 
                     %     continue;
                     % end
                     % 
                     % mi = obj.slopeIneq(edgeiNo);
                     % a1 = atan(mi);
                     % %if a1 < 0
                     % %    a1 = a1 + 2*pi;
                     % %end
                     % mj = obj2.slopeIneq(edgejNo);
                     % a2 = atan(mj);
                     % %if a2 < 0
                     % %    a2 = a2 + 2*pi;
                     % %end
                     % 
                     % if (a1+a2 > pi) | (a1+a2 < -pi)
                     %     continue
                     % end

%%%%%%%%%%%%%%%%

                    if nvi == 2   
                   %     disp('here2')
                   %     obj.isconvex (obj2, i, j, vxi(2), vyi(2)) 
                        if ~obj.isconvex (obj2, i, j, vxi(2), vyi(2))
                            continue
                        end

                     % edgeiNo = obj.getOtherEdgeAtVertex (i,[vxi(2),vyi(2)]);
                     % if edgeiNo == 0
                     %     continue
                     % end
                     % 
                     % if ~obj.ineqs(edgeiNo).isLinear 
                     %     continue;
                     % end
                     % edgejNo = obj2.getOtherEdgeAtVertex (j,[vxi(2),vyi(2)]);
                     % if edgejNo == 0
                     %     continue
                     % end
                     % 
                     % if ~obj2.ineqs(edgejNo).isLinear 
                     %     continue;
                     % end
                     % 
                     % mi = obj.slopeIneq(edgeiNo);
                     % a1 = atan(mi);
                     % %if a1 < 0
                     % %    a1 = a1 + 2*pi;
                     % %end
                     % mj = obj2.slopeIneq(edgejNo);
                     % a2 = atan(mj);
                     % %if a2 < 0
                     % %    a2 = a2 + 2*pi;
                     % %end
                     % 
                     % if (a1+a2 > pi) | (a1+a2 < -pi)
                     %     continue
                     % end
                    end

%%%%%%%%%%%%%%%%%

                         %atan(mi)+atan(mj)
                       
                         l = true;
                       n = n + 1;
                       marki(n) = i;
                       markj(n) = j;
                       break;
                     end
                 end  
             end
           end
           %if l
           %    break
           %end
         end
         if l
           % obj.print;
           % obj2.print;
           % marki
           % markj
           obj3 = obj;
           obj.ineqs(marki) = []; 
           obj2.ineqs(markj) = [];
           obj = obj+obj2;
           %obj.print
           obj = obj.simplifyUnboundedRegion;
           if isempty(obj)
               disp('empty')
               l = false;
               obj = obj3;
               % obj.print
           end
           %obj.print
         end

         end



     function l = hasNegativeIneqs(obj)
         l = true;
         for i = 1:size(obj.ineqs,2)
           for j = i+1:size(obj.ineqs,2)
               if (obj.ineqs(i) == -obj.ineqs(j))
                   return;
               end
           end
         end
         l = false;
     end

     end

     methods % limits
       function [l,limf] = limitOfFAtVertices (obj, f)
         vars = obj.vars;  
         vars2 = [vars(2),vars(1)];
         for j = 1: obj.nv
           l1 = f.limit(vars,[obj.vx(j),obj.vy(j)]);
           l2 = f.limit(vars2,[obj.vy(j),obj.vx(j)]);
           if (l1 == l2)
             l(j) = true;
             %limf(j) = double(l1.f);
             limf(j) = l1.f;
           else
             ldrx1(j) = false;
             limdrx1(j) = 0;
           end
         end
       end 
     end

     methods % normal cone
         function NC = getNormalConeVertex(obj, s1, s2)
             % obj = obj.envelope(i).d
            NC = sym(zeros(obj.nv,2));
            meanx = sum(obj.vx)/obj.nv;
            meany = sum(obj.vy)/obj.nv;
            
             for j = 1: obj.nv-1
                slope = obj.slope(j,j+1);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j);
                end
                
                if subs(eq,[s1,s2],[obj.vx(j+1),obj.vy(j+1)]) < 0
                    eq = -eq;
                end
                NC(j,1) = eq;
                if pslope ~= inf
                    q = obj.yIntercept (j+1,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j+1);
                end
                if subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) < 0
                    eq = -eq;
                end
                NC(j+1,2) = eq;
             end
             j = obj.nv;
             slope = obj.slope(j,1);
             pslope = -1/slope;
             if pslope == -inf
                    pslope = inf;
                end
                
             if pslope ~= inf
               q = obj.yIntercept (j,pslope);
               eq = s2 - pslope*s1 - q;
             else
               eq = s1 - obj.vx(j);
             end
             if subs(eq,[s1,s2],[obj.vx(1),obj.vy(1)]) < 0
               eq = -eq;
             end
              
             NC(j,1) = eq;
             if pslope ~= inf
               q = obj.yIntercept (1,pslope);
               eq = s2 - pslope*s1 - q;
             else
               eq = s1 - obj.vx(1);
             end
             if subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) < 0
               eq = -eq;
             end
                
             NC(1,2) = eq;
                
        end

        function [NC] = getNormalConeEdge(obj, s1, s2)
          NC = sym(zeros(obj.nv,2));
          for j = 1: obj.nv-1
            slope = obj.slope(j,j+1);
            pslope = -1/slope;
            if pslope == -inf
              pslope = inf;
            end
            if pslope ~= inf
              q = obj.yIntercept (j,pslope);
              eq = s2 - pslope*s1 - q;
            else
              eq = s1 - obj.vx(j);
            end
            if subs(eq,[s1,s2],[obj.vx(j+1),obj.vy(j+1)]) > 0    
              eq = -eq;
            end
            NC(j,1) = eq;
            if pslope ~= inf
              q = obj.yIntercept (j+1,pslope);
              eq = s2 - pslope*s1 - q;
            else
              eq = s1 - obj.vx(j+1);
            end
            if subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0   
              eq = -eq;
            end
            NC(j,2) = eq;
          end
          j = obj.nv;
          slope = obj.slope(j,1);
          pslope = -1/slope;
          if pslope == -inf
            pslope = inf;
          end
          if pslope ~= inf
            q = obj.yIntercept (j,pslope);
            eq = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(j);
          end
          if subs(eq,[s1,s2],[obj.vx(1),obj.vy(1)]) > 0   
            eq = -eq;
          end
          NC(j,1) = eq;
          if pslope ~= inf
            q = obj.yIntercept (1,pslope);
            eq  = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(1);
          end
          if subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0   
            eq = -eq;
          end
          NC(j,2) = eq;
        end
     end
     
        
end