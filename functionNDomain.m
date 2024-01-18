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
             disp("Function")
             obj.f.printL;
             disp("Domain")
             obj.d.print;
         end

         function printL(objL)
             for i = 1:size(objL,2)
                 objL(i).print;
             end
         end

         function printLatex (obj)
             disp(" ")
             disp("Function")
             disp(" ")
             obj.f.printLLatex;
             disp("Domain")
             disp(" ")
             obj.d.printLatex;
         end

         function printLLatex(objL)
             for i = 1:size(objL,2)
                 objL(i).printLatex;
             end
         end
         function printM(objL)
             colorList = ["red","blue","yellow","green","purple","cyan","orange","brown","crimson", "pink","tan","aquamarine","navy", "PaleGreen"];
             fprintf("display(inequal({");
             f = objL(1).f.f;
             j = 1;
             for i = 1:size(objL,2)
                 if objL(i).f.f == f
                   objL(i).d.printMaple;
                   fprintf(",");
                   
                 else
                   fprintf("},x=-15..15,y=-15..15,color=[");
                   fprintf(colorList(j)) ;
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf("],nolines),inequal({")  ;
                   f = objL(i).f.f;
                   j = j+1;
                   objL(i).d.printMaple;
                   fprintf(","); 
                 end
                 
             end
             fprintf("},x=-15..15,y=-15..15,color=[");
                   fprintf(colorList(j)) ;
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf("],nolines))")  ;
             fprintf("\n");
             
         end

          function plotDomain(obj)
             figure;
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
             n = 0
             f = obj (1).f
             c = colors(mod(n,6)+1)
             for i =1:size(obj,2)
                i
                if (f.f ~= obj (i).f.f)
                  n = n + 1
                  c = colors(mod(n,6)+1)
                  f = obj (i).f
                end
                obj (i).d.plot;
                textR = "R"+num2str(i);
                textR="";
                obj (i).d.plotRegionC(textR,c);
             end
          end

          % for closed regions and objL1 = objL2
          function [objL,index] = times (objL1, objL2)
             n = 0;
             objL=functionNDomain.empty();
             for i = 1:size(objL1,2)
               for j = i+1:size(objL1,2)
                 rf = objL1(i).d + objL2(j).d;
                 if isempty(rf)
                   continue
                 end
                 rf = rf.simplify;
                 if rf.nv <= 2
                     rf = region.empty;
                 end
                 if isempty(rf)
                     disp("empty")
                   continue
                 end
                 n = n + 1;
                 objL(n) = functionNDomain([objL1(i).f(1), objL2(j).f(1)],rf);
                 index(n,1:2) = [i,j];

                 % objL1(i).d.print
                 % objL2(j).d.print
                 % rf.print
                 r = objL1(i).d - rf;
                 if ~ isempty(r)
                   n = n + 1;  
                   objL(n) = functionNDomain([objL1(i).f(1)],r);
                   index(n,1) = [i];
                 end

                 r = objL1(j).d - rf;
                 if ~ isempty(r)
                   n = n + 1;  
                   objL(n) = functionNDomain([objL2(j).f(1)],r);
                   index(n,1) = [j];
                 end
               end
             end
                      
          end


          function [objR,index2] = maximumPC(objL, index) %, f, r2)
           n = 0;
           for i = 1:size(objL,2)
               i,size(objL(i).f,2)
             if size(objL(i).f,2) == 1
               n = n + 1;
               objR(n) = objL(i);
               index2(n) = index(i);
               continue;
             end
             [l, fmax, ind, lSing] = objL(i).d.maximum(objL(i).f);
             if lSing
                continue
             end
             l
             if l
               n = n + 1;
               objR(n) = functionNDomain([fmax],objL(i).d);
               index2(n) = index(i,ind);
               continue
             end  
             ineqs = objL(i).d.splitmax3 (objL(i).f(1),objL(i).f(2));
             ineqs1 = sym.empty ;          
             for k = 1: size(objL(i).d.ineqs,2)
               ineqs1(k) = objL(i).d.ineqs(k).f;
             end
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(1).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1 = d1.simplify;
             n = n + 1;
             objR(n) = functionNDomain([objL(i).f(1)],d1);
             index2(n) = index(i,1);
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(2).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1 = d1.simplify;
             n = n + 1;
             objR(n) = functionNDomain([objL(i).f(2)],d1);
             index2(n) = index(i,2);
           end
           if n == 0
               objR = functionNDomain.empty();
              return
           end
           return
           
         end
          % will work only when entire R2 is covered
         function objL = mtimes (objL1, objL2)
             n = 0;
             objL=functionNDomain.empty();
             for i = 1:size(objL1,2)
               for j = 1:size(objL2,2)
                 rf = objL1(i).d + objL2(j).d;
                 if isempty(rf)
                   continue
                 end
                 rf = rf.simplifyOpenRegion;
                 if isempty(rf)
                     disp("empty")
                   continue
                 end
                 n = n + 1;
                 objL(n) = functionNDomain([objL1(i).f(1), objL2(j).f(1)],rf);
               end
             end
                      
         end

         function objR2 = maximumP(objL, lmerge) %, f, r2)
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
             d1.print
             d1 = d1.simplifyOpenRegion;
             disp("Further subdivision")
             objL(i).f(1).print
             objL(i).f(2).print
             d1.print
             d1.printMaple
             n = n + 1;
             objR(n) = functionNDomain([objL(i).f(1)],d1);
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(2).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1.print
             d1 = d1.simplifyOpenRegion;
             d1.print
             d1.printMaple
             n = n + 1;
             objR(n) = functionNDomain([objL(i).f(2)],d1);
           end
           if n == 0
              return
           end
           if ~lmerge
               objR2 = jSort(objR);
               return
            end
            % disp("b4 merge")
            % objR.printL
           objR2 = mergeL(objR);
            % disp("aft merge")
            % objR.printL
         end

         function [objL2,index,lCh] = maxEqDom(objL)  
           lCh = false;  
           ia(1) = 1;
           n = 0;
           for i = 1:size(objL,2)
              marked(i) = false;
           end
           ja = [];
           % ja has indices of all equal functions , ia by col no
           for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(objL(i).d == objL(j).d)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
            end
            ia
            ja
            for i = 1:size(objL,2)
              marked(i) = false;
            end
            m = 0;
            for i = 1:size(objL,2)
              if marked(i)
                continue;
              end
              if ia(i) == ia(i+1)
                m = m + 1;
                objL2(m) = objL(i);
                index(m) = i;
                marked(i) = true;
                continue;    
              end

              
              for j=ia(i):ia(i+1)-1
                [l, fmax, ind, lSing] = objL(i).d.maximum([objL(i).f, objL(ja(j)).f]);
                l
                lCh = true;
                marked(ja(j)) = true;
               if l
                 m = m + 1;
                 objL2(m) = functionNDomain(fmax,objL(i).d);
                 if ind == 1
                   index(m) = i;
                 else
                   index(m) = ja(j);  
                 end
                 
                 
               else
                 ineqs = objL(i).d.splitmax2 (objL(i).f, objL(ja(j)).f);
                 ineqs1 = sym.empty;
                 for k = 1: size(objL(i).d.ineqs,2)
                    ineqs1(k) = objL(i).d.ineqs(k).f;
                 end
                            
                 ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(1);
                 d1 = region(ineqs1,objL(i).d.vars);
                 d1 = d1.simplify; 
                 m = m + 1;
                 objL2(m) = functionNDomain(objL(i).f,d1);
                 index(m) = i;
                 
                 ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(2);
                 d1 = region(ineqs1,objL(i).d.vars);
                 d1 = d1.simplify; 
                 m = m + 1;
                 objL2(m) = functionNDomain(objL(ja(j)).f,d1);
                 index(m) = ja(j);
                 
               end  
                  
                 
              end
            end
             
          end
          
      
         

         function objL2 = jSort(objL)  
          ia(1) = 1;
          n = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          ja = [];
          % ja has indices of all equal functions , ia by col no
          for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(objL(i).f.f == objL(j).f.f)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
          %ia
          %ja
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          m = 0;
          for i = 1:size(objL,2)
            if marked(i)
                continue;
            end
            m = m + 1;
            objL2(m) = objL(i);
            marked(i) = true;
            for j=ia(i):ia(i+1)-1
                m = m + 1;
                objL2(m) = objL(ja(j));
                marked(ja(j)) = true;
            end
           end
             
          end
          
          function [objL2,index] = mergeL(objL)  
          ia(1) = 1;
          n = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          ja = [];
          % ja has indices of all equal functions , ia by col no
          for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(objL(i).f.f == objL(j).f.f)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
          %ia
          %ja
          m = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          for i = 1:size(objL,2)
            if  marked(i)
                continue
            end
            if (ia(i) == ia(i+1))
                if marked(i)
                   continue
                end
                 
                m = m + 1;
           %     objL(i).d.print
                objL2(m) = objL(i);
                marked(i) = true;
                index(m) = i;
            else
                % get common boundary and merge
                % make groups and add 
               r = objL(i).d;
             %  disp ('first r')
             %  r.print
               lmerge = true;
               while lmerge
                 lmerge = false;
                 for j=ia(i):ia(i+1)-1
                   
                   if marked(ja(j))
                       continue
                   end
                   [l,r] = r.merge (objL(ja(j)).d);
                   if l
                     marked(ja(j)) = true;
                     lmerge = true;
              %       ja(j)
              %       objL(ja(j)).d.print
              %       r.print
                   end
                 end
               end
               m = m + 1;
               objL2(m) = functionNDomain([objL(i).f],r);
              % r.print
               marked(i) = true;
               index(m) = i;
               for j=ia(i):ia(i+1)-1
                 if marked(ja(j))
                   continue
                 end
                 r = objL(ja(j)).d ;
               %  disp('r2')
               %  ja(j)
               %  r.print
                 lmerge = true;
                 while lmerge
                   lmerge = false;
                   for k=j+1:ia(i+1)-1
                     if marked(ja(k))
                       continue
                     end
                     [l,r] = r.merge (objL(ja(k)).d);
                     
                     if l
                       marked(ja(k)) = true;
                       lmerge = true;
                %       ja(k)
                %       objL(ja(k)).d.print
                %       r.print
                     end
                   end
                 end
                 m = m + 1;
                % r.print
                 objL2(m) = functionNDomain([objL(i).f],r);
                 marked(i) = true;
                 index(m) = i;
                 marked(ja(j)) = true;

                 
               end

            end
        end
      end
     end
end 