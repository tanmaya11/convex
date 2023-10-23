classdef plq
  properties
      nPieces = 0;
      pieces = plq_1piece.empty();
      nmaxf
      maxf=functionF.empty();
      maxd = region.empty();
  end
  methods
      function obj = plq(ps)
          obj.nPieces = 0;
          for i = 1:size(ps,2)
             obj.nPieces = obj.nPieces + 1; 
             obj.pieces(i) = ps(i);
          end
      end

      % io
      function print(obj)
          disp("")
          disp("")
          disp("")
          disp("")
          disp("")
          for i = 1:obj.nPieces
%if i > 1
%                  return
%              end
              disp(["Piece ", num2str(i)])
              obj.pieces(i).print;
              
              %disp("extra")
              %obj.pieces(i).conjf.printL
            end
      end

      function obj = rdMaxfd (obj, uNo)
            s1 = sym('s1');
            s2 = sym('s2');
            n = 0;
            line = fgetl(uNo);
            lfunc = 0;
            lregion = 0;
            nV = 0;
            
            while (ischar(line))
                if strfind(line, "Max Piece") 
                    lfunc = 1;
                else if lfunc == 1
                    n = n + 1;
                    obj.maxf(n,1) = functionF(str2sym(line));
                    disp("Function")
                    obj.maxf(n,1).print
                    lfunc = 0;
                    
               else if lregion == 1
                       
                  if nV == 0
                      nV = str2num(line);
                      iV = 0;
                      nEq = 0;
                  else if (iV < nV)
                          iV = iV + 1;
                  else if isempty(line)
                      r = region(ineq, [s1,s2]);
                      obj.maxd(n,1) = r; 
                      lregion = 0;
                      nEq = 0;
                      ineq = sym.empty;
                      nV = 0;
                      %r.print
                  else
                      nEq = nEq+1;
                      ineq(nEq) = str2sym(line);
                  end
                  end
                        
                  end
                else if isempty(line)
                   lregion = 1;
                end
                end
                end
                
                end
                line = fgetl(uNo);
                
                
                %expr = sym(exprString)
            end
            disp('out')
           % r = region(ineq, [s1,s2]);
           % obj.maxd(n,1) = r; 
           % r.print
           n
           size(obj.maxd)
            
            
        end

      function obj = rdMaxfd2 (obj, uNo)
            s1 = sym('s1');
            s2 = sym('s2');
            %uNo = fopen('../data/max3d.m','r')
            n = 0;
            line = fgetl(uNo);
            lfunc = 0;
            lregion = 0;
            nV = 0;
            
            while (ischar(line))
                if strfind(line, "Max Piece") 
                    lfunc = 1;
                else if lfunc == 1
                    nf = str2num(line);
                    n = n + 1;
                    obj.nmaxf(n) = nf;
                    for i = 1:nf
                      line = fgetl(uNo);  
                      obj.maxf(n,i) = functionF(str2sym(line));
                      %disp("Function")
                      %obj.maxf(n,i).print
                    end
                    lfunc = 0;
                    
               else if lregion == 1
                       
                  if nV == 0
                      nV = str2num(line);
                      iV = 0;
                      nEq = 0;
                  else if (iV < nV)
                      iV = iV + 1;
                  else if isempty(line)
                      r = region(ineq, [s1,s2]);
                      obj.maxd(n,1) = r; 
                      lregion = 0;
                      nEq = 0;
                      ineq = sym.empty;
                      nV = 0;
                      %r.print
                  else
                      nEq = nEq+1;
                      ineq(nEq) = str2sym(line);
                  end
                  end
                        
                  end
                else if isempty(line)
                   lregion = 1;
                end
                end
                end
                
                end
                line = fgetl(uNo);
                
                
                %expr = sym(exprString)
            end

          
            %obj.maxf(n,1).print
            %obj.maxd(n,1).print

            if nEq ~= 0
              r = region(ineq, [s1,s2]);
              obj.maxd(n,1) = r; 
              r.print
            end


            
            fclose(uNo);
            return
            size(obj.maxf,1)
            size(obj.maxd,1)
            rm = [];
            for i = 1:size(obj.maxf,1)
              if obj.maxd(i,1).hasNegativeIneqs
                  obj.maxd(i,1).print
                  
                  rm = [rm,i]
              end
            end
            rm
            %size(obj.maxf)
            obj.maxf(rm,:) = [];
            obj.maxd(rm,:) = [];
            obj.nmaxf(rm) = [];
            %obj.maxd(18,1).print
            
            %size(obj.maxf)
        end

        
        function plot(obj)
           for i = 1:obj.nPieces
               obj.pieces(i).plot
           end
      end

      function plotMaxd(obj)
          figure;
          for i =1:size(obj.maxf,1)
            obj.maxd(i,1).plotRegion;
          end
      end
      function plotDomain(obj)
           %figure;
           for i = 1:obj.nPieces
               obj.pieces(i).plotDomain
           end
      end

      
      function obj = convexEnvelope(obj)
          for i = 1:obj.nPieces
              %i
              %if i == 1 
              %    continue;
              %end
              %if i >= 2 
              %    continue;
              %end
              obj.pieces(i)=obj.pieces(i).convexEnvelope;
              %if i ==  1
              %    return
              %end
              %disp('end')
          end
      end
       
      
      
      
      function obj = conjugate(obj)
          for i = 1:obj.nPieces
             % disp('conjugate')
             % i
             % obj.pieces(i).envf.printL

              obj.pieces(i)=obj.pieces(i).conjugate;
             % disp('done')
             % size(obj.pieces(i).envf)
             % obj.pieces(i).envf.printL
             % size(obj.pieces(i).conjf)
              %obj.pieces(i).conjf.printL
              %if i == 1
              %    return
              %end
          end
      end

      function obj = intersectionConjugateDomain (obj)
        %n = 0;
        
        for i = 1:obj.nPieces
            %disp('piece')
            %i
            %obj.pieces(i).conjfia
           obj.pieces(i) = obj.pieces(i).intersectionConjugateDomain;
        end


%           for j1 = 1:size(obj.pieces(i1).conjfia,2)-1
%              for k1 = obj.pieces(i1).conjfia(j1):obj.pieces(i1).conjfia(j1+1)-1
%               
%                 for i2 = 1:obj.nPieces
%                   for j2 = 1:size(obj.pieces(i2).conjfia,2)-1
%                       j2
%                       obj.pieces(i2).conjfia(j2),obj.pieces(i2).conjfia(j2+1)
%                     for k2 = obj.pieces(i2).conjfia(j2):obj.pieces(i2).conjfia(j2+1)-1
%                        if(i1 == i2 & k2 < obj.pieces(i1).conjfia(j1+1))
%                          continue
%                        end 
%                  
%                      %k1
%                      %k2
%                      if(k1 <= k2)
%                        %disp('Conjugate Domain Intersection')
%                        %disp([k1,k2])
%                        %if (k1 == 2  & k2 == 9)
%                        %  [l,r1] = intersection3(obj.pieces(i1).conjd(k1), obj.pieces(i2).conjd(k2), true);
%                        %else
%                            [l,r1] = intersection3(obj.pieces(i1).conjd(k1), obj.pieces(i2).conjd(k2), false);
%                        %end
%                        %if k2 == 9
%                        %    obj.pieces(i1).conjd(k1).print
%                        %    obj.pieces(i2).conjd(k2).print
%                        %    l
%                        %end
%                        if l
%                            k1
%                            k2
%                       %   disp("Conjugate Intersection")
%                          % fill vertices of r
%                          %obj.pieces(i1).conjd(k1).print
%                          %obj.pieces(i2).conjd(k2).print
%                          f1 = obj.pieces(i1).conjf(k1);
%                          f2 = obj.pieces(i1).conjf(k2);
%                  %        size(r1,2)
%                          for ir = 1:size(r1,2)
%                       %       disp('r1')
%                       %       r1(ir).print
%                            r1(ir) = r1(ir).getVertices();  
%                        %    if k1 == 2 & k2 == 5
%                        %        r1(ir) = r1(ir).getVertices(true);  
%                        %        r1(ir).print
%                        %    end 
%                            % Removing regions which are points
%                            if r1(ir).nv == 1            % problem detecting R2R1 in example1
%                        %        continue;
%                            end 
%                            %disp('Feasible region')
%                            %disp(n)
%                            %r1(ir).isFeasibleWBPts
%                            n = n + 1;
%                        %    disp("k2k2")
%                        %    k1
%                        %    k2
%                            obj.pieces(ir)maxd(n) = r1(ir);
%                            %r(n) = r(n).getVertices();
%                     %       disp(n)
%                     %       r(n).print
% 
%                            obj.maxf(n,1) = f1;
%                            obj.max(n,2) = f2;
%                            %return
%                          end
%                          %return
%                          
%                        %if k1 == 2 & k2 == 5
%                        %    return
%                        %end
% 
%                        end
%                        %return
%                        %obj.pieces(i1).conjd(k1).print
%                        %disp('Conjugate Domain 2')
%               
%                        %obj.pieces(i2).conjd(k2).print
%                     end
%                   end
%                  end
%               end
%             end
%            
%           end
%           
%         end
%           
%           
%           
          %         obj.conjd(k1).intersection2(obj.conjd(k2))
          
        end

      function obj = maximum(obj) %, f, r2)



          for i=1:obj.nPieces
              obj.pieces(i) = obj.pieces(i).maximum;
          end
          return
          % 1 - eq constant
          % 2 - linear ineq const
          % 3 - linear


          % not used quadratic ineq - check that
          n = 0;
          for i = 1:size(obj.maxd,2)

              % do this if functions are linear
             disp(["i", num2str(i)])
             %r2(i).print

% check this

             %f1 = f(i,1);
             %f2 = f(i,2);
             % [l, fmax] = r2(i).maxArray (f1, f2) ;
%%%%

              [l, fmax, ind] = obj.maxd(i).maximum(obj.maxf(i,:));
               if l
                 n = n + 1;
                 maxf(n) = fmax;
                 %maxd(n) = r2(i);
                 maxd(n) = obj.maxd(i);
                 continue
               end  



             if f1.isConst & f2.isConst  
                 disp('cc')
                 disp ('should never come here')
             elseif f1.isConst & f2.isLinear
                disp('cl')
               
             elseif f1.isLinear & f2.isConst
                disp('lc')
               

             elseif f1.isLinear() & f2.isLinear() 
                disp('ll')
               
               % Find intersection point
               %[vx,vy] = solveF (f1, f2);
             elseif f1.isConst & f2.isQuad   

                 disp('cq')
               
               
             elseif f1.isQuad & f2.isConst
                 disp('qc')
               
               
             elseif f1.isLinear & f2.isQuad   
                 disp('lq')
               
               
             elseif f1.isQuad & f2.isLinear
                 disp('ql')
               
               
             elseif f1.isQuad & f2.isQuad
             disp('qq')  
               
             end 

             disp('reached here')

%               p1 = 0; 
%               p2 = 0;
%           
%               
%               if f1.isConst
%                   p1 = 1;
%               end
%               if f2.isConst
%                   
%                   p2 = 1;
%               end
%               
%               % 2. check f with ineq for obvious results
% 
%               if p1 == 0
%                   [l1,less,fIneq1, nf] = r2(i).funcIneq (f1)
%                   if l1 
%                   if nf
%                       p1 = 3;
%                   else    
%                       p1 = 2;
%                   end
%                   end
%               end
%               if p2 == 0
%                   [l2,less,fIneq1, nf] = r2(i).funcIneq (f2)
%                   if l2
%                       if nf
%                           p2 = 3;
%                       else
%                         p2 = 2;
%                       end
%                   end
%               end
%               disp(["p1", num2str(p1)])
%               disp(["p2", num2str(p2)])
%               
%               if (p1 == 1 & p2 == 1)
%                   n = n + 1;
%                   if (double(f1.f) < double(f2.f))
%                     maxf(n) = f2
%                     maxd(n) = r2(i)
%                   else
%                     maxf(n) = f1
%                     maxd(n) = r2(i)
%                   end
%                   return
%               elseif (p1 == 1 & p2 == 2)
%                   n = n + 1;
%                   if (double(f1.f) < fIneq1)
%                     maxf(n) = f2
%                     maxd(n) = r2(i)
%                   else
%                     maxf(n) = f1
%                     maxd(n) = r2(i)
%                   end
%               elseif (p1 == 2 & p2 == 3)
%                   fv1 = fIneq1;
%                   fv2 = r2(i).funcVertices (f2);
%                   lg = true;
%                   for i = 1:size(fv2,2)
%                       if (fv2(i) < fv1)
%                           lg = false;
%                           break
%                       end 
%                   end
%                   lg = true;
%                   if l1 & lg
%                       n = n + 1;
%                       disp('f2')
%                       maxf(n) = f2
%                       maxd(n) = r2(i);
%                   end
%                   for i = 1:size(fv2,2)
%                       if (fv2(i) > fv1)
%                           lg = false;
%                           break
%                       end 
%                   end
%                   
%                   if ~l1 & lg
%                       n = n + 1;
%                       maxf(n) = f1
%                       maxd(n) = r2(i);
%                   end
%               elseif (p1 == 3 & p2 == 3)    
%                 fv1 = r2(i).funcVertices (f1);
%                 fv2 = r2(i).funcVertices (f2);
%                 lg = true;
%                 for i = 1:size(fv2,2)
%                    
%                   if fv2(i).f < fv1(i).f
%                     lg = false;
%                     break
%                   end 
%                 end
%                 if lg
%                   n = n + 1;
%                   f2
%                   maxf(n) = f2;
%                   maxd(n) = r2(i);
%                 else
%                 lg = true;
%                 for i = 1:size(fv2,2)
%                   if fv1(i).f < fv2(i).f
%                     lg = false;
%                     break
%                   end 
%                 end
%                 if lg
%                   n = n + 1;
%                   f1
%                   maxf(n) = f1;
%                   maxd(n) = r2(i);
%                 end
%                 
% 
%                   % Evaluate f2 at vertices
%                   % If all max or min resolve
%               end
% 
% 
%               % if all vertices have same max f
%                             % check pointwise if one function is max
%                             % get intersection of functions in domain   
%                             
%                          
%           end
          end
          
      end

      function obj = maximumInFirstPairs(obj)
           
           for i=1:2
              for k1 = 1:size(obj.pieces(i).maxd,1)
                lc(i,k1) = false;
              end
           end
           n = 0;
           for i=1:2
             for k1 = 1:size(obj.pieces(i).maxd,1)
               for j=i+1:2
                 for k2 = 1:size(obj.pieces(j).maxd,1)
                   [l,rf] = intersection3(obj.pieces(i).maxd(k1), obj.pieces(j).maxd(k2), false);
                   if l
                       %k1, k2
                       %size(rf)
                       lc(i,k1) = true;
                       lc(j,k2) = true;
                       for irf=1:size(rf,2)
                         %  rf(i).print
                         n = n + 1;
                         obj.maxd(n,1) = rf(irf);
                         obj.maxf(n,1) = obj.pieces(i).maxf(k1);
                         obj.maxf(n,2) = obj.pieces(j).maxf(k2);
                         obj.nmaxf(n) = 2;
                         if n == 5
                         %    return
                         end

                       end

                       
                       
                   end
%               if ~l
%                   continue
                 end
               end
             end
           end
           for i=1:2
              for k1 = 1:size(obj.pieces(i).maxd,1)
                if lc(i,k1)
                    continue;
                end
                n = n + 1;
                obj.maxd(n,1) = obj.pieces(i).maxd(k1);
                obj.maxf(n,1) = obj.pieces(i).maxf(k1);
                obj.nmaxf(n) = 1;
                         
              end
           end
           
           
       end

        
        function obj = maximumInPairsAddi(obj, ind)

            for i =1:size(obj.maxf,1)
               lc(1,i) = false;
            end
           
            for k1 = 1:size(obj.pieces(ind).maxd,1)
                lc(2,k1) = false;
            end
            disp('size of p3')
            size(obj.pieces(ind).maxd,1)
           n = 0;
           for k1=1:size(obj.maxf,1)
             for k2 = 1:size(obj.pieces(ind).maxd,1)
                [l,rf] = intersection3(obj.maxd(k1,1), obj.pieces(ind).maxd(k2), false);
                
                if l
                   % k1, k2
                  lc(1,k1) = true;
                  lc(2,k2) = true;
                  %size(rf)
                  %obj.maxd(k1,1).print
                  %obj.pieces(ind).maxd(k2).print
                  
                  for irf=1:size(rf,2)
                    n = n + 1
                    maxd(n,1) = rf(irf);
                    if n == 56
                        rf(irf).isFeasibleWBPts
                   rf(irf).print
                    end
                    maxf(n,1) = obj.maxf(k1);
                    maxf(n,2) = obj.pieces(ind).maxf(k2);
                    nmaxf(n) = 2;
                    
                  end
                  %return
                end
              end
           end
           %lc(1,:)
           %lc(2,:)
           disp('n')
           n
           for i =1:size(obj.maxf,1)
               if lc(1,i) 
                 continue
               end  
               n = n + 1;
               maxd(n,1) = obj.maxd(i);
               maxf(n,1) = obj.maxf(i);
               nmaxf(n) = 1;
               
           end
           
            for k1 = 1:size(obj.pieces(ind).maxd,1)
               if lc(2,k1) 
                 continue
               end
               n = n + 1;
               maxd(n,1) = obj.pieces(ind).maxd(k1);
               maxf(n,1) = obj.pieces(ind).maxf(k1);
               nmaxf(n) = 1;
                
            end
           %Put explicit copy and check 
           obj.maxf=functionF.empty();
           obj.maxd = region.empty();
           obj.nmaxf = []

           for i = 1:n
             obj.nmaxf(i) = nmaxf(i);
             obj.maxd(i,1) = maxd(i,1);
             for k = 1:nmaxf(i)
               obj.maxf(i,k) = maxf(i,k);
             end
           end
           
           size(maxf)
           size(obj.maxf)
           size(obj.nmaxf)

           
       end


      


       % other part implemented - need to be integrated here and for
       % conjugate
       function obj = maximumP(obj) %, f, r2)
          

          n = 0;
          for i = 1:size(obj.maxf,1)
              %i, obj.nmaxf(i)

               % check size of obj.maxf(i,:) and fix
               if obj.nmaxf(i) == 1
                 n = n + 1;
                 maxf(n) = obj.maxf(i);
                 maxd(n) = obj.maxd(i,1);
                 continue;
               end
               [l, fmax, ind] = obj.maxd(i).maximum(obj.maxf(i,:));
               if l
                 n = n + 1;
                 maxf(n) = fmax;
                 %maxd(n) = r2(i);
                 maxd(n) = obj.maxd(i,1);
                 continue
               else  
                 disp('maximum : check if it reaches here')
                 i
                 % temp code
                 %if n == 4
                % figure;
                % obj.maxd(i,1).print;
                % obj.maxd(i,1).plot;
                % obj.maxd(i,1).plotRegion;
                % end
                 %continue
                 %obj.maxd(i,1).print
                 %obj.maxf(i,:).printL
                 %%%%%%%%%%

                          
                 ineqs = obj.maxd(i,1).splitmax3 (obj.maxf(i,1), obj.maxf(i,2));
                 ineqs1 = sym.empty ;          
                 for k = 1: size(obj.maxd(i,1).ineqs,2)
                   ineqs1(k) = obj.maxd(i).ineqs(k).f;
                 end
                 ineqs1(size(obj.maxd(i,1).ineqs,2)+1) = ineqs(1).f;
                 d1 = region(ineqs1,obj.maxd(i,1).vars);
                 d1 = d1.simplify(obj.maxd(i,1).vars);
                 d1.print
                 n = n + 1
                 maxf(n) = obj.maxf(i,1);
                 %maxf(n).print
                 maxd(n) = d1;

                 %if n == 5
                 %figure;
                 %d1.print;
                 %d1.plot;
                 %d1.plotRegion;
                 %end
                 %maxd(n).print
                 %maxf(n).print


                 ineqs1(size(obj.maxd(i,1).ineqs,2)+1) = ineqs(2).f;
                 %ineqs1(size(obj.maxd(i,1).ineqs,2)+1) = ineqs(2);
                 d1 = region(ineqs1,obj.maxd(i,1).vars);
                 d1 = d1.simplify(obj.maxd(i,1).vars);
                 d1.print
                 n = n + 1
                 maxf(n) = obj.maxf(i,2);
                 maxd(n) = d1;
                 %maxd(n).print
                 %maxf(n).print
                 %if n == 6
                 %figure;
                 %d1.print;
                 %d1.plot;
                 %d1.plotRegion;
                 %return
                 %end
                 
                 
                 
               end
          end
          if n == 0
              return
          end
%           disp('b4 merge' )
%           n
           [nmaxf,nmaxd] = obj.merge(maxf,maxd);
           obj.maxf=functionF.empty();
          obj.maxd = region.empty();
          %size(nmaxf,2)
%           for i =1:size(maxf,2)
%             obj.maxf(i,1) = maxf(i);
%             obj.maxd(i,1) = maxd(i);
%           end
          for i =1:size(nmaxf,2)
            obj.maxf(i,1) = nmaxf(i);
            obj.maxd(i,1) = nmaxd(i);
          end
%           
      end 

      

      function [nmaxf,nmaxd] = merge0(obj,maxf,maxd)
          ia(1) = 1;
          n = 0;
          for i = 1:size(maxf,2)
              marked(i) = false;
          end

          % ja has indices of all equal functions , ia by col no
          for i = 1:size(maxf,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(maxf,2)
                  if isAlways(maxf(i).f == maxf(j).f)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
          %nmaxf = [];
          %nmaxd = [];
          m = 0;
          for i = 1:size(maxf,2)
              marked(i) = false;
          end
        %  nmaxf= [];
        %  nmaxd= [];
          for i = 1:size(maxf,2)
              
            if  marked(i)
                continue
            end
            if (ia(i) == ia(i+1)) 
                m = m + 1;
                nmaxf(m) = maxf(i);
                nmaxd(m) = maxd(i);
            else
                % get common boundary and merge
                % make groups and add 
               r = maxd(i);
              
               for j=ia(i):ia(i+1)-1
                   marked(ja(j)) = true;
                   [l,r] = r.merge (maxd(ja(j)));
                   
                   if ~l
                     m = m + 1;
                     nmaxf(m) = maxf(i);
                     nmaxd(m) = maxd(ja(j));  
                   end
               end
               m = m + 1;
               nmaxf(m) = maxf(i);
               nmaxd(m) = r;  
                   
            end
          end
      end

      function [nmaxf,nmaxd] = merge(obj,maxf,maxd)
          disp('in merge')
          ia(1) = 1;
          n = 0;
          size(maxf,2)
          for i = 1:size(maxf,2)
              marked(i) = false;
          end

          % ja has indices of all equal functions , ia by col no
          for i = 1:size(maxf,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(maxf,2)
                  
                  if isAlways(maxf(i).f == maxf(j).f)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
          %nmaxf = [];
          %nmaxd = [];
          m = 0;
          for i = 1:size(maxf,2)
              marked(i) = false;
          end
         % ia
         % ja
        %  return
        %  nmaxf= [];
        %  nmaxd= [];
        for i = 1:size(maxf,2)
            i
            if  marked(i)
                continue
            end
            if (ia(i) == ia(i+1)) 
                m = m + 1;
                nmaxf(m) = maxf(i);
                nmaxd(m) = maxd(i);
            else
                % get common boundary and merge
                % make groups and add 
               r = maxd(i);
               if i == 57
               r.print
               end
               
               lmerge = true;
               while lmerge
                 lmerge = false;
                 for j=ia(i):ia(i+1)-1
                   
                   if marked(ja(j))
                       continue
                   end
                   if i == 57
               maxd(ja(j)).print
               end
               
         %          maxd(ja(j)).print
                   [l,r] = r.merge (maxd(ja(j)));
         %      if i == 10
         %      l
         %      end
                   
                   if l
                     marked(ja(j)) = true;
                     lmerge = true;
                   end
%                    if ~l
%                      m = m + 1;
%                      nmaxf(m) = maxf(i);
%                      nmaxe(m) = maxe(i);
%                      nmaxd(m) = maxd(ja(j));  
%                    end
                 end
               end
                 for j=ia(i):ia(i+1)-1
                   
                   if marked(ja(j))
                       continue
                   end
                   marked(ja(j)) = true;
                   m = m + 1;
                   nmaxf(m) = maxf(i);
                   nmaxd(m) = maxd(ja(j));  
                 end
               
               m = m + 1;
               nmaxf(m) = maxf(i);
               %r = r.getVertices();
               nmaxd(m) = r;  
                   
            end
          end
      end

      end

end
