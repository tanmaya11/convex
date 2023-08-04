classdef plq
  properties
      nPieces = 0;
      pieces = plq_1piece.empty();
  end
  methods
      function obj = plq(ps)
          obj.nPieces = 0;
          for i = 1:size(ps,2)
             obj.nPieces = obj.nPieces + 1; 
             obj.pieces(i) = ps(i);
          end
      end

      function print(obj)
          disp("")
          disp("")
          disp("")
          disp("")
          disp("")
          for i = 1:obj.nPieces
if i > 2
                  return
              end
              disp(["Piece ", num2str(i)])
              obj.pieces(i).print;
              
              %disp("extra")
              %obj.pieces(i).conjf.printL
            end
      end

      function plot(obj)
           for i = 1:obj.nPieces
               obj.pieces(i).plot
           end
      end

      
      function obj = convexEnvelope(obj)
          for i = 1:obj.nPieces
              i
              if i == 1
                  continue;
              end
              obj.pieces(i)=obj.pieces(i).convexEnvelope;
              if i == 2
                  return
              end
              disp('end')
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
             % obj.pieces(i).conjf.printL
              if i == 2
                  return
              end
          end
      end

      function [f,r] = intersectionConjugateDomain (obj)
        n = 0;
        
        for i1 = 1:obj.nPieces
          for j1 = 1:size(obj.pieces(i1).conjfia,2)-1
             for k1 = obj.pieces(i1).conjfia(j1):obj.pieces(i1).conjfia(j1+1)-1
              
                for i2 = 1:obj.nPieces
                  for j2 = 1:size(obj.pieces(i2).conjfia,2)-1
                    for k2 = obj.pieces(i2).conjfia(j2):obj.pieces(i2).conjfia(j2+1)-1
                       if(i1 == i2 & k2 < obj.pieces(i1).conjfia(j1+1))
                         continue
                       end 
                 
                     k1
                     k2
                     if(k1 <= k2)
                   %    disp('Conjugate Domain Intersection')
                   %    disp([k1,k2])
                       if (k1 == 6  & k2 == 10)
                         [l,r1] = intersection3(obj.pieces(i1).conjd(k1), obj.pieces(i2).conjd(k2), true);
                       else
                           [l,r1] = intersection3(obj.pieces(i1).conjd(k1), obj.pieces(i2).conjd(k2), false);
                       end
                       if l
              %           disp("Conjugate Intersection")
                         % fill vertices of r
                         %obj.pieces(i1).conjd(k1).print
                         %obj.pieces(i2).conjd(k2).print
                         f1 = obj.pieces(i1).conjf(k1);
                         f2 = obj.pieces(i1).conjf(k2);
                 %        size(r1,2)
                         for ir = 1:size(r1,2)
                             disp('r1')
                             r1(ir).print
                           r1(ir) = r1(ir).getVertices();  
                           % Removing regions which are points
                           if r1(ir).nv == 1
                               continue;
                           end 
                           %disp('Feasible region')
                           %disp(n)
                           %r1(ir).isFeasibleWBPts
                           n = n + 1;
                           r(n) = r1(ir);
                           %r(n) = r(n).getVertices();
                    %       disp(n)
                    %       r(n).print

                           f(n,1) = f1;
                           f(n,2) = f2;
                           %return
                         end
                         %return
                         
                       %if k1 == 2 & k2 == 5
                       %    return
                       %end

                       end
                       %return
                       %obj.pieces(i1).conjd(k1).print
                       %disp('Conjugate Domain 2')
              
                       %obj.pieces(i2).conjd(k2).print
                    end
                  end
                 end
              end
            end
           
          end
          
        end
          
          
          
          %         obj.conjd(k1).intersection2(obj.conjd(k2))
          
        end

      function [maxf,maxd] = maximum(obj, f, r2)
          % 1 - eq constant
          % 2 - linear ineq const
          % 3 - linear


          % not used quadratic ineq - check that
          n = 0;
          for i = 1:size(r2,2)

              % do this if functions are linear
             disp(["i", num2str(i)])
             %r2(i).print

% check this

             %f1 = f(i,1);
             %f2 = f(i,2);
             % [l, fmax] = r2(i).maxArray (f1, f2) ;
%%%%

              [l, fmax, ind] = r2(i).maximum(f(i,:));
               if l
                 n = n + 1;
                 maxf(n) = fmax;
                 maxd(n) = r2(i);
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

      function [nmaxf,nmaxd] = merge(obj,maxf,maxd)

          ia(1) = 1;
          n = 0;
          for i = 1:size(maxf,2)
              marked(i) = false;
          end
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
          ia
          ja
          nmaxf = [];
          nmaxd = [];
          m = 0
          for i = 1:n
            if (ia(i) == ia(i+1))
                m = m + 1
                nmaxf(m) = maxf(i);
                nmaxd(m) = maxd(i);
            else
                % get common boundary and merge
                % make groups and add 
                %r = maxd(i);
                %for j=ia(i),ia(i+1)-1
                %    [l,r1] = r.merge (maxd(ja(j)))
                %    if l
                %        r1 = r;
                %    end
                %end
            end
          end
      end
      end

end
