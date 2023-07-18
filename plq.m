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
if i > 1
                  return
              end
              disp(["Piece ", num2str(i)])
              obj.pieces(i).print;
              
              %disp("extra")
              %obj.pieces(i).conjf.printL
            end
      end

      function obj = convexEnvelope(obj)
          for i = 1:obj.nPieces
              obj.pieces(i)=obj.pieces(i).convexEnvelope;
              if i == 1
                  return
              end
          end
      end
       
      function obj = conjugate(obj)
          for i = 1:obj.nPieces
              obj.pieces(i)=obj.pieces(i).conjugate;
              %obj.pieces(i).conjf.printL
              if i == 1
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
                 
                     
                     if(k1 <= k2)
                %       disp('Conjugate Domain Intersection')
                %       disp([k1,k2])
                       if (k1 == 2 )
                         [l,r1] = intersection3(obj.pieces(i1).conjd(k1), obj.pieces(i2).conjd(k2), false);
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
                             %disp('r1')
                             %r1(ir).print
                           n = n + 1;
                           r(n) = r1(ir);
                           r(n) = r(n).getVertices();
                           %r(n).print

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
          for i = 1:7 %size(r2,2)

              % do this if functions are linear
             disp(["i", num2str(i)])
             %r2(i).print
             f1 = f(i,1);
             f2 = f(i,2);
             if f1.isConst & f2.isConst  
              %   disp('cc')
               if f1.f < f2.f
                 n = n + 1;
                 maxf(n) = f2;
                 maxd(n) = r2(i);
                 continue
               else
                 n = n + 1;
                 maxf(n) = f1;
                 maxd(n) = r2(i);
                 continue
               end
             elseif f1.isConst & f2.isLinear
                % disp('cl')
               [l, fmax] = r2(i).maxArray (f1, f2) ;
               if l
                 n = n + 1;
                 maxf(n) = fmax;
                 maxd(n) = r2(i);
                 continue
               end  
               

             elseif f1.isLinear() & f2.isLinear() 
                % disp('ll')
               % evaluate at boundary  if one is max on boundary we r done
               [l, fmax] = r2(i).maxArray (f1, f2) ;
               if l
                 n = n + 1;
                 maxf(n) = fmax;
                 maxd(n) = r2(i);
                 continue
               end  
               
               % Find intersection point
               %[vx,vy] = solveF (f1, f2);
             end 



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
end
end