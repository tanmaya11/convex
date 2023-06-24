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

      function intersectionConjugateDomain (obj)
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
                       disp('Conjugate Domain 1')
                       disp([k1,k2])
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

      function obj = maximum(obj)
          for i = 1:obj.nPieces
              obj.pieces(i)=obj.pieces(i).maxOfConjugate;
              if i == 1
                  return
              end
          end
      end
  end
end
