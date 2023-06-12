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
          for i = 1:obj.nPieces
              disp(["Piece ", num2str(i)])
              disp("Domain")
              obj.pieces(i).d.print
              fprintf("\n")
              disp("Function")
              obj.pieces(i).f.print
              fprintf("\n\n\n")
              
              disp("Convex Envelope")
              size(obj.pieces(i).envf)
              for j=1:size(obj.pieces(i).envf,2) 
                disp('Function')  
                obj.pieces(i).envf(j).print
                disp('Domain')
                obj.pieces(i).envd(j).print
              end

    
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
       
  end
end
