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
              obj.pieces(i).print;
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
              if i == 1
                  return
              end
          end
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
