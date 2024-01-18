classdef plq
  properties
      nPieces = 0;
      pieces = plq_1piece.empty();
      maxConjugate = functionNDomain.empty();
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
             disp(["Piece ", num2str(i)])
              obj.pieces(i).print;
            end
      end

      function printLatex(obj)
          disp("")
          disp("")
          disp("")
          disp("")
          disp("")
          for i = 1:obj.nPieces
             fprintf("Piece " + num2str(i) + "\n")
             disp(" ")
              obj.pieces(i).printLatex;
          end
          return
          obj.maxConjugate.printLatex
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

      


      function obj = maximum(obj)
      
        for i=1: obj.nPieces
          i
          obj.pieces(i)=obj.pieces(i).convexEnvelope;
          disp("ConvexEnvelope")
          return
         %   obj.pieces(i)=obj.pieces(i).conjugate;
            % disp("Conjugate")
            % obj.pieces(i) = obj.pieces(i).maximumConjugate;
            % disp("MaxConjugate")
            % 
            % 
        end
        return
        obj = obj.maximumConjugate;
     %   obj.maxConjugate.printM;
       % obj.plotMaxConjugateDomain;
      end

      function printDomainMaple(obj)
          for i=1: obj.nPieces
            obj.pieces(i).Mprint;
          end
          obj.maxConjugate.printM;
       
      end

      function obj = maximumConjugate(obj)
          if obj.nPieces < 1
              return;
          end
          for k = 1:size(obj.pieces(1).maxConjugate,2) 
             obj.maxConjugate(k) = obj.pieces(1).maxConjugate(k);
          end
          
          for j = 2:obj.nPieces
              obj.maxConjugate = obj.maxConjugate * obj.pieces(j).maxConjugate;
              obj.maxConjugate = obj.maxConjugate.maximumP(false);
          end
      end

      function plotMaxConjugateDomain(obj)
            
             figure;
         
         
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
             n = 0
             f = obj.maxConjugate (1).f
             c = colors(mod(n,6)+1)

             for i =1:size(obj.maxConjugate,2)
                i
                if (f.f ~= obj.maxConjugate (i).f.f)
                  n = n + 1
                  c = colors(mod(n,6)+1)
                  f = obj.maxConjugate (i).f
                end
                obj.maxConjugate (i).d.plot;
                textR = "R"+num2str(i);
                textR="";
                obj.maxConjugate (i).d.plotRegionC(textR,c);
             end
         end

  end


end
