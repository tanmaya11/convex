classdef plq
  properties
      nPieces = 0;
      pieces = plq_1piece.empty();
      maxConjugate = functionNDomain.empty();
      biconjugate = functionNDomain.empty();
      biconjfia = [];
  end
  % 11 methods
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
          obj.maxConjugate.printL
          obj.biconjugate.printL
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
          %return
          obj.maxConjugate.printLLatex
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
      
        for i=1: obj.nPieces
          obj.pieces(i)=obj.pieces(i).convexEnvelope;

          
        end
        return
      end

      
      function obj = conjugate(obj)
      
        for i=1: obj.nPieces
          obj.pieces(i)=obj.pieces(i).convexEnvelope;

          obj.pieces(i)=obj.pieces(i).conjugate;
        end
        return
      end



      function obj = maximum(obj)
      
        for i=1: obj.nPieces
          i
           % if i == 2
           %     return
           % end
         tic
          obj.pieces(i)=obj.pieces(i).convexEnvelope;
          toc
%return
          disp("ConvexEnvelope")
          obj.pieces(i).print
      %    obj.pieces(i).envelope.printL
        %  continue
      %    return
        tic
             obj.pieces(i)=obj.pieces(i).conjugate;
             %continue %return
             toc
              disp("Conjugate")
              tic
           %   obj.pieces(i).conjfia
              obj.pieces(i).conjugates.printL
              %return
              obj.pieces(i) = obj.pieces(i).maximumConjugate;
              toc
              disp("MaxConjugate") 
              %obj.pieces(i).maxConjugate.printL
              %obj.pieces(i).maxConjugate.printM
             % return
             
           %       if i == 2
           %     return
           % end

            % 
            % 
            %return
        end
        %return
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Combining conj 1 and conj4 to check that parabolic region
        % disappears

        % obj.pieces(1).conjfia(1):obj.pieces(1).conjfia(2)-1 
        % 
        %   for k = obj.pieces(1).conjfia(1):obj.pieces(1).conjfia(2)-1 
        %      obj.pieces(2).maxConjugate(k) = obj.pieces(1).conjugates(k);
        % 
        %    end
        %    obj.pieces(1).conjugates(obj.pieces(1).conjfia(1):obj.pieces(1).conjfia(2)-1 ).printM2
        %   % obj.conjugates(obj.conjfia(2):obj.conjfia(3)-1 ).printM2
        %   % disp("mConjM")
        %   %     obj.maxConjugate.printM
        %   obj.pieces(2).conjfia(2):obj.pieces(2).conjfia(3)-1 
        %   for i = 2:size(obj.pieces(2).envelope,2)+1 %size(obj.conjfia,2)-1
        %       obj.pieces(2).conjfia(i):obj.pieces(2).conjfia(i+1)-1
        %       %obj.maxConjugate.printL
        %       obj.pieces(2).maxConjugate = obj.pieces(2).maxConjugate * obj.pieces(2).conjugates(obj.pieces(2).conjfia(i):obj.pieces(2).conjfia(i+1)-1);
        %       disp("mConj")
        %       obj.pieces(2).conjugates(obj.pieces(2).conjfia(i):obj.pieces(2).conjfia(i+1)-1).printM2
        %       %obj.maxConjugate.printL
        %       %return
        %       obj.pieces(2).maxConjugate.printM2
        %       obj.pieces(2).maxConjugate = obj.pieces(2).maxConjugate.maximumP(true);
        %       disp("temp ans")
        %       obj.pieces(2).maxConjugate.printL
        %       %obj.maxConjugate.printM
        %   end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %obj.pieces(1) = obj.pieces(1).maximumConjugateTemp;
        tic
        obj = obj.maximumConjugate;
        toc
     %   obj.maxConjugate.printM;
       % obj.plotMaxConjugateDomain;
      end

      function printDomainMaple(obj)
          for i=1: obj.nPieces
            obj.pieces(i).Mprint;
          end
          %return
          obj.maxConjugate.printM;
          obj.biconjugate.printM;
          %obj.maxConjugate.printM2;
       
      end

      function obj = maximumConjugate(obj)
          % if obj.nPieces < 1
          %     return;
          % end

          for k = 1:size(obj.pieces(1).maxConjugate,2) 
             obj.maxConjugate(k) = obj.pieces(1).maxConjugate(k);
             
             % return
          end
          %disp('piece1')
          %obj.maxConjugate.printL
          %disp('piece2')
          %obj.pieces(2).maxConjugate.printL
          %return
          for j = 2:obj.nPieces
              obj.maxConjugate = obj.maxConjugate * obj.pieces(j).maxConjugate;
              %obj.maxConjugate.printL
              %disp('in maxC')
              %obj.maxConjugate.printL
              %obj.maxConjugate.printM2

              obj.maxConjugate = obj.maxConjugate.maximumP(true);
              %disp("print in maxConjugate")
              %obj.maxConjugate.printL
              %obj.maxConjugate.printM2

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

  methods % biconjugate
   % 
   %      % linear / quadratic
    function obj = biconjugateF(obj)
        % for i = 1:obj.nPieces
        %    obj.pieces(i).biconjugateP
        % end
           bc = obj.maxConjugate.conjugateOfPiecePoly
           bc.printL
           obj.biconjugate = bc.addEq;
           obj.biconjugate.printL
      
    end
  
   %      % fractional
  end
   
end
