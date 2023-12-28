classdef testMaxMultiRegion < matlab.unittest.TestCase

    properties
        PTri
        PRect
    end

    methods (TestMethodSetup)
       function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=functionF(x*y);
            
            d(1)=domain([-5,-4;0,-4;1,3;-5,5],x,y);
            d(2)=domain([0,-4;2,0;2,1;1,3],x,y); 
            p(1) = plq_1piece(d(1),f);
            p(2) = plq_1piece(d(2),f);
            testCase.PRect = plq(p);

  %          testCase.PS.plotDomain
%             
       end
    end

    methods (Test)

        function testMax (testCase)
            for i=1:2
                i
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).convexEnvelope;
            disp("ConvexEnvelopw")
             testCase.PRect.pieces(i)=testCase.PRect.pieces(i).conjugate;
             disp("Conjugate")
             testCase.PRect.pieces(i)=testCase.PRect.pieces(i).intersectionConjugateDomain;
             disp("domain intersection")
             testCase.PRect.pieces(i)=testCase.PRect.pieces(i).maximum;
             disp("max")
             testCase.PRect.pieces(i).print
             %testCase.PRect.pieces(i).plot

            end
%             return
%             figure;
%              colors = ['b', 'r', 'g', 'm', 'c', 'y'];
%               n = 0
%               ii = 1
%               f = testCase.PRect.pieces(ii).maxf(1)
%               c = colors(mod(n,6)+1)
% 
%                for i =1:size(testCase.PRect.pieces(ii).maxf,1)
%                      i
%                    if (f.f ~= testCase.PRect.pieces(ii).maxf(i).f)  
%                        n = n + 1
%                        c = colors(mod(n,6)+1)
%                        f = testCase.PRect.pieces(ii).maxf(i)
%                    end
% %                    testCase.PRect.maxf(i,1).print;
% %                    testCase.PRect.maxd(i,1).print;
%                    testCase.PRect.pieces(ii).maxd(i).plot;
%                    textR = "R"+num2str(i);
%                    textR="";
%                    testCase.PRect.pieces(ii).maxd(i).plotRegionC(textR,c);
%                    
% 
% %           
%              end
% 
%            return
%            disp("rect")
%             testCase.PRect.pieces(i).maxd(1).print
%             testCase.PRect.pieces(i).maxf(1)
%             testCase.PRect.pieces(i).maxd(7).print
%             testCase.PRect.pieces(i).maxf(7)
%             testCase.PRect.nPieces=2;
             testCase.PRect = testCase.PRect.maximumInFirstPairs;
%             
             testCase.PRect=testCase.PRect.maximumP;
%            % return
%             size(testCase.PRect.maxf)
%             disp("final")
            for i = 1:size(testCase.PRect.maxf)
                i
            testCase.PRect.maxd(i).print
            testCase.PRect.maxf(i)
            end
            figure;
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
              n = 0
              f = testCase.PRect.maxf(1,1)
              c = colors(mod(n,6)+1)

               for i =1:size(testCase.PRect.maxf,1)
                     i
                   if (f.f ~= testCase.PRect.maxf(i,1).f)  
                       n = n + 1
                       c = colors(mod(n,6)+1)
                       f = testCase.PRect.maxf(i,1);
                   end
                   testCase.PRect.maxf(i,1).print;
                   testCase.PRect.maxd(i,1).print;
                   testCase.PRect.maxd(i,1).plot;
                   textR = "R"+num2str(i);
                   textR="";
                   testCase.PRect.maxd(i,1).plotRegionC(textR,c);
                   

%           
             end

           return            
            
            
             
            % write routine to check all domain and functions
        end

 end

    
end