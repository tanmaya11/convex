classdef testPlqTri < matlab.unittest.TestCase

    properties
        PS
        PS1
    end

    methods (TestMethodSetup)
       function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=functionF(x*y);
            d(1)=domain([0,0;2,0;1,1],x,y);
            d(2)=domain([2,1;2,0;1,1],x,y);
            d(3)=domain([-5,5;-1,0;-5,-4],x,y);
            d(4)=domain([-5,5;1,3;-1,0],x,y);
            d(5)=domain([-1,0;0,0;0,3/2],x,y);
            d(6)=domain([1,3;1,1;0,3/2],x,y);
            d(7)=domain([1,1;0,0;0,3/2],x,y);
            d(8)=domain([-5,-4;-1,0;0,-4],x,y);
            d(9)=domain([-1,0;0,-4;2,0],x,y);
            d(10)=domain([1,1;2,1;1,3],x,y);
            d(11)=domain([0,0;1,0;1,1],x,y);
            %d3=domain([0,0;1,2;2,2;2,1],x,y);
            % make this a loop idiot
%             p(1)=plq_1piece(d1,f);
%             p(2)=plq_1piece(d2,f);
%             p(3)=plq_1piece(d3,f);
%             p(4)=plq_1piece(d4,f);
%             p(5)=plq_1piece(d5,f);
%             p(6)=plq_1piece(d6,f);
%             p(7)=plq_1piece(d7,f);
%             p(8)=plq_1piece(d8,f);
             
            for i =1:11
                p(i)=plq_1piece(d(i),f);
            end
            testCase.PS = plq(p);
  %          testCase.PS.plotDomain
%             
       end
    end

    methods (Test)
        function testCreation (testCase)
             testCase.verifyEqual(testCase.PS.nPieces,4);
             %testCase.PS.pieces(1).print
             testCase.verifyEqual(testCase.PS.pieces(1).checkPiece1, true);
         end
% 
        function testConvexEnvelope (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            testCase.PS.pieces(1)=testCase.PS.pieces(1).convexEnvelope;
            disp("")
            testCase.PS.pieces(1).print
            testCase.verifyEqual(testCase.PS.pieces(1).checkconvexEnvelope1, true);
        end
% 
        function testConjugate (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            %testCase.PS = testCase.PS.conjugate();
            testCase.PS.pieces(1)=testCase.PS.pieces(1).convexEnvelope;
            testCase.PS.pieces(1)=testCase.PS.pieces(1).conjugate;
           

            testCase.verifyEqual(testCase.PS.pieces(1).checkconjugate1, true);
        end
% 
         function testMaximum (testCase)
%             testCase.PS = testCase.PS.intersectionConjugateDomain;
%             testCase.PS = testCase.PS.maximum;
%             testCase.verifyEqual(true, testCase.PS.pieces(1).checkMaximum1);
              testCase.PS.pieces(1)=testCase.PS.pieces(1).convexEnvelope;
              testCase.PS.pieces(1)=testCase.PS.pieces(1).conjugate;
              testCase.PS.pieces(1)=testCase.PS.pieces(1).intersectionConjugateDomain;
              testCase.PS=testCase.PS.maximum;
             
               testCase.PS.pieces(1).print
                testCase.PS.pieces(1).plot
            
         end

         function testCreation2 (testCase)
            testCase.verifyEqual(testCase.PS.pieces(2).checkPiece2, true);
        end

        function testConvexEnvelope2 (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            testCase.PS.pieces(2)=testCase.PS.pieces(2).convexEnvelope;
            testCase.PS.pieces(2).print
            testCase.verifyEqual(testCase.PS.pieces(2).checkconvexEnvelope2, true);
        end
        function testCreation3 (testCase)
            testCase.verifyEqual(testCase.PS.pieces(3).checkPiece3, true);
        end

        function testEta3 (testCase)
            testCase.verifyEqual(testCase.PS.pieces(3).checkEta3, true);
        end
        
        
        function testConvexEnvelopei (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            % 1 checked
            % 2 checked
            i = 8;
            %testCase.PS.pieces(i).print
            for i = 1:1
              testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
              testCase.PS.pieces(i)=testCase.PS.pieces(i).conjugate;
              testCase.PS.pieces(i)=testCase.PS.pieces(i).intersectionConjugateDomain;
              testCase.PS.pieces(i)=testCase.PS.pieces(i).maximum;
              testCase.PS.pieces(i).print
            end  
             %   testCase.PS.pieces(1).plot
            testCase.verifyEqual(true, true);
        end

        function testConjugateForGivenEx (testCase)
            i = 11;
            %testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
            %testCase.PS.pieces(i).print
            testCase.PS.pieces(i)=testCase.PS.pieces(i).setconvexEnvelope();
            testCase.PS.pieces(i).print
            testCase.PS.pieces(i)=testCase.PS.pieces(i).conjugate;
            testCase.PS.pieces(i).print
            
            testCase.verifyEqual(true, true);
        end


        function testMax (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            % 1 checked
            % 2 checked
            %testCase.PS.pieces(i).print
            for i = 1:3
              testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
              testCase.PS.pieces(i)=testCase.PS.pieces(i).conjugate;
              
              testCase.PS.pieces(i)=testCase.PS.pieces(i).intersectionConjugateDomain;
              testCase.PS.pieces(i)=testCase.PS.pieces(i).maximum;
              if i == 3
                %testCase.PS.pieces(i).print
                %testCase.PS.pieces(i).plot
                %testCase.PS.pieces(i).plotRegion;
              end 
            end 
            
            %return

            testCase.PS.nPieces=3;
            disp('maximumInPairs')
            testCase.PS = testCase.PS.maximumInFirstPairs
            testCase.PS=testCase.PS.maximumP;
            %testCase.PS.plotMaxd;
            

%             % write to file
%             uNo = fopen('op','w');
%             for i =1:size(testCase.PS.maxf,1)
%                    testCase.PS.maxf(i,1).fprint(uNo);
%                    testCase.PS.maxd(i,1).fprint(uNo);
%              end
%             fclose(uNo);




             for i = 3:testCase.PS.nPieces
                 testCase.PS = testCase.PS.maximumInPairsAddi (i)
                 disp('maximumInPairsAddi')
                 size(testCase.PS.maxf,1)
                 testCase.PS=testCase.PS.maximumP;
                 disp('max')
                 size(testCase.PS.maxf,1)
                 
             end
             %return
             figure;
             for i =1:size(testCase.PS.maxf,1)
                     i
                   testCase.PS.maxf(i,1).print;
                   testCase.PS.maxd(i,1).print;
                   testCase.PS.maxd(i,1).plot;
                   testCase.PS.maxd(i,1).plotRegion;
%           
             end
             return
%       
             %   testCase.PS.pieces(1).plot
            
            testCase.verifyEqual(true, true);
        end

        function testMaxRecreate (testCase)
            testCase.verifyEqual(true, true);
        end
    end

    
end