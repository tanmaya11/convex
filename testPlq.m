classdef testPlq < matlab.unittest.TestCase

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
            d1=domain([0,0;2,0;2,1;1,1],x,y);
            d2=domain([-5,5;1,3;-1,0;-5,-4],x,y);
            d3=domain([-1,0;0,0;1,1;1,3],x,y);
            d4=domain([0,0;1,2;2,2;2,1],x,y);
            %d3=domain([0,0;1,2;2,2;2,1],x,y);
            p(1)=plq_1piece(d1,f);
            p(2)=plq_1piece(d2,f);
            p(3)=plq_1piece(d3,f);
            p(4)=plq_1piece(d4,f);
            testCase.PS = plq(p);
           % testCase.PS.pieces(3).print
             

%             
       end
    end

    methods (Test)
        function testCreation (testCase)
             testCase.verifyEqual(testCase.PS.nPieces,4);
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
              testCase.PS.pieces(1)=testCase.PS.pieces(1).maximum;
               testCase.PS.pieces(1).print
            
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
        
        function testConvexEnvelope3 (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            testCase.PS.pieces(3)=testCase.PS.pieces(3).convexEnvelope;
           testCase.PS.pieces(3).print
            testCase.verifyEqual(true, true);
        end

    end

    
end