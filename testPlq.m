classdef testPlq < matlab.unittest.TestCase

    properties
        PS
    end

    methods (TestMethodSetup)
        function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=functionF(x*y);
            d1=domain([0,0;2,0;2,1;1,1],x,y);
            d2=domain([-5,5;1,3;-1,0;-5,-4],x,y);
            d3=domain([0,0;1,2;2,2;2,1],x,y);
            p(1)=plq_1piece(d1,f);
            p(2)=plq_1piece(d2,f);
            p(3)=plq_1piece(d3,f);
            testCase.PS = plq(p);
       end
    end

    methods (Test)
        function testCreation (testCase)
%             x=sym('x');
%             y=sym('y');
%             f=functionF(x*y);
%             d1=domain([0,0;2,0;2,1;1,1],x,y);
%             d2=domain([-5,5;1,3;-1,0;-5,-4],x,y);
%             d3=domain([0,0;1,2;2,2;2,1],x,y);
%             p(1)=plq_1piece(d1,f);
%             p(2)=plq_1piece(d2,f);
%             p(3)=plq_1piece(d3,f);
%             PS = plq(p);

            testCase.verifyEqual(testCase.PS.nPieces,3);
            testCase.verifyEqual(testCase.PS.pieces(1).checkPiece1, true)
            
        end

        function testConvexEnvelope (testCase)
            

            
            
            testCase.PS = testCase.PS.convexEnvelope();
            testCase.verifyEqual(testCase.PS.pieces(1).checkconvexEnvelope1, true)
            
            
        end
    end

    
end