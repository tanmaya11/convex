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
            d(1)=domain([0,0;2,0;2,1;1,1],x,y);
            d(2)=domain([-5,5;1,3;-1,0;-5,-4],x,y);
            d(3)=domain([-1,0;0,0;1,1;1,3],x,y);
            d(4)=domain([0,0;1,2;2,2;2,1],x,y);
            d(5)=domain([0,0;1,1;0,2;-1,1],x,y);
            d(6)=domain([-1,0;0,0;1,1;1,2],x,y);
            d(7)=domain([0,0;2,0;1,1],x,y);
            d(8)=domain([2,1;2,0;1,1],x,y);
            d(9)=domain([-5,5;-1,0;-5,-4],x,y);
            d(10)=domain([-5,5;1,3;-1,0],x,y);
            d(11)=domain([-1,0;0,0;1,1;],x,y);
            d(12)=domain([-1,0;1,1;1,3],x,y);
            d(13)=domain([0,0;1,1;1,3],x,y);
            d(14)=domain([-1,0;0,0;1,3],x,y);
            d(15)=domain([-1,0;0,0;0,3/2],x,y);
            d(16)=domain([1,3;1,1;0,3/2],x,y);
            d(17)=domain([1,1;0,0;0,3/2],x,y);
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
             
            for i =1:17
                p(i)=plq_1piece(d(i),f);
            end
            testCase.PS = plq(p);
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
           %testCase.PS.pieces(3).print
            testCase.verifyEqual(true, true);
        end

        function testConvexEnvelope5 (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            testCase.PS.pieces(5)=testCase.PS.pieces(5).convexEnvelope;
           testCase.PS.pieces(5).print
            testCase.verifyEqual(true, true);
        end

        function testConvexEnvelope6 (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            testCase.PS.pieces(6)=testCase.PS.pieces(6).convexEnvelope;
           testCase.PS.pieces(6).print
            testCase.verifyEqual(true, true);
        end

        function testConvexEnvelope7 (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            testCase.PS.pieces(7)=testCase.PS.pieces(7).convexEnvelope;
           testCase.PS.pieces(7).print
            testCase.verifyEqual(true, true);
        end

        function testConvexEnvelopei (testCase)
            %testCase.PS = testCase.PS.convexEnvelope();
            i = 17
            testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
           testCase.PS.pieces(i).print
            testCase.verifyEqual(true, true);
        end
    end

    
end