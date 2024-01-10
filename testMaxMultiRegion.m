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
       end
    end

    methods (Test)

        function testMax (testCase)
            testCase.PRect = testCase.PRect.maximum
        end
    end

    
end