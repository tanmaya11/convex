classdef testMaxMultiRegion < matlab.unittest.TestCase

    properties
        PTri
        PRect
        PRect2
        Poly
    end

    methods (TestMethodSetup)
       function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=symbolicFunction(x*y);
            
            d(1)=domain([-5,-4;0,-4;1,3;-5,5],x,y);
            d(2)=domain([0,-4;2,0;2,1;1,3],x,y); 
            d(3)=domain([-1,1;-3,-3;-4,-3],x,y);
            d(4)=domain([1,0;3,1;2,2;0,1],x,y);
            d(5)=domain([-5,-4;0,-4;2,0;2,1;1,3;-5,5],x,y);
            
            p(1) = plq_1piece(d(1),f);
            %f=symbolicFunction(x^2-y^2);
            p(2) = plq_1piece(d(2),f);
            testCase.PRect = plq(p);
            
            testCase.PTri = plq([plq_1piece(d(3),symbolicFunction(x^2-y^2))]);

            testCase.PRect2 = plq([plq_1piece(d(4),symbolicFunction(x^2-y^2))]);

            p(3) = plq_1piece(d(5),f);
            testCase.Poly = plq(p(3));
       end
    end

    methods (Test)

        function testMax (testCase)
            testCase.PRect = testCase.PRect.maximum
            testCase.PRect.printDomainMaple
            %% 
            testCase.PRect.printLatex
            %testCase.PRect.print
           
        end

        function testMaxP (testCase)
            testCase.Poly = testCase.Poly.maximum
             testCase.Poly.print
            % return
          
            %% 
           % testCase.Poly.printLatex
            testCase.Poly.printDomainMaple
           
           
        end

        
        % function testMax2 (testCase)
        %     disp('here')
        %     testCase.PTri = testCase.PTri.maximum
        % end

        function testMax3 (testCase)
            testCase.PRect2 = testCase.PRect2.maximum
            %testCase.PRect.printDomainMaple
            %testCase.PRect2.printLatex
           
        end

    end

    
end