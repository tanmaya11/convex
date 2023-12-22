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
            p(1) = plq_1piece(d(1),f);
            testCase.PRect = plq(p);

  %          testCase.PS.plotDomain
%             
       end
    end

    methods (Test)

        function testMax (testCase)
            i = 1
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).convexEnvelope;
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).conjugate;
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).intersectionConjugateDomain;
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).maximum;
            disp("rect")
%             testCase.PRect.pieces(i).maxd(1).print
%             testCase.PRect.pieces(i).maxf(1)
%             testCase.PRect.pieces(i).maxd(7).print
%             testCase.PRect.pieces(i).maxf(7)
            for i = 1:size(testCase.PRect.pieces(i).maxf,1)
                i
            testCase.PRect.pieces(1).maxd(i).print
            testCase.PRect.pieces(1).maxf(i)
           end
           return            
            
            
             
            % write routine to check all domain and functions
        end

 end

    
end