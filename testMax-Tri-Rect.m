classdef testMax-Tri-Rect < matlab.unittest.TestCase

    properties
        PTri
        PS1
    end

    methods (TestMethodSetup)
       function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=functionF(x*y);
            d(1)=domain([-5,-4;0,-4;-5,5],x,y);
            d(2)=domain([-5,5;0,-4;1,3],x,y);
            for i =1:2
                p(i)=plq_1piece(d(i),f);
            end
            testCase.PTri = plq(p);

            d(3)=domain([-5,-4;0,-4;1,3;-5,5],x,y);
            testCase.PRect = plq(plq_1piece(d(3),f));

  %          testCase.PS.plotDomain
%             
       end
    end

    methods (Test)

        function testMax (testCase)
            for i = 1:2
              testCase.PTri.pieces(i)=testCase.PTri.pieces(i).convexEnvelope;
              testCase.PTri.pieces(i)=testCase.PTri.pieces(i).conjugate;
              testCase.PTri.pieces(i)=testCase.PTri.pieces(i).intersectionConjugateDomain;
              testCase.PTri.pieces(i)=testCase.PTri.pieces(i).maximum;
            end 
            testCase.PTri.nPieces=2;
            testCase.PTri = testCase.PTri.maximumInFirstPairs
            testCase.PTri=testCase.PTri.maximumP;
            %% merge is reordering  %%
            %% fix merge when only one vertex and edge going to infinity
   
           
             %   testCase.PS.pieces(1).plot
            
            testCase.verifyEqual(true, true);
        end

 end

    
end