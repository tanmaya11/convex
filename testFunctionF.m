classdef testFunctionF < matlab.unittest.TestCase

    methods (Test)
        function noParameter (testCase)
          f1 = functionF();
          testCase.verifyEqual(f1.getNum(),1);
          testCase.verifyEqual(f1.getDen(),1);
          testCase.verifyEqual(f1.f,1);
        end
        function oneParameter (testCase)
          x = sym('x'); 
          f1 = functionF(x^2);
          testCase.verifyEqual(f1.getNum(),x^2);
          testCase.verifyEqual(f1.getDen(),1);
          testCase.verifyEqual(f1.f,x^2);
        end
        function twoParameter (testCase)
          x = sym('x'); 
          f1 = functionF(x^2,x+1);
          testCase.verifyEqual(f1.getNum(),x^2);
          testCase.verifyEqual(f1.getDen(),x+1);
          testCase.verifyEqual(f1.f,x^2/(x+1));
        end
        
    end

    
end