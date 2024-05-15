classdef testFunctionF < matlab.unittest.TestCase

    methods (Test)
        function noParameter (testCase)
          f1 = symbolicFunction();
          
          testCase.verifyEqual(f1.getNum(),0);
          testCase.verifyEqual(f1.getDen(),1);
          testCase.verifyEqual(f1.f,0);
        end
        function oneParameter (testCase)
          x = sym('x'); 
          f1 = symbolicFunction(x^2);
          testCase.verifyEqual(f1.getNum(),x^2);
          testCase.verifyEqual(f1.getDen(),sym(1));
          testCase.verifyEqual(f1.f,x^2);
        end
        function twoParameter (testCase)
          x = sym('x'); 
          f1 = symbolicFunction(x^2,x+1);
          testCase.verifyEqual(f1.getNum(),x^2);
          testCase.verifyEqual(f1.getDen(),x+1);
          testCase.verifyEqual(f1.f,x^2/(x+1));
        end
        
    end

    
end