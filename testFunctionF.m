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

        function testTangent(testCase)
            x = sym('x');
            y = sym('y');
            f2 = symbolicFunction(148*x - 196*y + (x + 7*y)^2 - 684)
            %f2.tangent(3.159091,-1.022727)
            f3 = symbolicFunction(x - (9*y)/5 - 5)
            [tx,ty] = solve(f3.f==0,f2.f==0)
            
            t = f2.tangent(tx,ty)
            t = t.normalize1
            f2 = -f2  
            t = f2.tangent(tx,ty)
            t = t.normalize1
        end
        
    end

    
end