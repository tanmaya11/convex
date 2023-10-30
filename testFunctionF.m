classdef testFunctionF < matlab.unittest.TestCase

    methods (Test)
        function testCreation (testCase)
            x=sym('x');
            y=sym('y');
            f=functionF(x*y);
            testCase.verifyEqual(isAlways(f.getF==x*y),true);
        end
        function testRat (testCase)
            x=sym('x');
            y=sym('y');
            num = x*y;
            den = x^2+y^2;
            f=functionF(num,den);
            testCase.verifyEqual(isAlways(f.getF==num/den),true);
            testCase.verifyEqual(isAlways(f.getNum==num),true);
            testCase.verifyEqual(isAlways(f.getDen==den),true);
        end
        function testOperations (testCase)
            x=sym('x');
            y=sym('y');
            f1 = functionF(x*y);
            f2 = functionF(x,y);
            f3 = functionF(x^2);
            f4 = functionF(x^2,y);
            g = f1+f1;
            testCase.verifyEqual(isAlways(g.getF==2*x*y),true);

            g = f1+f2;
            testCase.verifyEqual(isAlways(g.getF==x*y+x/y),true);
            testCase.verifyEqual(isAlways(g.getNum==x*y^2+x),true);
            testCase.verifyEqual(isAlways(g.getDen==y),true);
            
            g = f3+f4;
            testCase.verifyEqual(isAlways(g.getF==x^2+x^2/y),true);
            testCase.verifyEqual(isAlways(g.getNum==x^2*y+x^2),true);
            testCase.verifyEqual(isAlways(g.getDen==y),true);
            
            g = f1-f1;
            testCase.verifyEqual(isAlways(g.getF==0),true);
            testCase.verifyEqual(isAlways(g.getNum==0),true);
            testCase.verifyEqual(isAlways(g.getDen==1),true);

            g = f1-f2;
            testCase.verifyEqual(isAlways(g.getF==x*y-x/y),true);
            testCase.verifyEqual(isAlways(g.getNum==x*y^2-x),true);
            testCase.verifyEqual(isAlways(g.getDen==y),true);
            
            g = f3-f4;
            testCase.verifyEqual(isAlways(g.getF==x^2-x^2/y),true);
            testCase.verifyEqual(isAlways(g.getNum==x^2*y-x^2),true);
            testCase.verifyEqual(isAlways(g.getDen==y),true);

            % *  does not simplify
            g = f1*f1;
            testCase.verifyEqual(isAlways(g.getF==x*y*x*y),true);
            testCase.verifyEqual(isAlways(g.getNum==x^2*y^2),true);
            testCase.verifyEqual(isAlways(g.getDen==1),true);
            

            g = f1*f2;
            testCase.verifyEqual(isAlways(g.getF==x*y*x/y),true);
            testCase.verifyEqual(isAlways(g.getNum==x^2*y),true);
            testCase.verifyEqual(isAlways(g.getDen==y),true);

            g = f3*f4;
            testCase.verifyEqual(isAlways(g.getF==x^4/y),true);
            testCase.verifyEqual(isAlways(g.getNum==x^4),true);
            testCase.verifyEqual(isAlways(g.getDen==y),true);

            g = -f1;
            testCase.verifyEqual(isAlways(g.getF== -x*y),true);
            testCase.verifyEqual(isAlways(g.getNum== -x*y),true);
            testCase.verifyEqual(isAlways(g.getDen==1),true);

            g = -f3;
            testCase.verifyEqual(isAlways(g.getF== -x^2),true);
            testCase.verifyEqual(isAlways(g.getNum== -x^2),true);
            testCase.verifyEqual(isAlways(g.getDen==1),true);

        end
    end
  
end