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

        function testDerivatives (testCase)
            x=sym('x');
            y=sym('y');
            f1 = functionF(x*y);
            testCase.verifyEqual(isAlways(f1.dfdx(x).getF == y),true);
            testCase.verifyEqual(isAlways(f1.dfdx(y).getF == x),true);
            f2 = functionF(x,y);
            testCase.verifyEqual(isAlways(f2.dfdx(x).getF == 1/y),true);
            testCase.verifyEqual(isAlways(f2.dfdx(y).getF == -x/y^2),true);
            f3 = functionF(x^2);
            testCase.verifyEqual(isAlways(f3.dfdx(x).getF == 2*x),true);
            testCase.verifyEqual(isAlways(f3.dfdx(y).getF == 0),true);
            f4 = functionF(x^2,y);
            testCase.verifyEqual(isAlways(f4.dfdx(x).getF == 2*x/y),true);
            testCase.verifyEqual(isAlways(f4.dfdx(y).getF == -x^2/y^2),true);
            testCase.verifyEqual(isAlways(f1.limit(x, 1).getF == y),true);
            
        end

        % error in isConst - check where it is used then fix
        function testInquiry (testCase)
            x=sym('x');
            y=sym('y');
            f1 = functionF(x*y);
            testCase.verifyEqual(f1.isPolynomial,true);
            testCase.verifyEqual(f1.isQuad,true);
            testCase.verifyEqual(f1.isLinear,false);
            testCase.verifyEqual(f1.isConst,false);

            f1 = functionF(x,y);
            testCase.verifyEqual(f1.isPolynomial,false);
            testCase.verifyEqual(f1.isQuad,false);
            testCase.verifyEqual(f1.isLinear,false);
            testCase.verifyEqual(f1.isConst,false);

            f1 = functionF(x^3);
            testCase.verifyEqual(f1.isPolynomial,true);
            testCase.verifyEqual(f1.isQuad,false);
            testCase.verifyEqual(f1.isLinear,false);
            testCase.verifyEqual(f1.isConst,false);

             f1 = functionF(x+7*y);
            testCase.verifyEqual(f1.isPolynomial,true);
            testCase.verifyEqual(f1.isQuad,false);
            testCase.verifyEqual(f1.isLinear,true);
            testCase.verifyEqual(f1.isConst,false);

            f1 = functionF(7);
            testCase.verifyEqual(f1.isPolynomial,true);
            testCase.verifyEqual(f1.isQuad,false);
            testCase.verifyEqual(f1.isLinear,false);
            %testCase.verifyEqual(f1.isConst,true);
        end
    end
  
end