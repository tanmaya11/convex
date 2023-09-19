classdef testRegion < matlab.unittest.TestCase

    methods(TestMethodSetup)
        % Setup for each test
        function setUpTestData(testCase)
          
        end
    end

    methods(Test)
        % Test methods

        function testMinus(testCase)
            x = sym('x');
          y = sym('y');
          l1 = [-x,-y,x+y-1];
          l2 = [-y,x+y-1,y-x];
          r = region(l1,[x,y]);
          s = region(l2,[x,y]);
          %r.print;
          %s.print;
            t = r - s;
           % t.print;
            testCase.verifyEqual(isequal(t.ineqs(1).f,-x), true);
            testCase.verifyEqual(isequal(t.ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(t.ineqs(3).f,x-y), true);
        end

        function testMinus2(testCase)
            x = sym('x');
          y = sym('y');
          l1 = [-x,-y,x+y-1];
          l2 = [-x,-y,x+y-1,y-x-0.5,-y+x-0.5];
          r = region(l1,[x,y]);
          s = region(l2,[x,y]);
          r.print;
          s.print;
            t = r - s;
            size(t)
            t(1).print;
            t(2).print;
            testCase.verifyEqual(isequal(t(1).ineqs(1).f,-x), true);
            testCase.verifyEqual(isequal(t(1).ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(t(1).ineqs(3).f,x-y+1/2), true);
            testCase.verifyEqual(isequal(t(2).ineqs(1).f,-y), true);
            testCase.verifyEqual(isequal(t(2).ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(t(2).ineqs(3).f,y-x+1/2), true);
        end
    end

end