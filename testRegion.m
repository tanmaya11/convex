classdef testRegion < matlab.unittest.TestCase

    properties
        r
        s
        t
        u
        v
        w
    end
    
    methods(TestMethodSetup)
        % Setup for each test
        function setUpTestData(testCase)
          x = sym('x');
          y = sym('y');
          l1 = [-x,-y,x+y-1,x-2,y-1];
          testCase.r = region(l1,[x,y]);
          l2 = [-y,x+y-1,y-x,x-1,y-0.5];
          testCase.s = region(l2,[x,y]);
          l3 = [y-x,y-x^2,-y,y+x^2];
          testCase.t = region(l3,[x,y]);
          l4 = [y-x^2,y^2-x];
          testCase.u = region(l4,[x,y]);
          testCase.v = region([-x,x+1],[x,y]);

          %l5 = [-x - 7*y - 4, x + 7*y - 10, 148*x - 196*y + (x + 7*y)^2 - 684,-x - 4,(9*y)/5 - x + 5] 

          %l5 = [-x - 7*y - 4, x + 7*y - 10, 148*x - 196*y + (x + 7*y)^2 - 684,-(9*y)/5 + x - 5] ;

          l5 = [x+2*y+4, x+5*y/2-1, x+7*y-46]

        
          testCase.w = region(l5,[x,y]);
        end
    end

    methods(Test)
        % Test methods

        function testCreation(testCase)
           x = sym('x');
           y = sym('y');
          % l1 = [-x,-y,x+y-1];
          % l2 = [-y,x+y-1,y-x];
          % r = region(l1,[x,y]);
          % s = region(l2,[x,y]);
          % % r.print;
          % s.print;
          % r.printMaple;
          % s.printMaple;

          testCase.verifyEqual(isequal(testCase.r.ineqs(1).f,-x), true);
            testCase.verifyEqual(isequal(testCase.r.ineqs(2).f,-y), true);
            testCase.verifyEqual(isequal(testCase.r.ineqs(3).f,x+y-1), true);
            testCase.verifyEqual(isequal(testCase.r.vx,[0,0,1]), true);
            testCase.verifyEqual(isequal(testCase.r.vy,[0,1,0]), true);

            testCase.verifyEqual(isequal(testCase.s.ineqs(1).f,-y), true);
            testCase.verifyEqual(isequal(testCase.s.ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(testCase.s.ineqs(3).f,y-x), true);
            testCase.verifyEqual(isequal(testCase.s.vx,[1,0,0.5]), true);
            testCase.verifyEqual(isequal(testCase.s.vy,[0,0,0.5]), true);
        end

        function testslopeAtVertex(testCase)
            m = testCase.r.slopeAtVertex([1,2],[0,0])
          
            m = testCase.s.slopeAtVertex([1,3],[0,0])

            testCase.t.print
            m = testCase.t.slopeAtVertex([1,2],[0,0])
            m = testCase.t.slopeAtVertex([1,2],[1,1])
            
        end
        function testsimplifyUnboundedRegion(testCase)
            %Fix when slopes are equal
            % testCase.r.print;
            % testCase.r = testCase.r.simplifyUnboundedRegion;
            % testCase.r.print;

            % testCase.s.print;
            % testCase.s = testCase.s.simplifyUnboundedRegion;
            % testCase.s.print;
            % 
            % testCase.t.print;
            % testCase.t = testCase.t.simplifyUnboundedRegion;
            % testCase.t.print;
            % 
            % 
            % testCase.u.print;
            % testCase.u = testCase.u.simplifyUnboundedRegion;
            % testCase.u.print;
            % 
            % testCase.v.print;
            % testCase.v = testCase.v.simplifyUnboundedRegion;
            % testCase.v.print;
            % 
            testCase.w.print;

            testCase.w = testCase.w.simplifyUnboundedRegion;
            testCase.w.print;

        end
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