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

          l6 = [x + 2*y + 4, x + (5*y)/2 - 1, 46 - 7*y - x, x - 2*y + 44 ] 

          testCase.w = region(l6,[x,y]);

          %l7 = [-x - 7*y - 4, x + 7*y-10, 196*y-148*x-(x+7*y)^2 + 684, 4-2*y-x ] 

          %l7 = [x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, 4-2*y-x, x - y/3-14/3,5*y/7-x+25/7 ] 

          %l7 = [-x - 7*y - 4, x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, x-9*y/5-5, x-y/3-14/3 ] 

          %l7 = [-x - 7*y - 4, x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, -x+9*y/5+5, -x-4 ] 

          %l7 = [x-y/3-14/3, 10-7*y-x]
          %l7 = [ -196*y+148*x+(x+7*y)^2 - 684, -x+9*y/5+5 ] 

          l7 = [-x-7*y-4, x+7*y-10, 148*x-196*y+(x+7*y)^2-684, -x-2*y-4, x+2*y-4, (x+2*y)^2-24*y-8*x+32, 2*y-x-44] 

          testCase.w = region(l7,[x,y]);
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

        function testremoveTangent (testCase)
            x = sym('x');
            y = sym('y');
            l7 = [ -196*y+148*x+(x+7*y)^2 - 684, x-9*y/5-5,-x-7*y-4, x+7*y-10 ] 
            testCase.w = region(l7,[x,y]);
            testCase.w.print;

            
             testCase.w = testCase.w.removeTangent(testCase.w.nv,testCase.w.vx,testCase.w.vy);
             testCase.w.print;
             return
             % tangent removed

            l7 = [ -196*y+148*x+(x+7*y)^2 - 684, x-9*y/5-5 ] 
            testCase.w = region(l7,[x,y]);
            testCase.w.print;

             testCase.w = testCase.w.removeTangent(testCase.w.nv,testCase.w.vx,testCase.w.vy);
             testCase.w.print;
             % tangent not removed


             l7 = [ 196*y-148*x-(x+7*y)^2 + 684, x-9*y/5-5 ] 
            testCase.w = region(l7,[x,y]);
            testCase.w.print;

             testCase.w = testCase.w.removeTangent(testCase.w.nv,testCase.w.vx,testCase.w.vy);
             testCase.w.print;
             % tangent not removed

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

        function testMerge(testCase)
          x = sym('x');
          y = sym('y');
          l1 = [-x - 7*y - 4, x + 7*y-10, 196*y-148*x-(x+7*y)^2 + 684, 4-2*y-x ] 
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion
                  
          l2 = [(y/3)-x+ 14/3, 10-7*y-x, y-2]
          r2 = region(l2,[x,y]);
          r2 = r2.simplifyUnboundedRegion
          l3 = [x+7*y+4,4-2*y-x]
          r3 = region(l3,[x,y]);
          r3 = r3.simplifyUnboundedRegion
          l4 = [-y/3 +x - 14/3, 10 -7*y - x, y-2, 5*y/7 -x + 25/7]
          r4 = region(l2,[x,y]);
          r4 = r4.simplifyUnboundedRegion
         
          l5 = [x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, 4-2*y-x, x - y/3-14/3,5*y/7-x+25/7 ] 
          r5 = region(l5,[x,y]);
          r5 = r5.simplifyUnboundedRegion
          r1.print
          r2.print
          [l,r] = merge(r1,r2)
          r.print

          r3.print
          [l,r] = merge(r,r3)
          r.print

           r4.print
          [l,r] = merge(r,r4)
          r.print

           r5.print
          [l,r] = merge(r,r5)
          r.print
          r.printMaple
        end

        function testMerge2(testCase)
          x = sym('x');
          y = sym('y');

          l0 = [x + 7*y - 10 , 148*x - 196*y + (x + 7*y)^2 - 684 , 4 - 2*y - x , (5*y)/7 - x + 25/7 ]
          r0 = region(l0,[x,y]);
          r0 = r0.simplifyUnboundedRegion;
          
          l1 = [x-y/3-14/3, 10-7*y-x, y-2, (5*y)/7-x+25/7] 
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          r0.print
          r1.print
          [l,r] = merge(r0,r1)
          l
          r.print

          return
          
        end

        function testMerge3(testCase)
          x = sym('x');
          y = sym('y');

          l0 = [x+2*y+4, x-(9*y)/5-5, -y-5, x+7*y-46];
          r0 = region(l0,[x,y]);
          r0 = r0.simplifyUnboundedRegion;
          
          l1 = [x-y/3-14/3, x-2*y+44, 46-7*y-x] ;
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          r0.print
          r1.print
          [l,r] = merge(r0,r1)
          l
          r.print

          return
          
       end
    end

end