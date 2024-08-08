classdef testMaxMultiRegion < matlab.unittest.TestCase

    properties
        PTri
        PRect
        PRect2
        PRect3
        Poly
        PTri2
        PThesis
        POpen
    end

    methods (TestMethodSetup)
       function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=symbolicFunction(x*y);
            
            d(1)=domain([-5,-4;0,-4;1,3;-5,5],x,y);
            d(2)=domain([0,-4;2,0;2,1;1,3],x,y); 
            d(3)=domain([-1,1;-3,-3;-4,-3],x,y);
            d(4)=domain([1,0;3,1;2,2;0,1],x,y);
            d(5)=domain([-5,-4;0,-4;2,0;2,1;1,3;-5,5],x,y);
            

            
            p(1) = plq_1piece(d(1),f);
            %f=symbolicFunction(x^2-y^2);
            p(2) = plq_1piece(d(2),f);
            testCase.PRect = plq(p);
            
            testCase.PTri = plq([plq_1piece(d(3),symbolicFunction(x^2-y^2))]);

            testCase.PRect2 = plq([plq_1piece(d(4),symbolicFunction(x^2-y^2))]);

            p(3) = plq_1piece(d(5),f);
            testCase.Poly = plq(p(3));

            d(6)=domain([0,0;1,1;2,1;2,0],x,y);
            testCase.PRect3 = plq([plq_1piece(d(6),symbolicFunction(x*y))]);

            % Expt1
            d(7)=domain([0,0;1,1;2,0],x,y);
            d(8)=domain([1,1;2,1;2,0],x,y);
            
            % Expt8
            d(7)=domain([-5,0;2,0;0,-4;-5,-4],x,y);
            d(8)=domain([-5,5;1,3;2,1;2,0;-5,0],x,y);
            
            q(1) = plq_1piece(d(7),f);
            q(2) = plq_1piece(d(8),f);
            
            testCase.PTri2 = plq(q);

            %cube1
            d(9)=domain([-1,-1;-1,1;1,1;1,-1],x,y);
            q(3) = plq_1piece(d(9),f);
            testCase.PThesis = plq(q(3));

            %cube4


       end
    end

    methods (Test)

        function testMax (testCase)
            testCase.PRect = testCase.PRect.maximum
           % return
             testCase.PRect.print
             testCase.PRect.printDomainMaple
            %% 
            %testCase.PRect.printLatex
            
           
        end

        function testBiconjugate (testCase)
            s_1 = sym('s1');
            s_2 = sym('s2');
            f = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = s_1 + 2*s_2 + 4;
            ineq(2) = s_1 - (9*s_2)/5 - 5;
            ineq(3) = - s_2 - 5 ;
            ineq(4) = s_1 + 7*s_2 - 46;
            d = region(ineq,[s_1,s_2]);
            d = d.removeInfV;
            
            d = d.poly2orderUnbounded;
            edgeNo = d.getEdgeNosInf(d.vars)
            d.print
            return
            
            d.ineqs(edgeNo) = d.ineqs;
            d.print
            return
            p = functionNDomain([f], [d]);
            p.printL
            pc = p(1).conjugateOfPiecePoly ;

        end


        function testBiconjugate2 (testCase)
            s_1 = sym('s1');
            s_2 = sym('s2');
            f1 = symbolicFunction(-5*s_1 - 4*s_2 - 20);
            ineq(1) = s_1 + 4;
            ineq(2) =  s_2 + 5 ;
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);
            %d.print
            %return
            f2 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = s_1 - (9*s_2)/5 - 5;
            ineq(2) =  -s_2 - 5 ;
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs; 
            
            p(2) = functionNDomain(f2,d2);

            f3 = symbolicFunction(-4*s_2);
            ineq(1) = - s_1 - 4;
            ineq(2) =  (9*s_2)/5 - s_1 + 5;
            d3 = region(ineq,[s_1,s_2]);
            d3 = d3.removeInfV;
            d3 = d3.poly2orderUnbounded;
            edgeNo = d3.getEdgeNosInf(d3.vars)
            d3.ineqs(edgeNo) = d3.ineqs; 
            p(3) = functionNDomain(f3,d3);

            
            
            p.printL
            %return
            pc = p.conjugateOfPiecePoly ;

            d = pc(1).d+pc(2).d;
            d.print
            d = d + pc(3).d;





        end

        function testBiconjugate3 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s1');
            s_2 = sym('s2');
            f1 = symbolicFunction(s_1 + 3*s_2 - 1);
             

            
            ineq(1) = s_1 - 2*s_2 - 1;
            ineq(2) =  s_2/3 - s_1 + 13/3 ;
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);
            %d.print
            %return
            
            f2 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = s_1 - (5*s_2)/7 - 25/7;
            ineq(2) =  s_1 - s_2/3 - 13/3 ;

            
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs; 
            
            p(2) = functionNDomain(f2,d2);

            f3 = symbolicFunction(2*s_1+s_2-2);

           
             
            ineq(1) = -s_1 + 2*s_2 + 1;
            ineq(2) =  -s_2 + 2;
            d3 = region(ineq,[s_1,s_2]);
            d3 = d3.removeInfV;
            d3 = d3.poly2orderUnbounded;
            edgeNo = d3.getEdgeNosInf(d3.vars)
            d3.ineqs(edgeNo) = d3.ineqs; 
            p(3) = functionNDomain(f3,d3);

            f4 = symbolicFunction(2*s_1);
            
            ineq(1) = (5*s_2)/7 - s_1 + 25/7;
            
            ineq(2) =  s_2 - 2;
            d4 = region(ineq,[s_1,s_2]);
            d4 = d4.removeInfV;
            d4 = d4.poly2orderUnbounded;
            edgeNo = d4.getEdgeNosInf(d4.vars)
            d4.ineqs(edgeNo) = d4.ineqs; 
            p(4) = functionNDomain(f4,d4);
            
            
            p.printL
            %return
            pc = p.conjugateOfPiecePoly ;
           
            d = pc(1).d+pc(2).d;
            disp('1+2')
            d.print
            disp('1+2+3')
            d = d+pc(3).d;
            d.print
            disp('1+2+3+4')
            d = d+pc(4).d;
            d.print
            return




        end

        function testBiconjugate4 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s1');
            s_2 = sym('s2');

            % V1
            f1 = symbolicFunction(2*s_1);
            ineq(1) = (5*s_2)/7 - s_1 + 25/7;
            ineq(2) =  4 - 2*s_2 - s_1 ;
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);

            % V2
            f2 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = - s_1 - 2*s_2 - 4 ;
            ineq(2) =  s_1 + 2*s_2 - 4 ;
            ineq(3) = 48*s_1 - 56*s_2 + (s_1 + 2*s_2)^2 - 184;
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs;
            p(2) = functionNDomain(f2,d2);

            % E1
            f3 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq = sym.empty() 
            ineq(1) = s_1 + 2*s_2 + 4 ;
            ineq(2) =  s_1 - (9*s_2)/5 - 5 ;
            
            
            d3 = region(ineq,[s_1,s_2]);
            d3 = d3.removeInfV;
            d3 = d3.poly2orderUnbounded;
            edgeNo = d3.getEdgeNosInf(d3.vars)
            d3.ineqs(edgeNo) = d3.ineqs;
            p(3) = functionNDomain(f3,d3);

            % E2
            f4 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq = sym.empty() 
            ineq(1) = s_1 - (5*s_2)/7 - 25/7 ;
            ineq(2) =  4 - 2*s_2 - s_1 ;
            
            d4 = region(ineq,[s_1,s_2]);
            d4 = d4.removeInfV;
            d4 = d4.poly2orderUnbounded;
            edgeNo = d4.getEdgeNosInf(d4.vars)
            d4.ineqs(edgeNo) = d4.ineqs;
            p(4) = functionNDomain(f4,d4);

            % V3
            f5 = symbolicFunction(-4*s_2);
            ineq = sym.empty() 


            ineq(1) = s_1 + 2*s_2 + 4 ;
            ineq(2) = (9*s_2)/5 - s_1 + 5;
            
            d5 = region(ineq,[s_1,s_2]);
            d5 = d5.removeInfV;
            d5 = d5.poly2orderUnbounded;
            edgeNo = d5.getEdgeNosInf(d5.vars)
            d5.ineqs(edgeNo) = d5.ineqs;
            p(5) = functionNDomain(f5,d5);
            

            % E3
            f6 = symbolicFunction(0.125*s_1^2 + 0.5*s_1*s_2 + s_1 + 0.5*s_2^2 + -2*s_2 + 2);
            ineq = sym.empty() 


            ineq(1) = - s_1 - 2*s_2 - 4 ;
            ineq(2) = s_1 + 2*s_2 - 4;
            ineq(3) = 56*s_2 - 48*s_1 - (s_1 + 2*s_2)^2 + 184;


            
            d6 = region(ineq,[s_1,s_2]);
            d6 = d6.removeInfV;
            d6 = d6.poly2orderUnbounded;
            edgeNo = d6.getEdgeNosInf(d6.vars)
            d6.ineqs(edgeNo) = d6.ineqs;
            p(6) = functionNDomain(f6,d6);

            p.printL
            
            pc = p.conjugateOfPiecePoly ;
            %pc = pc.mergeL;
            pc.printL
return
            d = pc(1).d+pc(6).d;
            disp('1+6')
            d.print
            
            % disp('1+2+8')
            %  d = d+pc(8).d;
            % 
            %  d.print
              disp('1+2+10')
             d = d+pc(10).d;

             d.print
             return
            disp('1+2+6+8')
             d = d+pc(8).d;
             d.print
             return
            % d = d+pc(5).d;
            % d.print
            d = d+pc(6).d;
            d.print
            d = d+pc(7).d;
            d.print
            d = d+pc(8).d;
            d.print
            d = d+pc(9).d;
            d.print
            d = d+pc(10).d;
            d.print
            return
            
        end

        function testConvex (testCase)
            %testCase.PRect3.print
            testCase.PRect3 = testCase.PRect3.convexEnvelope
            
            %% 
            %testCase.PRect.printLatex
            %testCase.PRect3.print
           
        end
        function testConjugate (testCase)
            %testCase.PRect3.print
            testCase.PRect3 = testCase.PRect3.conjugate
            
            %% 
            %testCase.PRect.printLatex
            %testCase.PRect3.print
           
        end

        function testMaxR3 (testCase)
            testCase.PRect3.print
            testCase.PRect3 = testCase.PRect3.maximum
            
            %% 
            %testCase.PRect.printLatex
           % testCase.PRect3.print
           % testCase.PRect3.printDomainMaple
        end
        
        function testMaxT (testCase)
            %testCase.PRect3.print
            testCase.PTri2 = testCase.PTri2.maximum
            
            %% 
            %testCase.PRect.printLatex
            %testCase.PTri2.print
            testCase.PTri2.printDomainMaple
        end
        
        function testMaxP (testCase)
            testCase.Poly.print
            testCase.Poly = testCase.Poly.maximum
             testCase.Poly.print
           %  return
          
            %% 
           % testCase.Poly.printLatex
            testCase.Poly.printDomainMaple
           
           
        end

        
        % function testMax2 (testCase)
        %     disp('here')
        %     testCase.PTri = testCase.PTri.maximum
        % end

        function testMax3 (testCase)
            testCase.PRect2 = testCase.PRect2.maximum
            %testCase.PRect.printDomainMaple
            %testCase.PRect2.printLatex
           
        end

        function testMaxThesis (testCase)
            %testCase.PRect3.print
            warning('off','all') 
            testCase.PThesis = testCase.PThesis.maximum
            
            %% 
            %testCase.PRect.printLatex
           % testCase.PThesis.print
            testCase.PThesis.printDomainMaple
        end

        function testMaxThesis2 (testCase)
            %testCase.PRect3.print
            warning('off','all')

            n = 2;
            x = sym('x');
            y = sym('y');
            f=symbolicFunction(x*y);
            
            div = symbolicFunction(2,x);
            div = subs(div.f,x,n);
            del = symbolicFunction(x);
            del = subs(del.f,x,-1);
            
            %del = -1;
            for i=1:n+1
                A(i) = del;
                del = div + del;
                
            end
            A
            %return

             
            m = 0;
            dom = [];
            for i = 1:n
                for j = 1:n
                    for ki = i:i+1
                        for kj = j:j+1
                          dom = [dom;[A(ki),A(kj)]];
                        end
                    end
                    temp = dom(3,1:2);
                    dom(3,1:2) = dom(4,1:2);
                    dom(4,1:2) = temp;
                    %dom
                    m = m + 1;
                    d(m) = domain(dom,x,y);
                    %d(m).print
                    q(m) = plq_1piece(d(m),f);
                    dom = [];
                    %break
                end
                %break
            end
            testCase.PThesis = plq(q);

            testCase.PThesis = testCase.PThesis.maximum
            
            %% 
            %testCase.PRect.printLatex
            testCase.PThesis.print
            testCase.PThesis.printDomainMaple
        end


        function testOpenconvex (testCase)
            x=sym('x');
            y=sym('y');
            d = domain();
            d = d.domainEdge([y,-x-1,x-1],[x,y]);
            %d.print
            testCase.POpen = plq([plq_1piece(d,symbolicFunction(x*y))]);
            testCase.POpen = testCase.POpen.convexEnvelope
            testCase.POpen.print
        end

        
    end

    
end