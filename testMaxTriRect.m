classdef testMaxTriRect < matlab.unittest.TestCase

    properties
        PTri
        PRect
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
            q(1) = plq_1piece(d(3),f);
            testCase.PRect = plq(q);

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
              testCase.PTri.pieces(i).print
              %testCase.PTri.pieces(i).plot
            end 
            %return
            testCase.PTri.nPieces=2;
            testCase.PTri = testCase.PTri.maximumInFirstPairs;
            
            testCase.PTri=testCase.PTri.maximumP;
           % return
            size(testCase.PTri.maxf)
            disp("tri")
            for i = 1:size(testCase.PTri.maxf)
                i
            testCase.PTri.maxd(i).print
            testCase.PTri.maxf(i)
            end
            
            return
            %% merge is reordering  %%
            %% fix merge when only one vertex and edge going to infinity
   
           
             %   testCase.PS.pieces(1).plot
            i = 1
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).convexEnvelope;
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).conjugate;
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).intersectionConjugateDomain;
            testCase.PRect.pieces(i)=testCase.PRect.pieces(i).maximum;
            disp("rect")
%             testCase.PRect.pieces(i).maxd(1).print
%             testCase.PRect.pieces(i).maxf(1)
%             testCase.PRect.pieces(i).maxd(7).print
%             testCase.PRect.pieces(i).maxf(7)
            for i = 1:size(testCase.PRect.pieces(i).maxf,1)
                i
            testCase.PRect.pieces(1).maxd(i).print
            testCase.PRect.pieces(1).maxf(i)
           end
return            
            
            
            % temp code in case something merges
            %[nmaxf,nmaxd] = testCase.PRect.merge(testCase.PRect.maxf,testCase.PRect.maxd)   
            
            % Currently merge needs to be fixed so checking directly 

            size(testCase.PRect.pieces(1).maxf,1)
            size(testCase.PTri.maxf,1)
            %testCase.verifyEqual(size(testCase.PRect.pieces(1).maxf,1)==size(testCase.PTri.maxf,1), true);
%             for i = 1: size(testCase.PRect.pieces(1).maxf,1)
%                 mark1(i) = false;
%                 mark2(i) = false;
%             end
%             for i = 1: size(testCase.PRect.pieces(1).maxf,1)
%               if mark1(i)
%                   continue;
%               end
%               for j = 1: size(testCase.PRect.pieces(1).maxf,1)
%                   if mark2(j)
%                     continue;
%                   end
%                   if ~ (testCase.PRect.pieces(1).maxf(i) == testCase.PTri.maxf(j))
%                       continue
%                   end
%                   if ~ (testCase.PRect.pieces(1).maxd(i) == testCase.PTri.maxd(j))
%                       continue
%                   end
%                   mark1(i) = true;
%                   mark2(j) = true;
%               end
%             end
%             nSingleton = 0;
%             for i = 1: size(testCase.PRect.pieces(1).maxf,1)
%                 if mark1(i)
%                     continue
%                 end
%                 if (testCase.PRect.pieces(1).maxd(i).nv==1)
%                     nSingleton = nSingleton + 1;
%                     mark1(i)=true;
%                 end
%             end
%             for i = 1: size(testCase.PRect.pieces(1).maxf,1)
%                 if mark2(i)
%                     continue
%                 end
%                 if (testCase.PTri.maxd(i).nv==1)
%                     nSingleton = nSingleton - 1;
%                     mark2(i)=true;
%                 end
%             end
%             testCase.verifyEqual(nSingleton==0, true);
%             for i = 1: size(testCase.PRect.pieces(1).maxf,1)
%                 if mark1(i)
%                     continue
%                 end
%                  i, 1
%                  testCase.PRect.pieces(1).maxf(i).print
%                  testCase.PRect.pieces(1).maxd(i).print
%                 
%             end
%             for i = 1: size(testCase.PRect.pieces(1).maxf,1)
%                  if mark2(i)
%                      continue
%                  end
%                  i, 2
% 
%                  testCase.PTri.maxf(i).print
%                 testCase.PTri.maxd(i).print
%                 
%             end
% %return
%              mark1
%              mark2
%             %testCase.verifyEqual(all(mark1)&all(mark2), true);
%             %return
%             for i =1:size(testCase.PTri.maxf,1)
%                 testCase.PTri.maxf(i).print
%                 testCase.PTri.maxd(i).print
%             end
            %return
              figure;
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
              n = 0
              f = testCase.PTri.maxf(1,1)
              c = colors(mod(n,6)+1)

               for i =1:size(testCase.PTri.maxf,1)
                     i
                   if (f.f ~= testCase.PTri.maxf(i,1).f)  
                       n = n + 1
                       c = colors(mod(n,6)+1)
                       f = testCase.PTri.maxf(i,1);
                   end
                   testCase.PTri.maxd(i,1).plot;
                   textR = "R"+num2str(i);
                   textR="";
                   testCase.PTri.maxd(i,1).plotRegionC(textR,c);
                   

%           
             end
             
            % write routine to check all domain and functions
        end

 end

    
end