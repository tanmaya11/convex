classdef testPlqTri2 < matlab.unittest.TestCase

    properties
        PS
        PS1
    end

    methods (TestMethodSetup)
       function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=functionF(x*y);
            d(1)=domain([-5,-4;0,-4;-5,5],x,y);
            d(2)=domain([-5,5;0,-4;1,3],x,y);
            d(3)=domain([0,-4;2,1;1,3],x,y);
            d(4)=domain([0,-4;2,0;2,1],x,y);
             
            for i =1:4
                p(i)=plq_1piece(d(i),f);
            end
            testCase.PS = plq(p);
  %          testCase.PS.plotDomain
%             
       end
    end

    methods (Test)
%         function testCreation (testCase)
%              testCase.verifyEqual(testCase.PS.nPieces,4);
%              %testCase.PS.pieces(1).print
%              testCase.verifyEqual(testCase.PS.pieces(1).checkPiece1, true);
%          end
% % 
%         function testConvexEnvelope (testCase)
%             %testCase.PS = testCase.PS.convexEnvelope();
%             testCase.PS.pieces(1)=testCase.PS.pieces(1).convexEnvelope;
%             disp("")
%             testCase.PS.pieces(1).print
%             testCase.verifyEqual(testCase.PS.pieces(1).checkconvexEnvelope1, true);
%         end
% % 
%         function testConjugate (testCase)
%             %testCase.PS = testCase.PS.convexEnvelope();
%             %testCase.PS = testCase.PS.conjugate();
%             testCase.PS.pieces(1)=testCase.PS.pieces(1).convexEnvelope;
%             testCase.PS.pieces(1)=testCase.PS.pieces(1).conjugate;
%            
% 
%             testCase.verifyEqual(testCase.PS.pieces(1).checkconjugate1, true);
%         end
% % 
%          function testMaximum (testCase)
% %             testCase.PS = testCase.PS.intersectionConjugateDomain;
% %             testCase.PS = testCase.PS.maximum;
% %             testCase.verifyEqual(true, testCase.PS.pieces(1).checkMaximum1);
%               testCase.PS.pieces(1)=testCase.PS.pieces(1).convexEnvelope;
%               testCase.PS.pieces(1)=testCase.PS.pieces(1).conjugate;
%               testCase.PS.pieces(1)=testCase.PS.pieces(1).intersectionConjugateDomain;
%               testCase.PS=testCase.PS.maximum;
%              
%                testCase.PS.pieces(1).print
%                 testCase.PS.pieces(1).plot
%             
%          end
% 
%          function testCreation2 (testCase)
%             testCase.verifyEqual(testCase.PS.pieces(2).checkPiece2, true);
%         end
% 
%         function testConvexEnvelope2 (testCase)
%             %testCase.PS = testCase.PS.convexEnvelope();
%             testCase.PS.pieces(2)=testCase.PS.pieces(2).convexEnvelope;
%             testCase.PS.pieces(2).print
%             testCase.verifyEqual(testCase.PS.pieces(2).checkconvexEnvelope2, true);
%         end
%         function testCreation3 (testCase)
%             testCase.verifyEqual(testCase.PS.pieces(3).checkPiece3, true);
%         end
% 
%         function testEta3 (testCase)
%             testCase.verifyEqual(testCase.PS.pieces(3).checkEta3, true);
%         end
%         
%         
%         function testConvexEnvelopei (testCase)
%             %testCase.PS = testCase.PS.convexEnvelope();
%             % 1 checked
%             % 2 checked
%             i = 8;
%             %testCase.PS.pieces(i).print
%             for i = 1:1
%               testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
%               testCase.PS.pieces(i)=testCase.PS.pieces(i).conjugate;
%               testCase.PS.pieces(i)=testCase.PS.pieces(i).intersectionConjugateDomain;
%               testCase.PS.pieces(i)=testCase.PS.pieces(i).maximum;
%               testCase.PS.pieces(i).print
%             end  
%              %   testCase.PS.pieces(1).plot
%             testCase.verifyEqual(true, true);
%         end
% 
%         function testConjugateForGivenEx (testCase)
%             i = 11;
%             %testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
%             %testCase.PS.pieces(i).print
%             testCase.PS.pieces(i)=testCase.PS.pieces(i).setconvexEnvelope();
%             testCase.PS.pieces(i).print
%             testCase.PS.pieces(i)=testCase.PS.pieces(i).conjugate;
%             testCase.PS.pieces(i).print
%             
%             testCase.verifyEqual(true, true);
%         end


        function testMax (testCase)
            uNo = fopen('op','w');
            %testCase.PS.pieces(1).d.polygon.fprint(uNo)
            for i = 4:4
              fprintf(uNo, "Piece " + num2str(i) + "\n")
              testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
              %testCase.PS.pieces(i).print
              return
              fprintf(uNo, "convexEnvelope " + num2str(i) + "\n")
               
              
              testCase.PS.pieces(i)=testCase.PS.pieces(i).conjugate;
              fprintf(uNo, "conjugate " + num2str(i) + "\n")
                 testCase.PS.pieces(i).print
                testCase.PS.pieces(i).plot
              continue
              testCase.PS.pieces(i)=testCase.PS.pieces(i).intersectionConjugateDomain;
              fprintf(uNo, "intersectionConjugateDomain " + num2str(i) + "\n")
              testCase.PS.pieces(i)=testCase.PS.pieces(i).maximum;
              fprintf(uNo, "maximum " + num2str(i) + "\n")
              %if i == 3
                testCase.PS.pieces(i).print
                testCase.PS.pieces(i).plot
                testCase.PS.pieces(i).plotRegion;
              %end 
            end 
            
            return

            testCase.PS.nPieces=2;
            disp('maximumInPairs')
            testCase.PS = testCase.PS.maximumInFirstPairs
            
 %           return
            testCase.PS=testCase.PS.maximumP;
            %testCase.PS.plotMaxd;
            

%             % write to file
%             uNo = fopen('op','w');
%             for i =1:size(testCase.PS.maxf,1)
%                    testCase.PS.maxf(i,1).fprint(uNo);
%                    testCase.PS.maxd(i,1).fprint(uNo);
%              end
%             fclose(uNo);

             fprintf(uNo, "Max12\n ")
             for i =1:size(testCase.PS.maxf,1)
                % s = num2str(i);
                % s = "Max Piece " + num2str(i) + "\n";
                   fprintf(uNo, "Max Piece " + num2str(i) + "\n");
                   testCase.PS.maxf(i,1).fprint(uNo);
                   fprintf(uNo, "\n")
                   testCase.PS.maxd(i,1).fprint(uNo);
                   fprintf(uNo, "\n")
                 end 
             fclose(uNo)
             return
             for i = 3:testCase.PS.nPieces
                 testCase.PS = testCase.PS.maximumInPairsAddi (i)
                 disp('maximumInPairsAddi')
                 size(testCase.PS.maxf,1)
                 testCase.PS=testCase.PS.maximumP;
                 disp('max')
                 size(testCase.PS.maxf,1)
                 fprintf(uNo, "Max " + num2str(i) + "\n")
                 for i =1:size(testCase.PS.maxf,1)
                   testCase.PS.maxf(i,1).fprint(uNo);
                   fprintf(uNo, "\n")
                   testCase.PS.maxd(i,1).fprint(uNo);
                   fprintf(uNo, "\n")
                 end
             end
             %return
             figure;
             for i =1:size(testCase.PS.maxf,1)
                     i
                   %testCase.PS.maxf(i,1).print;
                   %testCase.PS.maxd(i,1).print;
                   testCase.PS.maxd(i,1).plot;
                   testCase.PS.maxd(i,1).plotRegion;
%           
             end
             
           
             %   testCase.PS.pieces(1).plot
            
            testCase.verifyEqual(true, true);
        end

        % check then split into read expr + read region
        function testMaxRecreate (testCase)
           %if false 
%            uNo = fopen('../data/max12r.m','r')
%            %uNo = fopen('../data/max123.m','r')
%            
%             uNo = fopen('../data/max9r.m','r')
%             testCase.PS = testCase.PS.rdMaxfd(uNo);
%             fclose(uNo);


           i = 10;
           %fprintf(uNo, "Piece " + num2str(i) + "\n")
           testCase.PS.pieces(i)=testCase.PS.pieces(i).convexEnvelope;
          
           %fprintf(uNo, "convexEnvelope " + num2str(i) + "\n")
           testCase.PS.pieces(i)=testCase.PS.pieces(i).conjugate;
           %fprintf(uNo, "conjugate " + num2str(i) + "\n")
           testCase.PS.pieces(i)=testCase.PS.pieces(i).intersectionConjugateDomain;
           %fprintf(uNo, "intersectionConjugateDomain " + num2str(i) + "\n")
           testCase.PS.pieces(i)=testCase.PS.pieces(i).maximum;
            testCase.PS.pieces(i).print
            testCase.PS.pieces(i).plot
            return
%            %fprintf(uNo, "maximum " + num2str(i) + "\n")
           testCase.PS.nPieces=i;
           %testCase.PS = testCase.PS.rdMaxfd;


           testCase.PS = testCase.PS.maximumInPairsAddi (i)
           disp('maximumInPairsAddi')
           
           uNo = fopen('op3d','w');
           for i =1:size(testCase.PS.maxf,1)
              fprintf(uNo, "Max Piece " + num2str(i) + "\n")
              fprintf(uNo, num2str(testCase.PS.nmaxf(i)) + "\n")
              for k = 1:testCase.PS.nmaxf(i)
                testCase.PS.maxf(i,k).fprint(uNo);
              end
              fprintf(uNo, "\n")
              testCase.PS.maxd(i,1).fprint(uNo);
              fprintf(uNo, "\n")
           end
           
           fclose(uNo);
            testCase.verifyEqual(true, true);
           return
             
           
        end

         function testMaxRecreate2 (testCase)
            
           
          
           uNo = fopen('../data/max10dr.m','r')
           %uNo = fopen('../data/max3r.m','r')
           testCase.PS = testCase.PS.rdMaxfd2(uNo);
           disp('rd done')
           
           size(testCase.PS.maxf,1)
           testCase.PS=testCase.PS.maximumP;
           %return

           disp('max')
           size(testCase.PS.maxf,1)
           uNo = fopen('op3','w');
           
           %fprintf(uNo, "Max " + num2str(i) + "\n")
           for i =1:size(testCase.PS.maxf,1)
              fprintf(uNo, "Max Piece " + num2str(i) + "\n")
              testCase.PS.maxf(i,1).fprint(uNo);
              fprintf(uNo, "\n")
              testCase.PS.maxd(i,1).fprint(uNo);
              fprintf(uNo, "\n")
           end
           fclose(uNo);
             
           testCase.verifyEqual(true, true);
        end
    
        function testMaxRecreate3 (testCase)
            
           

          %%%%
          
           %uNo = fopen('../data/max3r.m','r')
           %uNo = fopen('../data/max4r.m','r')
           uNo = fopen('../data/max10r.m','r')
           %uNo = fopen('../data/max12r.m','r')
           %uNo = fopen('../data/max12t.m','r')
             
           testCase.PS = testCase.PS.rdMaxfd(uNo);
           fclose(uNo);
           %%%%
           colors = ['b', 'r', 'g', 'm', 'c', 'y'];
           n = 0
           f = testCase.PS.maxf(1,1)
           c = colors(mod(n,6)+1)
           for i =1:size(testCase.PS.maxf,1)
           end  
           figure;
             for i =1:size(testCase.PS.maxf,1)
                     i
                   if (f.f ~= testCase.PS.maxf(i,1).f)  
                       n = n + 1
                       c = colors(mod(n,6)+1)
                       f = testCase.PS.maxf(i,1);
                   end
                   testCase.PS.maxf(i,1).print;
                   testCase.PS.maxd(i,1).print;
                   testCase.PS.maxd(i,1).plot;
                   textR = "R"+num2str(i);
                   textR="";
                   testCase.PS.maxd(i,1).plotRegionC(textR,c);
                   

%           
             end
             n
           testCase.verifyEqual(true, true);
        end
    end

    
end