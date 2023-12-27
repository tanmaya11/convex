classdef region
    % ineqs always stored as <= 0
    properties
        ineqs=functionF;
        nv;
        vx;
        vy;
        vars;
    end

    methods  % testing
         function l = checkPiece1 (obj)
             l = false;
             
             x = sym('x');
             y = sym('y');
             if obj.vars ~= [sym('x'), sym('y')]
                 return
             end
             if obj.nv ~= 4
                 return
             end
             if obj.vx ~= [1,2,2,0]
                 return
             end
             if obj.vy ~= [1,1,0,0]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==y-1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==x-2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==-y);
                 return
            end
            if ~ isAlways (obj.ineqs(4).f==y-x);
                 return
            end
             l = true;
             return



         end

         function l = checkPiece2 (obj)
             l = false;
             
             x = sym('x');
             y = sym('y');
             %obj.print
             if obj.vars ~= [sym('x'), sym('y')]
                 return
             end
             if obj.nv ~= 4
                 return
             end
             if obj.vx ~= [-5  1  -1  -5]
                 return
             end
             if obj.vy ~= [5  3  0  -4]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==x/3 + y - 10/3);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==(3*x)/2 - y + 3/2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==x - y + 1);
                 return
            end
            if ~ isAlways (obj.ineqs(4).f==- x - 5);
                 return
            end
             l = true;
             return



         end

         function l = checkPiece3 (obj)
             l = false;
             
             x = sym('x');
             y = sym('y');
             if obj.vars ~= [sym('x'), sym('y')]
                 return
             end
             if obj.nv ~= 4
                 return
             end
             if obj.vx ~= [ 1  1  0  -1]
                 return
             end
             if obj.vy ~= [3  1  0  0  
]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==x - 1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==x - y);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==-y);
                 return
            end
            if ~ isAlways (obj.ineqs(4).f==y - (3*x)/2 - 3/2);
                 return
            end
             l = true;
             return



         end

         function l = checkConvexDomain11 (obj)
             l = false;
             
             x = sym('x');
             y = sym('y');
             if obj.vars ~= [x, y]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [0,2,0]
                 return
             end
             if obj.vy ~= [0,0,1]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==-y);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==y-x);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==x+y-2);
                 return
            end
             l = true;
             return
         end

         function l = checkConvexDomain12 (obj)
             l = false;
             %obj.print
             x = sym('x');
             y = sym('y');
             if obj.vars ~= [x, y]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [1  2  2 ]
                 return
             end
             if obj.vy ~= [1  1  0 ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==y - 1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==2 - y - x);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==x - 2);
                 return
            end
             l = true;
             return
         end

         function l = checkConvexDomain21 (obj)
             l = false;
             %obj.print
             x = sym('x');
             y = sym('y');
             if obj.vars ~= [x, y]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [1  -5  -1 ]
                 return
             end
             if obj.vy ~= [3  5  0  ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==x + 3*y - 10);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==3*x - 2*y + 3);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==- 5*x - 4*y - 5);
                 return
            end
             l = true;
             return
         end

         function l = checkConvexDomain22 (obj)
             l = false;
             %obj.print
             x = sym('x');
             y = sym('y');
             if obj.vars ~= [x, y]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [-5  -1  -5 ]
                 return
             end
             if obj.vy ~= [-4  0  5   ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==x - y + 1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==- x - 5);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==5*x + 4*y + 5);
                 return
            end
             l = true;
             return
         end

         
         function l = checkConjugateDomain11 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 5
                 return
             end
             if obj.vx ~= [0,-intmax,-intmax,0,-intmax]
                 return
             end
             if obj.vy ~= [0,intmax,-intmax,-intmax,0]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s1+s2);
                 return
            end
             l = true;
             return
         end

         function l = checkConjugateDomain12 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [0.5,0,intmax]
                 return
             end
             if obj.vy ~= [1.5,0,-intmax]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1+s2-2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==-s1-s2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==(s1+s2)^2/8-s1);
                 return
            end
            
            l = true;
             return
         end

         function l = checkConjugateDomain13 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [5.000000e-01  2147483647  2147483647]
                 return
             end
             if obj.vy ~= [1.500000e+00  2147483647  1.500000e+00]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==2 - s2 - s1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s2 - s1 - 1);
                 return
            end
            
            l = true;
             return
         end


         function l = checkConjugateDomain14 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [0  2147483647  0]
                 return
             end
             if obj.vy ~= [0  -2147483647  -2147483647]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==-s1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s1 + s2);
                 return
            end
            
            l = true;
             return
         end
         function l = checkConjugateDomain15 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 2
                 return
             end
             if obj.vx ~= [5.000000e-01  5.000000e-01]
                 return
             end
             if obj.vy ~= [1.500000e+00  2147483647]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==2 - s2 - s1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s1 - s2 + 1);
                 return
            end
            
            l = true;
             return
         end

         function l = checkConjugateDomain16 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [5.000000e-01  0  -2147483647]
                 return
             end
             if obj.vy ~= [1.500000e+00  0  2147483647]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 + s2 - 2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==- s1 - s2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==s1 - (s1 + s2)^2/8);
                 return
            end
            
            l = true;
             return
         end

         function l = checkConjugateDomain17 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 4
                 return
             end
             if obj.vx ~= [1  2147483647  1  2147483647]
                 return
             end
             if obj.vy ~= [2  2147483647  2147483647  2]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==1 - s1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==2 - s2);
                 return
            end
            
            l = true;
             return
            
         end

         function l = checkConjugateDomain18 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 4
                 return
             end
             if obj.vx ~= [1  -2147483647  1  -2147483647]
                 return
             end
             if obj.vy ~= [2  2147483647  2147483647  2]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 - s2 + 1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s1 - 1);
                 return
            end
            
            l = true;
             return
            
         end

         function l = checkConjugateDomain19 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 5
                 return
             end
             if obj.vx ~= [1  2147483647  -2147483647  1  2147483647]
                 return
             end
             if obj.vy ~= [2  -2147483647  -2147483647  -2147483647  2]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s2 - 2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s2 - s1 - 1);
                 return
            end
            
            l = true;
             return



            
         end

         function l = checkMaxDomain11 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [-5.000000e-01  -2147483647  -2147483647]
                 return
             end
             if obj.vy ~= [5.000000e-01  2147483647  5.000000e-01]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 + s2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s1 - s2 + 1);
                 return
            end
            
            l = true;
             return
         end


         function l = checkMaxDomain12 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 5
                 return
             end
             if obj.vx ~= [0  -5.000000e-01  -2147483647  0  -5.000000e-01]
                 return
             end
             if obj.vy ~= [0  5.000000e-01  -2147483647  -2147483647  -2147483647]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s1 + s2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==s2 - s1 - 1);
                 return
            end
            l = true;
             return
         end

         function l = checkMaxDomain13 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 1
                 return
             end
             if obj.vx ~= [5.000000e-01]
                 return
             end
             if obj.vy ~= [1.500000e+00 ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 + s2 - 2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==- s1 - s2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==s1 - s2 + 1);
                 return
            end
            if ~ isAlways (obj.ineqs(4).f==(s1 + s2)^2/8 - s1);
                 return
            end
            l = true;
             return
         end

         function l = checkMaxDomain14 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [0  5.000000e-01  2147483647]
                 return
             end
             if obj.vy ~= [0  1.500000e+00  -2147483647 ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 + s2 - 2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==- s1 - s2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==s2 - s1 - 1 );
                 return
            end
            if ~ isAlways (obj.ineqs(4).f==(s1 + s2)^2/8 - s1);
                 return
            end
            l = true;
             return
         end

         function l = checkMaxDomain15 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [1  2147483647  2147483647 ]
                 return
             end
             if obj.vy ~= [2  2147483647  2 ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s2 - s1 - 1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==2 - s2);
                 return
            end
           
            l = true;
             return
         end

         function l = checkMaxDomain16 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [0  2147483647  0]
                 return
             end
             if obj.vy ~= [0  -2147483647  -2147483647 ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==-s1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==s1 + s2);
                 return
            end
           
            l = true;
             return
         end

         function l = checkMaxDomain17 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 2
                 return
             end
             if obj.vx ~= [1  1]
                 return
             end
             if obj.vy ~= [2  2147483647   ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 - s2 + 1);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==1 - s1);
                 return
            end
           
            l = true;
             return
         end

         function l = checkMaxDomain18 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [-2147483647  -5.000000e-01  5.000000e-01]
                 return
             end
             if obj.vy ~= [2147483647  5.000000e-01  1.500000e+00   ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 + s2 - 2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==- s1 - s2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==s1 - s2 + 1 );
                 return
            end
            if ~ isAlways (obj.ineqs(4).f==(s1 + s2)^2/8 - s1);
                 return
            end
            l = true;
             return
         end

         function l = checkMaxDomain19 (obj)
             l = false;
             
             s1 = sym('s1');
             s2 = sym('s2');
             if obj.vars ~= [s1, s2]
                 return
             end
             if obj.nv ~= 3
                 return
             end
             if obj.vx ~= [-5.000000e-01  0  5.000000e-01]
                 return
             end
             if obj.vy ~= [5.000000e-01  0  1.500000e+00   ]
                 return
             end
             
            
             if ~ isAlways (obj.ineqs(1).f==s1 + s2 - 2);
                 return
             end
            if ~ isAlways (obj.ineqs(2).f==- s1 - s2);
                 return
            end
            if ~ isAlways (obj.ineqs(3).f==s2 - s1 - 1 );
                 return
            end
            if ~ isAlways (obj.ineqs(4).f==s1 - (s1 + s2)^2/8);
                 return
            end
            l = true;
             return
         end
    end


%  
     methods
         function obj = region(fs, vars) %, not)
        
            if nargin == 0
                return
            end
            m = size(fs,1);
            n = size(fs,2);
            nineq = 0;
            for i = 1:m
              for j = 1:n
                f = functionF(fs(i,j));        
                if (~f.isZero)
                  nineq = nineq + 1;  
                  obj.ineqs(nineq) = f;
                end  
              end
            end
            % add code to get unique ineqs
            obj.vars = vars;  
            %disp('b4')
            %obj.print
            obj = obj.normalize1;
            obj = obj.unique;
            %disp('aft')
            %obj.print
            %obj = obj.removeDenominator;
            %disp('b4 get vertices')
            obj = obj.getVertices  ;
            if obj.nv == 0
                obj = region.empty;
            end
            % put simplify in here
         end

         
         function obj = regionWPts(obj, vx, vy, x, y)
             n = 0;
             for i = 1: size(vx,2)-1
                 if vx(i) == vx(i+1)
                   n = n + 1;
                   ineq(n) = x-vx(i)  ;
                 else
                   m = (vy(i+1)-vy(i))/(vx(i+1)-vx(i));
                   c = vy(i)-m*vx(i);
                   n = n + 1;
                   ineq(n) = y-m*x-c
                 end 
             end
             
             if vx(1) == vx(end)
               n = n + 1;
               ineq(n) = x-vx(1)  
             else
               m = (vy(end)-vy(1))/(vx(end)-vx(1));
               c = vy(1)-m*vx(1);
               n = n + 1;
               ineq(n) = y-m*x-c
             end
             meanx = sum(vx)/3;
             meany = sum(vy)/3;
             for i = 1:n
                 if double(subs(ineq(i),[x,y],[meanx,meany])) > 0
                     ineq(i) = -ineq(i);
                 end

             end
             obj = region(ineq,[x,y]);
         end

         function f = plus(obj1,obj2)
            l = []; 
            for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
            end
            %obj2 = obj2.removeDenominator;
            for i = 1:size(obj2.ineqs,2)
               lunique = true;
               for j = 1:size(obj1.ineqs,2)
                 if (obj1.ineqs(j) == obj2.ineqs(i))
                     lunique = false;
                     break;
                 end
                 if (obj1.ineqs(j).f == -obj2.ineqs(i).f)
                     f = region.empty;
                     return;
                 end
               end
               if lunique
                   l = [l,obj2.ineqs(i).f];
               end
            end 
            f = region(l,obj1.vars);
         end


         function v = getEndpoints (obj, i)
            % obj.print
             n = 0;
          %  disp('start')
             for j = 1:obj.nv
              %   double(subs(obj.ineqs(i).f,obj.vars,[obj.vx(j),obj.vy(j)]))
               if abs(double(subs(obj.ineqs(i).f,obj.vars,[obj.vx(j),obj.vy(j)]))) <= 1.0e-6
                   n = n + 1;
                   v(n,1) = obj.vx(j);
                   v(n,2) = obj.vy(j);
               end
             end
            
             v = sortrows(v);
         end


         % doesnt work when returning more than one region 
         % fix it
         % to detect - end point interval contained in bigger interval
         function f = minus(obj1,obj2)
             %disp('in minus')
             if (obj1 == obj2)
                 f = [region.empty];
                 return
             end     
            l = []; 
            for i = 1:size(obj1.ineqs,2)
                %disp('i')
                %i, size(obj1.ineqs,2)
                %ep = obj1.getEndpoints(i)
                v1(i,:,:) = obj1.getEndpoints(i);
                l = [l,obj1.ineqs(i).f];
            end
            l2 = [];
            rm = [];
            for i = 1:size(obj2.ineqs,2)
              v2(i,:,:) = obj2.getEndpoints(i);
              ladd = true;  
              for j = 1:size(obj1.ineqs,2)
                if (obj2.ineqs(i) == obj1.ineqs(j)) 
                    % replace with all
                    lv = true;
                    for i1 = 1:size(v1,2)
                        for i2 = 1:size(v1,3)
                            lv = lv & (v1(j,i1,i2) == v2(i,i1,i2));
                        end
                    end
                    if lv
                      rm = [rm,j];
                    end  
                    ladd = false;
                    break;
                end
              end
              if ladd
                  l2 = [l2,obj2.ineqs(i).f];
              end
            end
            l(rm) = [];
            for i = 1: size(l2,2)
                l = [l,-l2(i)];
            end
            if (size(l,2)) == 0
                f = [region.empty];
                return
            end

            % original code
%             for i = size(l,2)
%                 l(i)
%             end
            f0 = region(l,obj1.vars);
            if isempty(f0)
                f = [f0];
                return
            end
%            f0 = f0.getVertices;
%            disp('here')
            %f0.print
            %f0.nv
            if f0.nv >= 3
              %f = [f0.simplify(obj1.vars)];
              f = [f0.simplify];
              return
            else  
              f1 = f0.divideRegions(obj1);
              f = [];
              for i = 1: size(f1,2)
                  if f1(i).nv <= 2
                      continue
                  end
                  
                  %f0 = f1(i).simplify(obj1.vars);
                  %disp('divided region')
                  %f1(i).print
                  f0 = f1(i).simplify;
                  %disp('divided region after simplify')
                  %f0.vars
                  %f0.print
                  if f0 == obj1
                      continue
                  end
                  f = [f,f0];
              end
            end

            
            % create multiple regions here
            
            
         end
         
         function f = minus0(obj1,obj2)
%              if obj1.not | obj2.not
%                  disp("Implement not in minus")
%              end
            l = []; 
            for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
            end
            %n = size(obj1.ineqs,2);
            l2 = [];
            for i = 1:size(obj2.ineqs,2)
              ladd = true;  
              for j = 1:size(obj1.ineqs,2)
                if (obj2.ineqs(i) == obj1.ineqs(j))
                    ladd = false;
                    break;
                    
                end

                
              end
              if ladd
                  l2 = [l2,obj2.ineqs(i).f];
              end
            end
            if (size(l2,2) ~= 0)
              mult = (-1) ^ size(l2,2);
              for i = 1: size(l2,2)
                mult = mult * l2(i);
              end
              l = [l,mult];
            end
            f = region(l,obj1.vars);
         end
         
         
         function l = eqVertices(obj1,obj2)
             l = false;
             %obj1.nv
             if obj1.nv ~= obj2.nv
                return
             end
             for i = 1:obj1.nv
                 l1(i) = false;
                 l2(i) = false;
             end

             for i = 1:obj1.nv
               for j = 1:obj2.nv
                   if l2(j)
                       continue
                   end
                   if (abs(obj1.vx(i) - obj2.vx(j))<1.0d-15) & (abs(obj1.vy(i) - obj2.vy(j)) < 1.0d-15)
                       l1(i) = true;
                       l2(j) = true;
                       break;
                   end
               end
             end
             if (all(l1)==true )
                 l = true;
             else
                 l = false;
             end
         end

         function l = eq(obj1,obj2)
             l = false;
             if (size(obj1.ineqs,2)~=size(obj2.ineqs,2))
                 return;
             end
             if ~ obj1.eqVertices(obj2)
                 return
             end
             l = true;
         end

         function obj = unique(obj)
             n = 0;
             
            duplicates = [];
            for i = 1:size(obj.ineqs,2)
                for j = i+1:size(obj.ineqs,2)
                    if (obj.ineqs(i) == obj.ineqs(j))
                        n = n + 1;
                        duplicates(n) = j;
                        lMark(j) = 1;
                    end
                end

            end
           
            obj.ineqs(duplicates) = [];
         end

         function m = slope (obj,i,j)
          m = (obj.vy(i)-obj.vy(j))/(obj.vx(i)-obj.vx(j));
         end

         function m = slopeIneq(obj,i)
             if obj.ineqs(i).isLinear
                 c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                 m = -c(1);
             else
                 m = -intmax;
                 disp("to be implemented in slopeIneq")
             end
         end

          function c = yIntercept (obj,i,m)
          c = obj.vy(i)-m*obj.vx(i);   
          end

         function set = ithSet (obj, logical_indices)
           set = true;
           if all(logical_indices == false)
             set = false;
           end
                   
           for j = 1:size(logical_indices,2)
              if (logical_indices(j))
                set = set & obj.ineqs(j).f<=0;
              end
           end
            
         end
         
         function obj = removeSuperSet(obj, vars)
           % powerset indices
           %set_size=size(obj.ineqs,2);
           %logical_indices = mat2cell(double(dec2bin(0:2^set_size-1,set_size))==48,ones(2^set_size,1));
           logical_indices = supersetIndices(size(obj.ineqs,2))
           
           for i = 1:size(logical_indices,1)
              lS(i) = false  ;
           end
           disp('superset')
           size(logical_indices,1)
           %for i1 = 1:size(logical_indices,1)
             i1 = 1;  
             set1 = ithSet (obj, logical_indices{i1})

             s1 = solve(set1,vars)
             %if islogical(set1) 
             %    continue
             %end
             %if (~ set1) 
             %    continue
             %end
             %assume(set1);
               
             for i2 = 1:size(logical_indices,1)
                 if i1 == i2
                     continue
                 end
               set2 = ithSet (obj, logical_indices{i2})
               s2 = solve(set2,vars)
               if isequal(s1,s2)
                   disp('found')
               end
               continue;
               assume(set1);
               l2 = isAlways(set2)
               if ( l2)
                 assume(set2);
                 l1 = isAlways(set1)
               end
               
               continue
               
               if isAlways(set2)
                   lS(i2) = true;
                   disp('found')

                   logical_indices{i1}
                   logical_indices{i2}
              %     duplicates = [];
              %     for j = 1:size(logical_indices{i2},2)
              %         logical_indices{i1}(j)
              %         logical_indices{i2}(j)
              %         logical_indices{i1}(j)  ~= logical_indices{i2}(j)
              %         if (logical_indices{i1}(j)  ~= logical_indices{i2}(j) )
              %             duplicates = [duplicates,j]
              %         end
              %     end
              %     obj.ineqs(duplicates) = [];
                   
                   %return
               end
               %return
               
             %end
             %unassume(set1);

           end
           lS
         end

         function print(obj)
             disp("Variables")
             obj.vars
             disp(["nVertices = ", num2str(obj.nv)]);
             fprintf("vx =  ")
             fprintf("%d  ", obj.vx);
             fprintf("\n")
             fprintf("vy =  ")
             fprintf("%d  ", obj.vy);
             fprintf("\n\n")
%              disp("not")
%              obj.not
%              if obj.not
%                disp("Not of union of following ineqs")
%                obj.ineqs.printLIneq;
%              else
               disp("Intersection of following ineqs")
               obj.ineqs.printLIneq;
%              end
         %    disp(["Vertices = ", num2str(obj.nv)]);
         %    disp(obj.vx)
         %    disp(obj.vy)
         end

         function fprint(obj, uNo)

             %fprintf(uNo, obj.vars);
             fprintf(uNo, num2str(obj.nv)+"\n");
             for i = 1:obj.nv
               fprintf(uNo, num2str(obj.vx(i)) + "  " + num2str(obj.vy(i)) + "\n")  
             end
             %fprintf(uNo, num2str(obj.vx));
             %fprintf(uNo, num2str(obj.vy));
             obj.ineqs.fprintLIneq(uNo);
         end

         
         function plot (obj)
             
             l1 = min(min(obj.vx),min(obj.vy));
             if l1 < -10
                 l1 = -10;
             end
             l2 = max(max(obj.vx),max(obj.vy));
             if l2 > 10
                 l2 = 10;
             end
             l1 = -10;
             l2 = 10;
           obj.ineqs.plotLIneq (obj.vars, [l1,l2])   ;

         end

         function plotL (obj)
             figure;
             size(obj)
             for i = 1:size(obj,2)
                 i
                 obj(i).plot
                 obj(i).plotRegion;
             end
         end


       
         function value = getGlobalParameter(obj)
 %           persistent icolor;
 %           value = icolor
             persistent icolor;
            if isempty (icolor)
                icolor=1;
            else
                icolor=icolor+1;
            end
            value = icolor;
         end

         function [vx,vy] = plotRegionC (obj, textR, c)
            limitsx = [-6,6];
            limitsy = [-15,20];
            
            pts = 75;
            colors = ['b', 'r', 'g', 'm', 'c', 'y'];
            stepx = (limitsx(2)-limitsx(1))/pts;
            stepy = (limitsy(2)-limitsy(1))/pts;
            n = 0;
            vx = [];
            vy = [];
            ci = limitsx(1);
            for i = 1:pts
                cj = limitsy(1);
                for j = 1:pts
                  if obj.ptFeasible (obj.vars,[ci,cj]);
                      n = n+1;
                      vx(n) = ci;
                      vy(n) = cj;
                  end
                  cj = cj+stepy;    
                end
                ci = ci+stepx;    
            end
            %c = colors(1+mod(obj.getGlobalParameter,6))
            if n == 0
                disp ('region not displayed')
            end
            avx = sum(vx)/n;
            avy = sum(vy)/n;
            m = 0;
            for i = 1:n
                if (abs(vx(i) -avx)<0.1 & abs(vy(i) -avy)<0.1)
                    continue
                end    
                m = m+1;
                vx1(m) = vx(i);
                vy1(m) = vy(i);
            end
            %text = "R";
            text(avx,avy,textR,'FontSize',12, 'FontWeight', 'bold','Color', 'k')
            if m==0
                fill(vx, vy, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            else
                fill(vx1, vy1, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            end
            
         end

         function [vx,vy] = plotRegion (obj, textR)
            limitsx = [-6,6];
            limitsy = [-15,20];
            
            pts = 75;
            colors = ['b', 'r', 'g', 'm', 'c', 'y'];
            stepx = (limitsx(2)-limitsx(1))/pts;
            stepy = (limitsy(2)-limitsy(1))/pts;
            n = 0;
            vx = [];
            vy = [];
            ci = limitsx(1);
            for i = 1:pts
                cj = limitsy(1);
                for j = 1:pts
                  
                  if obj.ptFeasible (obj.vars,[ci,cj]);
                      n = n+1;
                      vx(n) = ci;
                      vy(n) = cj;
                  end
                  cj = cj+stepy;    
                end
                ci = ci+stepx;    
            end
            c = colors(1+mod(obj.getGlobalParameter,6))
            if n == 0
                obj.print 
                obj.ptFeasible (obj.vars,[2,-9])
                disp ('region not displayed')
            end
            avx = sum(vx)/n;
            avy = sum(vy)/n;
            m = 0;
            for i = 1:n
                if (abs(vx(i) -avx)<0.1 & abs(vy(i) -avy)<0.1)
                    continue
                end    
                m = m+1;
                vx1(m) = vx(i);
                vy1(m) = vy(i);
            end
            %text = "R";
            text(avx,avy,textR,'FontSize',12, 'FontWeight', 'bold','Color', 'k')
            if m==0
                fill(vx, vy, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            else
                fill(vx1, vy1, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            end
            
         end

         function plotByVertex (obj)
             %figure;
             limit = 6;
             for i = 1:obj.nv
                 vertices_ineq1(i, 1) = obj.vx(i);
                 if (vertices_ineq1(i, 1) > limit) 
                     vertices_ineq1(i, 1) = limit
                 end
                 if (vertices_ineq1(i, 1) < -limit) 
                     vertices_ineq1(i, 1) = -limit
                 end

                 vertices_ineq1(i, 2) = obj.vy(i);

                 if (vertices_ineq1(i, 2) > limit) 
                     vertices_ineq1(i, 2) = limit
                 end
                 if (vertices_ineq1(i, 2) < -limit) 
                     vertices_ineq1(i, 2) = -limit
                 end
                 
             end
             fill(vertices_ineq1(:, 1), vertices_ineq1(:, 2), 'b', 'FaceAlpha', 0.5);
         end
        
         function m = slopes (obj)
           n = 0 ;  
           for i = 1:size(obj.ineqs,2)
              if obj.ineqs(i).isLinear
                  c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                  n = n + 1;
                  m(n) = -c(1)/c(2);
              end
           end
         end
         % max over a region
         function [l, fmax, index, lsing] = maxArray (obj, f1, f2) 
          lsing = false;  
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          %obj.print
          %disp("Slopes")
          m = obj.slopes();
          for i = 1:size(fv1,2)
              sv1(i) = fv1(i).f;
              sv2(i) = fv2(i).f;
              %if (isnan(sv1(i)))
              %    sv1(i) = 0;
              %    sv2(i) = 0;
              %elseif (isnan(sv2(i)))
              %    sv1(i) = 0;
              %    sv2(i) = 0;
              %end
          end
          l = true;
          %sv1
          %obj.print;
          %double(sv1);
          %double(sv2);
          if all(abs(double(sv1 - sv2))< 1.0d-14)
              %disp("FIX IN MAXIMUM")
              %size(sv1)
              if size(sv1) == 1
                  %disp("Singleton")
                  fmax = f1;
                  index=1;
                %  return
              end
              nx = 0;
              for iv = 1:obj.nv
                  if (abs(obj.vx(iv))==intmax|abs(obj.vy(iv))==intmax)
                      continue
                  end
                  nx = nx + 1;
                  vx0(nx) = obj.vx(iv);
                  vy0(nx) = obj.vy(iv);
              end
              
              if nx > 1
                sx = mean(vx0);
                sy = mean(vy0);
                [l,fmax,index] = maxFromPt(obj, [sx,sy], [f1,f2]);
                if l 
                    return
                end
              end
              %disp("Here")
              % use slope mid pt to get directions
              for i = 1:size(m,2)
                
                for j = i+1: size(m,2)
                  if (abs(m(i))~= inf) & (abs(m(j))~= inf)
                  d =  (m(i)+m(j) )/2;
                  else
                      disp('infinity')
                  if (abs(m(i))==inf)
                        d = tan((pi/2 + atan(m(j)))/2);
                    else
                        d = tan((pi/2 + atan(m(i)))/2);
                  end
                  end
                  c = vy0(1) - d * vx0(1);
                  px = vx0(1) + 0.1;
                  py = d*px+c;
                  %double(px),double(py)
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                  %l
                  if l 
                    return
                  end
                  px = vx0(1) - 0.1;
                  py = d*px+c;
                  
                  %double(px),double(py)
                %  obj.ptFeasible(obj.vars, [px,py])
                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
                 % l
                  if l 
                    return
                  end
                  
                end
              end

                  
%                  px = vx0(1) + 0.1
%                  py = vy0(1)
%                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
%                  if l 
%                    return
%                  end
%                  px = vx0(1) - 0.1
%                  py = vy0(1)
%                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
%                  if l 
%                    return
%                  end
%                  px = vx0(1) 
%                  py = vy0(1)+ 0.1
%                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
%                  if l 
%                    return
%                  end
%                  px = vx0(1) 
%                  py = vy0(1) - 0.1
%                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
%                  if l 
%                    return
%                  end
%                  px = vx0(1) + 0.1
%                  py = vy0(1) + 0.1
%                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
%                  if l 
%                    return
%                  end
%                  px = vx0(1) - 0.1
%                  py = vy0(1) - 0.1
%                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
%                  if l 
%                    return
%                  end
%                  px = vx0(1) - 2/9
%                  py = vy0(1) - 1/9
%                  [l,fmax,index] = maxFromPt(obj, [px,py], [f1,f2]);
%                  if l 
%                    return
%                  end
%               
              disp("SINGLETON REGION")
              lsing = true;
              
              end
              

          if all(double(sv1) <= double(sv2))
              fmax = f2;
              index=2;
          elseif all(double(sv2) <= double(sv1))
              fmax = f1;
              index=1;
          else
              l = false;
              fmax = 0;
              index=0;
          end
         
        end
        
        function [l,fmax,index] = maxFromPt(obj, s, f)
          l = false;
          fmax = f(1);
          index=1;
          if obj.ptFeasible(obj.vars,s)
            fv01 = f(1).subsF(obj.vars,s).f;
            fv02 = f(2).subsF(obj.vars,s).f;
            
          else
              return;
          end
          l = true;
          %abs(double(fv01 - fv02))
          if abs(double(fv01 - fv02))> 1.0d-14
            if double(fv01) < double(fv02)
              fmax = f(2);
              index=1;
            else
              fmax = f(1);
              index=2;
            end
            return
          end
          l = false;
          fmax = f(1);
          index=1;
        end

        function [r] = splitmax2 (obj, f1, f2) 
           % disp('splitmax2')
          %f1
          %f2
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          f = f1-f2;
          vars = f.getVars;
          x = sym('x');
          y = sym('y');
          f = subs(f.f, vars ,[x,y]);
          
          fx = [];
          fy = [];
          n = 0;
          fxy = [];
          for i = 1:size(obj.ineqs,2)
            g = subs(obj.ineqs(i).f, vars,[x,y]);
            s = solve ([f==0,g==0],[x,y]);
            if isempty(s)
              continue;
            end
            
            if isempty(s.x)
              continue;
            end
            sx = double(s.x);
            sy = double(s.y);
            if ~obj.ptFeasible (vars,[sx,sy])
                continue
            end
            fx = [fx,sx];
            fy = [fy,sy];
            n = n+1;
            fxy(n,1)= sx;
            fxy(n,2)= sy;
          end
          %fxy
          fxy = unique(fxy,"rows");
          fx = fxy(:,1);
          fy = fxy(:,2);
          %f1
          %f2
          %size(fx)
          if size(fx,1) == 2
              m = (fy(2)-fy(1))/(fx(2)-fx(1));
              c = fy(1) - m*fx(1);
              ineq = vars(2)-m*vars(1)-c;
          end
          for i = 1:size(fv1,2)
              if (double(fv1(i).f) > double(fv2(i).f))
                  if subs(ineq,vars,[obj.vx(i),obj.vy(i)]) <= 0
                    r = [ineq,-ineq];
                    return
                  else
                    r = [-ineq,ineq];
                    return  
                  end
              end

          end
          
        end

        function [r] = splitmax3 (obj, f1, f2) 
           % disp('splitmax2')

          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          f = f1-f2;
           vars = f.getVars;
%           x = sym('x');
%           y = sym('y');
%           f = subs(f.f, vars ,[x,y])
%           
%           fx = [];
%           fy = [];
%           n = 0;
%           fxy = [];
%           for i = 1:size(obj.ineqs,2)
%             g = subs(obj.ineqs(i).f, vars,[x,y]);
%             s = solve ([f==0,g==0],[x,y]);
%             if isempty(s)
%               continue;
%             end
%             
%             if isempty(s.x)
%               continue;
%             end
%             sx = double(s.x);
%             sy = double(s.y);
%             if ~obj.ptFeasible (vars,[sx,sy])
%                 continue
%             end
%             fx = [fx,sx];
%             fy = [fy,sy];
%             n = n+1;
%             fxy(n,1)= sx;
%             fxy(n,2)= sy;
%           end
%           %fxy
%           fxy = unique(fxy,"rows");
%           fx = fxy(:,1);
%           fy = fxy(:,2);
%           %f1
%           %f2
%           %size(fx)
%           if size(fx,1) == 2
%               m = (fy(2)-fy(1))/(fx(2)-fx(1));
%               c = fy(1) - m*fx(1);
%               ineq = vars(2)-m*vars(1)-c;
%           end
          ineq = f;
          vars
          f.print
          for i = 1:size(fv1,2)
              if (double(fv1(i).f) >= double(fv2(i).f))
                  if subs(ineq.f,vars,[obj.vx(i),obj.vy(i)]) <= 0
                    r = [ineq,-ineq];
                    return
                  else
                    r = [-ineq,ineq];
                    return  
                  end
              end

          end
          
        end

        
        function [nl, v1, v2] = splitmaxArray (obj, f1, f2) 
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          v1 = [];
          v2 = [];
          for i = 1:size(fv1,2)
              sv1(i) = fv1(i).f
              sv2(i) = fv2(i).f
              % replace with limit
              if (isnan(sv1(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              elseif (isnan(sv2(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              end
              if sv1(i) == sv2(i)
                  v2 = [v2,i];
                  v1 = [v1,i];
              elseif sv1(i) > sv2(i)
                  v1 = [v1,i];
              else    
                  v2 = [v2,i];
              end

          end
          
          if size(v1,2) == 3
              nl(1) = 1;
          else
              nl(1) = 0;
          end
          if size(v2) == 3
              nl(2) = 1;
          else
              nl(2) = 0;
          end
          
        end
        
        function [lelim, obj] = deleteIneq (obj, vars)
          lelim = false;

          for i = 1:size(obj.ineqs,2)
              %obj.ineqs(i).f
              l(i) = obj.ineqs(i).f;
              mark(i) = false;
              %obj.ineqs(i).print
              nPts = 0;

              for j = 1:obj.nv
                  if obj.ineqs(i).subsF(obj.vars,[obj.vx(j),obj.vy(j)]).isZero
                      nPts = nPts+1;
                  end
              end
              if nPts <= 1
                  mark(i)=true;
              end
              %i,nPts
          end
         % l

          for i = 1:size(obj.ineqs,2)
              if ~mark(i)
                  continue
              end
             l1 = l;
             l1(i) = [];
          %   l1
             r = region (l1,vars);
           %  disp('in delete')
           %  r.print
             if isempty(r)
                 return
             end
%              if r.nv == 0
%                  return
%              end
             if obj.eqVertices(r)
            %     disp('equal')
                 lelim = true;
                 obj = r;
                 return;
             end
          end


        end


        function obj = simplify (obj) %, vars)
            %disp('in simplify')
            %obj.print
          for i = 1:size(obj.ineqs,2)
            [lelim, obj] = obj.deleteIneq (obj.vars);
            %i
            %lelim
            %disp("in simplify")
            %r0.print
            %obj = r0;
            if ~lelim
                return
            end
          end

        end

         % removes redundant ineq by finding points of intersection
         % check this code again
         function obj = simplify0 (obj, vars, obj2)
           %  obj.print
           rem = [];
           n = 0;
           x = sym('x');
           y = sym('y');
           %intersectingPts = []
           %intersectingEdges = []
           
           for i = 1:size(obj.ineqs,2)
             %  obj.ineqs(i).f
               
             f1 = subs(obj.ineqs(i).f, vars,[x,y]);
             for j = i+1:size(obj.ineqs,2)
                 f2 = subs(obj.ineqs(j).f, vars,[x,y])  ;
                 s = solve ([f1==0,f2==0],[x,y]);
                 
                 %disp('s.x')
                 %f1
                 %f2
                 %s.x
                 %s.y
                 if isempty(s.x)
                     continue;
                 end
                 sx = double(s.x);
                 sy = double(s.y);
                 l =true;
                 if nargin > 2
                     l = obj2.ptFeasible (vars,[sx,sy]);
                 end
                 %if obj2.ptFeasible (vars,[s.x,s.y])
                 if l
                     pointExists = false;
                     for k = 1:n
                         
                         if (intersectingPts(k,1:2) == [sx,sy])
                             pointExists = true;
                             break;
                         end
                     end
                     if pointExists
                       intersectingEdges(k,nEdges(k)+1) = i;
                       intersectingEdges(k,nEdges(k)+2) = j;
                       nEdges(k) = nEdges(k) + 2;
                     else
                       n = n + 1;
                       intersectingPts(n,1:2) = [sx,sy];
                       intersectingEdges(n,1:2) = [i,j];
                       nEdges(n) = 2;
                     end
                 end
                 
             end 
           end
           intersectingPts
           intersectingEdges
           %nEdges
           for i = 1:size(intersectingEdges,1)
               lp(i) = false;
           end
           
           % Mark feasible points as true
           for i = 1:size(intersectingEdges,1)
               %disp('intersecting point')
               %intersectingPts{i}
               %obj.ptFeasible (vars,intersectingPts{i})
               if obj.ptFeasible (vars,intersectingPts(i,1:2))
                 lp(i) = true;
                 continue;
               end
           end
           lp
           for i = 1:size(obj.ineqs,2)
               ls(i) = false;
           end
           % Mark pair of edges for each feasible point
           % split into 2 loops - pair and more

           for i = 1:size(obj.ineqs,2)
             ledgesP(i) = false;
           end 
           for i = 1:size(intersectingEdges,1)
               if ~lp(i)
                   continue
               end
               if nEdges(i) == 2
                   ls(intersectingEdges(i,1:2)) = true; 
                   nEdgesP(i,1:2) = intersectingEdges(i,1:2);
                   ledgesP(intersectingEdges(i,1))=true;
                   ledgesP(intersectingEdges(i,2))=true;
               end
           end
           ls
           for i = 1:size(intersectingEdges,1)
               if ~lp(i)
                   continue
               end
               if nEdges(i) == 2
                   continue
               end
               %edges = [];
               for j = 1:nEdges(i) %:2
                   if (mod(j,2)==0)
                       continue
                   end
                   if (~        ls(intersectingEdges(i,j))|ls(intersectingEdges(i,j+1))  )  
                      %edges = [edges,[intersectingEdges(i,j),intersectingEdges(i,j+1)]];
                      ls(intersectingEdges(i,j:j+1)) = true; 
                      %break
                  end
               end
           end
           ls
           % nEdges(i) : length of intersectingEdges(i,:)
           % At this point all feasible vertices are covered but some have
           % more than 2 edges hence next few loops to remove extra edges
           for i = 1:size(obj.ineqs,2) %size(intersectingEdges,1)
             
             ledges(i) = false;
           end 
           for i = 1:size(intersectingEdges,1)
             nEdgesM(i) = 0;
             
           end 
           %ledges
           % nEdgesM : Number of edges thro feasible vertex
           % logical ledges used to not double count
           for i = 1:size(intersectingEdges,1)
               if ~lp(i)
                   continue
               end
               for j = 1:nEdges(i)
                   if ~ledges(intersectingEdges(i,j))
                   if (ls(intersectingEdges(i,j)))
                       ledges(intersectingEdges(i,j)) = true;
                       nEdgesM(i) = nEdgesM(i)+1;
                   end
                   end
                   
               end    
               for j = 1:nEdges(i)
                 ledges(intersectingEdges(i,j)) = false;
               end
           end
           nEdgesM
           ledges
           % putting neccesary edges in nEdgeP to decide which ones can be
           % removed
           for i = 1:size(intersectingEdges,1)
               if ~lp(i)
                   continue
               end
               if (nEdgesM(i) == 2 & nEdges(i) == 2)
                   
                   continue
               elseif (nEdgesM(i) == 2)
                   
                 n = 0;  
                 
                 for j = 1:nEdges(i)
                   if ledgesP(intersectingEdges(i,j))  
                       continue
                   end
                   if (ls(intersectingEdges(i,j)))
                       n = n+1;
                       nEdgesP(i,n) = intersectingEdges(i,j);
                   end
                 end    
               end
               
           end
           nEdgesP
%unmark array is used to remove some edges which were earlier added
           for i = 1:size(obj.ineqs,2)
               lsp(i) = false;
               unmark(i) = false;
           end

           for i = 1:size(intersectingEdges,1)
               unMark(i) = false;    
               if nEdgesM(i) ~= 2
                   continue
               end
               if nEdgesP(i,1) == 0   % edge was marked due to other vertex hence not recorded here
                   continue;
               end
               lsp(nEdgesP(i,1)) = true;
               lsp(nEdgesP(i,2)) = true;
           end 
           lsp
       %    size(lsp)
           
           for i = 1:size(intersectingEdges,1)
               if nEdgesM(i) <= 2
                   continue
               end
               n = 0;
               for j=1:nEdges(i)
                  if (lsp(intersectingEdges(i,j))) 
                      n = n + 1;
                  else
                      unmark(intersectingEdges(i,j)) = true;
                  end
               end
               if (n >= 2)
               for j=1:nEdges(i)
                  if (lsp(intersectingEdges(i,j))) 
               %       n = n + 1
                  else
                      unmark(intersectingEdges(i,j)) = true;
                  end
               end
               
               end
           end
           unmark
           for i = 1:size(obj.ineqs,2)
               if (ls(i) & ~lsp(i) & unmark(i))
                   ls(i) = false;
               end
           end

           ls

           rem = [];
           for i = 1:size(ls,2)
               
              if (~ls(i))
                rem = [rem,i];
              end
           end
           rem
           obj.ineqs(rem) = [];
           return
         end

         function obj = simplify2 (obj, vars)
           rem = [];
             for i = 1:size(obj.ineqs,2)
             for j = i+1:size(obj.ineqs,2)
                 s = solve ([obj.ineqs(i).f==0,obj.ineqs(j).f==0],vars)
                 if isempty(s.x)
                     continue;
                 end
                 if obj2.ptFeasible (vars,[s.x,s.y])
                 
                 %obj.ptFeasible (vars,[s.x,s.y])
                 if obj.ptFeasible (vars,[s.x,s.y])
                     continue;
                 end
                 rem = [rem,i];
                 rem = [rem,j];
                 end
                 
             end 
           end 
           
           obj.ineqs(rem) = [];
         end

         function l = ptFeasible(obj, vars, point)
           l = true;
           for i = 1:size(obj.ineqs,2)
               %subs ([obj.ineqs(i).f],vars,point)
               for j = 1:size(point,1)
%                  obj.ineqs(i).f
%                  point(j,:)
%                 double(subs ([obj.ineqs(i).f],vars,double(point(j,:))))
                if double(subs ([obj.ineqs(i).f],vars,double(point(j,:)))) > 1.0e-12
                   l = false;
                   return
               end
               end
           end  
         end


         function obj = removeDenominator(obj)
           for i = 1:size(obj.ineqs,2)
              
              obj.ineqs(i) = obj.ineqs(i).removeDenominator2;
           end 

         end

         function l = isFeasibleWBPts(obj)
             l = false;
             for i = 1:obj.nv
               if ~ obj.ptFeasible(obj.vars, [obj.vx(i),obj.vy(i)])
                   return
               end
             end
             l = true;
         end

         function l = isFeasible(obj)
             l = false;
             for i = 1:size(obj.ineqs,2)
                 %obj.ineqs(i).print
               for j = i+1:size(obj.ineqs,2)
                   %obj.ineqs(j).print
                 if ( obj.ineqs(i) == unaryminus(obj.ineqs(j)))
                 %if ( obj.ineqs(i) == -obj.ineqs(j))
                     %disp('negative')
                     %obj.ineqs(i).f
                     %obj.ineqs(j).f
                    
                     return
                 end
                 if ( obj.ineqs(i) == obj.ineqs(j))
                     continue
                 end
                 
                 s = solve(obj.ineqs(i).f<=0,obj.ineqs(j).f<=0);
                 if isempty(s)
                     
                     %obj.ineqs(i).f
                     %obj.ineqs(j).f
                     
                     return;
                 end
                 % put case with only 1 variable then activate this
               %  class(s)
               %  if ~s.isSymType
               %  if isempty(s.x)
               %      return;
               %  end
               %  if isempty(s.y)
               %      return;
               %  end
               %  end
                 %obj.ineqs(i).print
                 %obj.ineqs(j).print
               end 
 
             end 
             
             l = true;
         end

         % stupid way of doing this
         function obj = intersection (obj1, obj2)
%          if obj1.not | obj2.not
%                  disp("Implement not in intersection")
%              end
         l = [];
           for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
           end 
         for i = 1:size(obj2.ineqs,2)
             lf = true;
             for j = 1:size(obj1.ineqs,2)
               if obj1.ineqs(j)==obj2.ineqs(i)
                   lf = false;
                   break;
               end
             end
             if lf
               l = [l,obj2.ineqs(i).f];
             end
         end 
         obj = region(l, obj1.vars);
         if isempty(obj)
             return
         end
         %obj.print
         if (isFeasible(obj))
             %disp('feasible')
             obj = obj.unique;
             if obj.nv <= 2
                 %disp('degenerate ');
                 obj = region.empty;
             else
                 if obj.nv == 0
           disp('in intersection nv 0')
         end
         
             end
         else
             obj = region.empty;
         end

     end

     function [linter, objR] = intersection2(obj1, obj2, lprint)
         % split into linear and quad first
         %disp("Intersection2")

         
         linter = true;
         
         objR = obj1;
         obj = [obj1,obj2];
         nl = 0;
         nq = 0;
         vars = obj1.vars;
         
         for j = 1:size(obj,2)
           for i = 1:size(obj(j).ineqs,2)  
             if obj(j).ineqs(i).isLinear
                nl = nl+1; 
                l(nl) = obj(j).ineqs(i);
             else
                nq = nq+1; 
                q(nq) = obj(j).ineqs(i);
             end
           end
            
            
         end
      if lprint
          disp('nl nq')
          nl
          nq
      end

         if nl > 0
             
             %[notF,l] = l(1:nl).removeParallel(vars, lprint);
             %if notF
             %    linter = false  ;
                 %disp("Not feasible 1")
           
             %    return
             %end
             
             %if lprint
             %disp('removeParallel')
             %    l.printL
             %end

             % 21 DEC commenting next few lines to check 
%              [notF,l] = l.removeSum(lprint);
%              if notF
%                  linter = false  ;
%                  %disp("Not feasible 2")
%            
%                  return
%              end
%              if lprint
%              disp('removeSum')
%              l.printL
%              end
             l2 =[];
             for i = 1:size(l,2)
                 l2 = [l2,l(i).f];
             end
             lR = region(l2, obj1.vars);
             if lprint
                 disp('vertex simplify')
             lR = lR.getVertices();
             lR.isFeasibleWBPts
             lR.print;
             lR.plot;
             end
             %size(lR.ineqs,2)
             lR = lR.simplify;  %(vars);
             %disp('aft simplify')
             %lR.print;
             
             % this can be removed - check
             % check for empty region instead - basically clean up
             if size(lR.ineqs,2) == 0
                 linter = false  ;
                % disp("Infeasible 3")
                 return
             end
             %disp("Linear")
             %lR.print
             
         end
         %nl
         %nq
         if nq > 0
         %disp("Quadratic")
         
          %   q
             q2 =[];
             for i = 1:size(q,2)
                 q2 = [q2,q(i).f];
             end
             
             qR = region(q2,obj1.vars);
           %  disp("lR")
           %  lR.print
            % disp("qR")
           %  qR.print;
           %  lR.print;
             objR = lR + qR;
         else    
             objR = lR;
         end
         %disp("Not in intersection2")
         %objR.not
         objR = objR.getVertices();
         if lprint
             disp("objR b4 simplify")
             
             objR.print
             disp("objR")
             
             objR = objR.simplify; % (objR.vars)
             objR.print
         end
         
         if objR.nv == 0
           linter=false  ;
         else
           linter = isFeasible(objR);
         end
     end
     
     function [linter,obj] = intersection3(obj1, obj2, lprint)
         % split into linear and quad first
         %disp("Intersection3")
         obj = region.empty();
       %  obj1.ineqs.printL
       %  obj2.ineqs.printL
         linter = false;
         n = 0;
         [linter, inter] = intersection2(obj1, obj2,lprint);
         if linter
             n = n + 1  ;
             obj(n) = inter;
         else
             obj(1) = obj1;     
         end  
     end
     
     function obj = normalizeEdge(obj)
         for i = 1:size(obj.ineqs,2)  
             c = getLinearCoeffs (obj.ineqs(i),obj.vars);
             
             if (c(2) == 0)
                 % double * f not overloaded
                 obj.ineqs(i).f = (1/c(1)) * obj.ineqs(i).f;
             else
                 obj.ineqs(i).f = (1/c(2)) * obj.ineqs(i).f;
             end
         end
     end

     function obj = normalize1(obj)
         for i = 1:size(obj.ineqs,2)
             obj.ineqs(i) = obj.ineqs(i).normalize1;
         end
     end
     % wont work for degree > 2
     % removed complex vertices
     function obj = getVertices(obj)
         
       obj.nv=0;
       obj.vx=[];
       obj.vy=[];
       t1 = sym('t1');
       t2 = sym('t2');
       varsTemp = [t1,t2];
       for i = 1:size(obj.ineqs,2)  
           f1 = obj.ineqs(i);
           f1 = f1.subsF (obj.vars,varsTemp);
           for j = i+1:size(obj.ineqs,2)  
               f2 = obj.ineqs(j);
               f2 = f2.subsF (obj.vars,varsTemp);
               s = solve ([f1.f==0,f2.f==0],varsTemp);
               
               if isempty(s)
                   continue;
               elseif isempty(s.t1)
                   continue;
               elseif isempty(s.t2)
                   continue;
               end
               if ~ isreal(s.t1)
                   s.t1
                   continue;
               end
               
               if (obj.ptFeasible(obj.vars, [s.t1,s.t2]))
                   for k = 1:size(s.t1,1)
                   obj.nv=obj.nv+1;
                   obj.vx(obj.nv) = double(s.t1(k));
                   obj.vy(obj.nv) = double(s.t2(k));
                   end 
               end
               
           end
           
       end

         
       %[obj.vx,obj.vy] = poly2cw(obj.vx,obj.vy);

       % putting intmax for inf to avoid Nans 
       % intmax + intmax = intmax
       %disp('ptFeasile')
       %obj.print
       n = obj.nv;
       if obj.ptFeasible(obj.vars, [intmax,intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = intmax;
          obj.vy(obj.nv) = intmax;
          %disp("++")
       end
       if obj.ptFeasible(obj.vars, [intmax,-intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = intmax;
          obj.vy(obj.nv) = -intmax;
          %disp("+-")
       end
       if obj.ptFeasible(obj.vars, [-intmax,intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = -intmax;
          obj.vy(obj.nv) = intmax;
          %disp("-+")
       end
       if obj.ptFeasible(obj.vars, [-intmax,-intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = -intmax;
          obj.vy(obj.nv) = -intmax;
          %disp("--")
       end
       for i = 1:n
           if nargin == 2
               %intmax,obj.vy(i)
               obj.ptFeasible(obj.vars, [intmax,obj.vy(i)])
           end
           if obj.ptFeasible(obj.vars, [obj.vx(i),intmax])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = obj.vx(i);
             obj.vy(obj.nv) = intmax;
           end
           if obj.ptFeasible(obj.vars, [obj.vx(i),-intmax])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = obj.vx(i);
             obj.vy(obj.nv) = -intmax;
           end
           if obj.ptFeasible(obj.vars, [intmax,obj.vy(i)])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = intmax;
             obj.vy(obj.nv) = obj.vy(i);
           end
           if obj.ptFeasible(obj.vars, [-intmax,obj.vy(i)])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = -intmax;
             obj.vy(obj.nv) = obj.vy(i);
           end
           
       end
       if obj.nv == 0
           return;
       end
       for i = 1:obj.nv
           V(i,1) = obj.vx(i);
           V(i,2) = obj.vy(i);
       end 
       %disp('unique vertices')
       %obj.vx,obj.vy
       V = unique(V,"rows");
       if obj.nv ~= size(V,1)
         obj.nv = size(V,1);
        % obj.vx
        % obj.vy
         obj.vx = V(:,1)';
         obj.vy = V(:,2)';
       end

     end



      function [nv, vx, vy] = vertexOfEdge(obj,ind)
              nv = 0;
%              disp('in')
%              ind
%              obj.nv
             % obj.vx
%              disp('out')
%              size(obj.vx)
              vx = [];
              vy = [];
               for i = 1:obj.nv
                 f1 = obj.ineqs(ind).subsF (obj.vars,[obj.vx(i),obj.vy(i)]);
                 if f1.isZero
                     nv = nv + 1;
                     vx(nv) = obj.vx(i);
                     vy(nv) = obj.vy(i);
                 end
               end
          end
      function f = divideRegions(obj,obj2)
         
       obj.nv=0;
         t1 = sym('t1');
         t2 = sym('t2');
         varsTemp = [t1,t2];
         ir = 0;
         % in this loop we are getting all points of intersection of pair
         % of ineqs.
         % remove duplicate points 
         for i = 1:size(obj.ineqs,2)  
           f1 = obj.ineqs(i);
           f1 = f1.subsF (obj.vars,varsTemp);
           for j = i+1:size(obj.ineqs,2)  
               f2 = obj.ineqs(j);
               f2 = f2.subsF (obj.vars,varsTemp);
               s = solve ([f1.f==0,f2.f==0],varsTemp);
               
               if isempty(s)
                   continue;
               elseif isempty(s.t1)
                   continue;
               elseif isempty(s.t2)
                   continue;
               end
               lpt = false;
               for k  = 1:ir
                   if all(points(k,:) == [s.t1,s.t2])
                       lpt = true;
                       break;
                   end
               end
               if lpt
                   continue
               end
               if ~obj2.ptFeasible (obj2.vars,[s.t1,s.t2])
                 continue
               end
               
               ir = ir + 1;
               index(ir) = 0;
               points(ir,1:2) = [s.t1,s.t2]; 
               for k = 1:size(obj.ineqs,2)  
                 if double(subs ([obj.ineqs(k).f],obj.vars,[s.t1,s.t2])) < 1.0e-12
                     index(ir) = index(ir)+1;
                     r(ir,index(ir)) = obj.ineqs(k).f;
                 end
               end
               
           end
           
         end
         % at this point r(ir) corresponds to satisfying ineqs for each
         % possible vertex

          %ir
          %index
          %for i = 1:ir
          %    points(i,:)
          %end
%           for i = 1:ir
%               disp('list')
%               for j = 1:index(ir)
%                   r(i,j)
%               end
%           end
          is = 0;
         s = r;
         for i=1:ir
                 leq = false;
                 for j = 1:is
                     if indexS(j) ~= index(i)
                         continue;
                     end
                     leq = true;
                     for k = 1:index(i)
                         if ~ isequal(r(i,k),s(j,k))
                             leq = false;
                             break;
                         end
                     end
                     if leq
                       break;
                     end
                 end
                 if leq
                    continue;
                 end
                 is = is+1;
             
                 indexS(is) = index(i);
                 for j = 1:indexS(is)
                     s(is,j) = r(i,j);
                 end
             
         end

         f = [];
         for i = 1:is
           f = [f, region(s(i,:),obj.vars)]; 
         end
       
       end

     function obj = conjugate(obj)
     end

     % If function is linear - check with lin ineqs 
                             % if coef are same try to resolve  
     function [l,less, fIneq, nf] = funcIneq (obj, f)
         % fix sign
         less = true;
       if ~f.isLinear
         fIneq = 0;
         l = false;
         nf = true
         return
       end
       c0 = f.getLinearCoeffs (obj.vars)
       for i =1:size(obj.ineqs,2)
           if ~obj.ineqs(i).isLinear
               continue;
           end
           c1 = obj.ineqs(i).getLinearCoeffs (obj.vars)
           if sign(c1(1)) ~= sign(c0(1))
               less = false;
               c1 = -c1;
           end
           if (c1(1) ~= c0(1))
               c1 = (c0(1)/c1(1)) * c1;
           end
           if (c1(2) == 0 & c0(2) == 0)
               fIneq = -c1(3);
               l = true;
               nf = false
               return
           elseif (c1(1)/c0(1) == c1(2)/c0(2))
               fIneq = -double(c1(3)/c1(2))+double(c0(3)/c0(2));
               l = true;
               nf = false
               return
           end
       end
       fIneq = 0;
       l = true;
       nf = true
         
     end
     
     % return function values at vertices
     function fv = funcVertices (obj, f)
         n = 0;
         for i =1:obj.nv
             if (abs(obj.vx(i))==intmax)
                 continue
             end
             if (abs(obj.vy(i))==intmax)
                 continue
             end
             n = n + 1;
             fv(n) = f.subsF(obj.vars,[obj.vx(i),obj.vy(i)]);
             if (isnan(fv(n).f))
                 fv(n) = f.limit(obj.vars,[obj.vx(i),obj.vy(i)]);
             end
         end
     end


     function [l, maxf, index, lSing] = maximum(obj, f)
          %obj.print
          %f(1).print
          %f(2).print
          
          if f(1)==f(2)
              l = true;
              lSing = false;
              maxf = f(1);
              index = 1;
              return
          end
          [l,  maxf, index, lSing] = obj.maxArray (f(1), f(2)) ;
     end

     function [f2, fe, d] = splitMax(obj, f, expr)
       [nl, v1, v2] = obj.splitmaxArray (f(1), f(2)) 
       f2 = [];
       fe = [];
       d = [];
       v{1}=v1;
       v{2}=v2;
       for i = 1:2
           if nl(i) == 0
               continue
           end
           f2 = [f2,f(i)];
           fe = [fe,expr(i)];
           d0 = obj.regionWPts(obj.vx(v{i}), obj.vy(v{i}), obj.vars(1), obj.vars(2));
           d = [d,d0];
           disp('Regiopn')
           d0.print
% create region with points

       end
       
     end
     
     function e = getOtherEdgeAtVertex (obj, ind, vertex)
         e = 0;
         for i = 1:size(obj.ineqs,2)
            if i == ind
                continue
            end
            f1 = obj.ineqs(i).subsF (obj.vars,vertex);
            if f1.isZero
                
                e = i;
                break
            end

         end
     end

     % wont work cause of intersection vs union
     function [l,obj] = merge (obj, obj2)
         l = false;
         n = 0; 
%          lQuad = false;
%          for i =1:size(obj.ineqs,2)
%               if obj.ineqs(i).isQuad
%                   lQuad = true;
%                   %return
%               end
%          end
%          for i =1:size(obj2.ineqs,2)
%               if obj2.ineqs(i).isQuad
%                   lQuad = true;
%                   %return
%               end
%          end
          
         % to be tested and added
         lQuad1 = false;
         nmq1 = 0;
         for i =1:size(obj.ineqs,2)
             if obj.ineqs(i).isQuad
                 lQuad1 = true
                 nmq1 = nmq1 + 1;
                 mq1(nmq1) = i;
             end
         end
         lQuad2 = false;
         nmq2 = 0;
         
         for i =1:size(obj2.ineqs,2)
             if obj2.ineqs(i).isQuad
                 lQuad2 = true
                 nmq2 = nmq2 + 1;
                 mq2(nmq2) = i;
               
             end
         end
         lQuad1
         lQuad2
         if (lQuad1 & lQuad2)
             disp("Quad merge")
             obj.print
             obj2.print
             for i = 1:nmq1
               for j = 1:nmq2
                   if obj.ineqs(mq1(i)) == -obj.ineqs(mq2(j))
                       n = n + 1;
                       marki(n) = mq1(i);
                       markj(n) = mq2(j);
                   end
               end
             end
             l = true;
             obj.ineqs(marki) = []; 
             obj2.ineqs(markj) = [];
             obj = obj+obj2;
             return
%          elseif lQuad1
%              return
%          elseif lQuad2
%              return
          end
         lQuad = lQuad1 | lQuad2;
         % cant do + as it returns empty due to contadictory ineqs
         %obj = obj + obj2;
         % hence remove those then do +
         %obj.print;
       %  disp('in merge')
       %  obj.print
       %  obj2.print
         %obj = obj.unique;
         %disp("unique")
         %obj.print;
         marki = [];
         markj = [];
         for i =1:size(obj.ineqs,2)
           for j =1:size(obj2.ineqs,2)
               %j
             if obj.ineqs(i) == -obj2.ineqs(j)
                 [nvi, vxi, vyi] = obj.vertexOfEdge(i);
                 [nvj, vxj, vyj] = obj2.vertexOfEdge(j);
              %   obj.ineqs(i)
              %   nvi,nvj
                 if nvi ~= nvj
                     continue
                 end
%                  if nvi ~= 2
%                      continue
%                  end
                  if nvi == 1 & (~lQuad)
                     
                       l = true;
                       n = n + 1;
                       marki(n) = i;
                       markj(n) = j;
                       break;
                     
                 end
                 if nvi == 2
                 
                 [tvxi, sortedIndices] = sort(abs(vxi));
                 vxi = vxi(sortedIndices);
                 vyi = vyi(sortedIndices);
                 
                 [tvxj, sortedIndices] = sort(abs(vxj));
                 vxj = vxj(sortedIndices);
                 vyj = vyj(sortedIndices);
                 
                 % check same vertices and convex angles 
                 if all(vxi == vxj) & all(vyi == vyj)
                    
                     edgeiNo = obj.getOtherEdgeAtVertex (i,[vxi(1),vyi(1)]);
                     if edgeiNo == 0
                         continue
                     end
                     if ~obj.ineqs(edgeiNo).isLinear 
                         continue;
                     end
                     edgejNo = obj2.getOtherEdgeAtVertex (j,[vxi(1),vyi(1)]);
                     if edgejNo == 0
                         continue
                     end
                     
                     if ~obj2.ineqs(edgejNo).isLinear 
                         continue;
                     end
                     
                     mi = obj.slopeIneq(edgeiNo);
                     a1 = atan(mi);
                     %if a1 < 0
                     %    a1 = a1 + 2*pi;
                     %end
                     mj = obj2.slopeIneq(edgejNo);
                     a2 = atan(mi);
                     %if a2 < 0
                     %    a2 = a2 + 2*pi;
                     %end
                     
                     if (a1+a2 > pi) | (a1+a2 < -pi)
                         continue
                     end

%%%%%%%%%%%%%%%%
                     
                    if nvi == 2   
                     
                     edgeiNo = obj.getOtherEdgeAtVertex (i,[vxi(2),vyi(2)]);
                     if edgeiNo == 0
                         continue
                     end
                     
                     if ~obj.ineqs(edgeiNo).isLinear 
                         continue;
                     end
                     edgejNo = obj2.getOtherEdgeAtVertex (j,[vxi(2),vyi(2)]);
                     if edgejNo == 0
                         continue
                     end
                     
                     if ~obj2.ineqs(edgejNo).isLinear 
                         continue;
                     end
                     
                     mi = obj.slopeIneq(edgeiNo);
                     a1 = atan(mi);
                     %if a1 < 0
                     %    a1 = a1 + 2*pi;
                     %end
                     mj = obj2.slopeIneq(edgejNo);
                     a2 = atan(mi);
                     %if a2 < 0
                     %    a2 = a2 + 2*pi;
                     %end
                     
                     if (a1+a2 > pi) | (a1+a2 < -pi)
                         continue
                     end
                    end
                 
%%%%%%%%%%%%%%%%%

                         %atan(mi)+atan(mj)
                       l = true;
                       n = n + 1;
                       marki(n) = i;
                       markj(n) = j;
                       break;
                     end
                 end  
             end
           end
         end
         if l
          % obj.print;
          % obj2.print;
          % marki
          % markj
           obj.ineqs(marki) = []; 
           obj2.ineqs(markj) = [];
           obj = obj+obj2;
           %obj.print
         end
         
         end
         
      

     function l = hasNegativeIneqs(obj)
         l = true;
         for i = 1:size(obj.ineqs,2)
           for j = i+1:size(obj.ineqs,2)
               if (obj.ineqs(i) == -obj.ineqs(j))
                   return;
               end
           end
         end
         l = false;
     end

     function [nmaxf,nmaxd] = mergeL(maxd,maxf)  % (obj,maxf,maxd)
          %disp('in merge')
          ia(1) = 1;
          n = 0;
          %size(maxf,2)
          for i = 1:size(maxf,2)
              %maxf(i).print
              marked(i) = false;
          end

          % ja has indices of all equal functions , ia by col no
          for i = 1:size(maxf,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(maxf,2)
                  
                  if isAlways(maxf(i).f == maxf(j).f)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
          %nmaxf = [];
          %nmaxd = [];
          m = 0;
          for i = 1:size(maxf,2)
              marked(i) = false;
          end
          ia
%          ja
        %  return
        %  nmaxf= [];
        %  nmaxd= [];
        for i = 1:size(maxf,2)
            
            if  marked(i)
                continue
            end
            if (ia(i) == ia(i+1)) 
                m = m + 1;
                nmaxf(m) = maxf(i);
                nmaxd(m) = maxd(i);
            else
                % get common boundary and merge
                % make groups and add 
               r = maxd(i);
               lmerge = true;
               while lmerge
                 lmerge = false;
                 for j=ia(i):ia(i+1)-1
                   
                   if marked(ja(j))
                       continue
                   end
                   [l,r] = r.merge (maxd(ja(j)));
                   
                   if l
                     marked(ja(j)) = true;
                     lmerge = true;
                   end
                 end
               end
               % fix here to combine others - currently combining only with
               % first
                 for j=ia(i):ia(i+1)-1
                   
                   if marked(ja(j))
                       continue
                   end
                   marked(ja(j)) = true;
                   m = m + 1;
                   nmaxf(m) = maxf(i);
                   nmaxd(m) = maxd(ja(j));  
                 end
               
               m = m + 1;
               nmaxf(m) = maxf(i);
               %r = r.getVertices();
               nmaxd(m) = r;  
                   
            end
        end
       % disp("m")
       % m
      end

     end

     
     
end