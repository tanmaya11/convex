figure
ezplot(['x + 7*y - 10' ]) 
hold on
ezplot(['x + 2*y - 4' ]) 
hold on
ezplot(['48*x - 56*y + 4*x*y + x^2 + 4*y^2 - 184']) 
hold on
ezplot(['x - (9*y)/5 - 5' ]) 
hold on
ezplot(['0*x - y - 5 ']) 

return


figure
ezplot(['x + 2*y + 4'] ,[-25,5],[-6,11])   
hold on
ezplot(['x - (9*y)/5 - 5'] ,[-25,5],[-6,11])   
hold on
ezplot(['0*x -  y - 5'] ,[-25,5],[-6,11])   
hold on
ezplot(['x + 7*y - 46'] ,[-25,5],[-6,11])   
hold on


return

figure
ezplot('x + 7*y - 10' ,[-10,5],[-3,3])   
hold on
ezplot('x - (5*y)/7 - 25/7' ,[-10,5],[-3,3])   
hold on
ezplot('- x - 2*y - 4' ,[-10,5],[-3,3])   
hold on
ezplot('48*x - 56*y + 4*x*y + x^2 + 4*y^2 - 184' ,[-10,5],[-3,3])   
hold on
ezplot('x - (9*y)/5 - 5' ,[-10,5],[-3,3])   
return




figure
ezplot('-x-7*y-4',[-110,10],[-2,16])   
hold on
ezplot('56*y-48*x-4*x*y-x^2-4*y^2+184',[-110,10],[-2,16])   
hold on
ezplot('x+2*y-4',[-110,10],[-2,16])   
hold on
%ezplot('x - y/3 - 14/3',[-55,-20],[-6,11])   
hold on
return


figure
ezplot('-x-7*y+10',[-50,100],[-20,5])   
hold on
ezplot('-x-2*y-4',[-50,100],[-20,5])   
hold on
ezplot('-56*y+48*x+4*x*y+x^2+4*y^2-184',[-20,5],[-20,16])   
hold on
return
ezplot('0*x+y-2',[-110,10],[-2,16])   
hold on
%ezplot('x - y/3 - 14/3',[-55,-20],[-6,11])   
hold on
return

figure
ezplot('x+7*y-10',[-110,10],[-2,16])   
hold on
ezplot('5*y/7-x+25/7',[-110,10],[-2,16])   
hold on
ezplot('-x-2*y+4',[-110,10],[-2,16])   
hold on
%ezplot('x - y/3 - 14/3',[-55,-20],[-6,11])   
hold on
return







figure
ezplot('x + 2*y + 4',[-55,-20],[-6,11])   
hold on
%ezplot('x - (9*y)/5 - 5',[-55,-20],[-6,11])   
hold on
ezplot('0*x - y - 5',[-55,-20],[-6,11])   
hold on
%ezplot('x - y/3 - 14/3',[-55,-20],[-6,11])   
hold on
ezplot('x - 2*y + 44',[-55,-20],[-6,11])   
hold on
return

figure
ezplot('x + 7*y+ 4',[0,10],[-5,5])   
hold on
ezplot('x + 7*y - 10',[0,10],[-5,5])   
hold on
ezplot('148*x - 196*y + (x + 7*y)^2 - 684',[0,10],[-5,5])    
hold on
ezplot('x - 9/5*y - 5',[0,10],[-5,5])    
hold on
ezplot('x - y/3 - 14/3',[0,10],[-5,5])   
hold on
return

figure
ezplot('x + 7*y - 10',[0,10],[-5,5])   
hold on
ezplot('148*x - 196*y + (x + 7*y)^2 - 684',[0,10],[-5,5])    
hold on
ezplot('4 - 2*y - x ',[0,10],[-5,5])   
hold on
ezplot('x - y/3 - 14/3',[0,10],[-5,5])   
hold on
ezplot('(5*y)/7 - x + 25/7 ',[0,10],[-5,5])   
hold on
return

figure
%ezplot('-1+5/2*y+x',[-30,-20],[-5,5]) 
%hold on
ezplot('-x-7*y-4',[0,10],[-5,5]) 
hold on
ezplot('x+7*y-10',[0,10],[-5,5]) 
hold on
ezplot('4-2*y-x',[0,10],[-5,5]) 
hold on
ezplot('196*y - 148*x - (x + 7*y)^2 + 684',[0,10],[-5,5]) 
hold on

return

figure
ezplot('-1+5/2*y+x',[-30,-20],[5,15]) 
hold on
ezplot('-x-2*y-4',[-30,-20],[5,15]) 
hold on
ezplot('x-2*y+44',[-30,-20],[5,15]) 
hold on

return



ezsurf('x*y')
   xlim([-2, 2]); % Set x-axis range from 0 to 6
   ylim([-2, 2]); % Set y-axis range from 0 to 6
return

figure
ezplot('x + 5/2*y -1',[-30,30]) 
hold on
ezplot('-x - 2*y - 4',[-30,30]) 
hold on
ezplot('x + 7*y -46',[-30,30]) 
hold on
%ezplot('8*x + 24*y -(x+2*y)^2-32') 
%hold on
xlim([-100, 100]); % Set x-axis range from 0 to 6
ylim([-100, 100]); % Set y-axis range from 0 to 6
return


figure
ezplot('x+7*y+4') 
hold on
ezplot('x+7*y-10') 
hold on
ezplot('-x + (9*y)/5 + 5') 
hold on
ezplot('-148*x + 196*y - (x + 7*y)^2 + 684') 
hold on

return




%ezplot('-x+9*y/5+5
hold on
ezplot('0*x+y+5')  
hold on
ezplot('x+0*y+4')  
hold on

return




figure
ezplot('y-sqrt(x)')
hold on
ezplot('y-x/2')

hold on
return
figure
ezplot('0*x+y')  
hold on
ezplot('x+0*y')  
hold on
ezplot('x^2-y-1')  
hold on
return



figure
%ezplot('0*x+y-2') 
%hold on
ezplot('x+y-6') 
hold on
ezplot('0*x+y')  
hold on
ezplot('y^2-x')  
hold on
return
h = ezplot('-x+4*y-4')  
set(h, 'Color', 'g')
hold on
h = ezplot('-x+8*y-12')  
set(h, 'Color', 'r')
hold on
h = ezplot('x+8*y-20')  
set(h, 'Color', 'r')
hold on
return

ezplot('-5*y -4*x -20 ')
%ezplot('−4*x − 5*y − 20')
hold on
return

figure
ezplot('(5*x^2)/56 + (5*x)/7 - (5*y^2)/4 + 10/7')
hold on
ezplot('-x -7*y -4')
hold on
ezplot('x + 7*y -10')
hold on
ezplot('-x -2*y -4')
hold on
ezplot('x + 2*y -4')
hold on
ezplot('148*x- 196*y+(x+7*y)^2-684')
hold on
return

% ezplot('0*x-y-5')  
 hold on
% ezplot('x+7*y+4')  
%hold on
ezplot('-x+0*y-4')  
hold on
ezplot('9*y/5 -x + 5 ')
%ezplot('x+5*y/9 + 20/9 ')
return



ezplot('x - y/2 - 2')
hold on
ezplot('x - 2 + 0*y')
hold on
ezplot('(2*y)/5 - x + 8/5') 
xlim([-0.5, 2.5]); % Set x-axis range from 0 to 6
ylim([-4.5, 1.5]); % Set y-axis range from 0 to 6
%return



figure
ezplot('x + y/3 - 2')
hold on
ezplot('0*x -y')
hold on
ezplot('(2*y)/5 - x + 8/5') 
xlim([-0.5, 2.5]); % Set x-axis range from 0 to 6
ylim([-4.5, 1.5]); % Set y-axis range from 0 to 6
%return

figure
ezplot('x - y/2 - 2')
hold on
ezplot('0*x -y')
hold on
ezplot('(2*y)/5 - x + 8/5') 
xlim([-0.5, 2.5]); % Set x-axis range from 0 to 6
ylim([-4.5, 1.5]); % Set y-axis range from 0 to 6

%return

figure
ezplot('2 - y/3 - x')
hold on
ezplot('x - 2 + 0*y')
hold on
ezplot('(2*y)/5 - x + 8/5') 

xlim([-0.5, 2.5]); % Set x-axis range from 0 to 6
ylim([-4.5, 1.5]); % Set y-axis range from 0 to 6

return


% 
% figure
% y/7 - x + 4/7 <= 0 
% x + y/3 - 2 <= 0 
% -y <= 0 
% Variables
% 
% ans =
% 
% [x, y]
% 
%     "nVertices = "    "4"
% 
% vx =  5.714286e-01  1  1.600000e+00  1.818182e+00  
% vy =  0  3  0  5.454545e-01  
% 
% Intersection of following ineqs
% y/7 - x + 4/7 <= 0 
% x - (2*y)/5 - 8/5 <= 0 
% -y <= 0 
% x + y/3 - 2 <= 0 
% Variables
% 
% ans =
% 
% [x, y]
% 
%     "nVertices = "    "3"
% 
% vx =  2  1.818182e+00  1.600000e+00  
% vy =  0  5.454545e-01  0  
% 
% Intersection of following ineqs
% x + y/3 - 2 <= 0 
% -y <= 0 
% (2*y)/5 - x + 8/5 <= 0 
% return





figure;
ezplot(' x + y   ')
hold on
ezplot(' x - y   ')
hold on
ezplot(' 0.5 * x - y   ')
hold on
ezplot(' 0*x + y -1  ')
hold on
xlim([-2, 2]); % Set x-axis range from 0 to 6
ylim([-1, 2]); % Set y-axis range from 0 to 6

return


ezplot('x + 2*y + 4 ') 
hold on
ezplot('(9*y)/5 - x + 5  ')
hold on
ezplot('0*x + y + 5  ') 
hold on


figure;
ezplot('- x - 7*y - 4  ')
hold on
ezplot('x + 2*y + 4 ') 
hold on
ezplot('(9*y)/5 - x + 5  ')
hold on
ezplot('0*x + y + 5  ') 
hold on

return

figure;
ezplot('2*x^3-24*x^2+90*x+7 ')
return

figure;
ezplot('- x - 7*y - 4  ')
hold on
ezplot('x +7*y - 10  ')
hold on
ezplot('x + 2*y + 4 ') 
hold on
return


% - s1 - 7*s2 - 4 <= 0 
% s1 + 7*s2 - 10 <= 0 
% s1 + 2*s2 + 4 <= 0 
% 
% s1 + 7*s2 + 4 <= 0 
% s1 - (9*s2)/5 - 5 <= 0 
% - s2 - 5 <= 0 
% s1 + 2*s2 + 4 <= 0 

figure;
ezplot('- x - 7*y - 4  ')
hold on
ezplot('(9*y)/5 - x + 5  ')
hold on
ezplot('56*y-48*x-4*x*y-x^2-4*y^2+184') 
hold on
ezplot('x + 2*y - 4 ') 
hold on
ezplot('x + 2*y + 4 ') 
hold on
return




figure;
ezplot('- x - 7*y - 4  ')
hold on
ezplot('148*x - 196*y + (x + 7*y)^2 - 684 ')
hold on
ezplot('(9*y)/5 - x + 5  ')
hold on
ezplot('56*y-48*x-4*x*y-x^2-4*y^2+184') 
hold on
ezplot('x + 2*y - 4 ') 
hold on
return

figure;
ezplot('- x - 7*y - 4  ')
hold on
ezplot('x + 2*y + 4 ') 
hold on
ezplot('56*y-48*x-4*x*y-x^2-4*y^2+184') 
hold on
ezplot('x + 2*y - 4 ') 
hold on

return

figure;
ezplot('- x - 7*y - 4  ')
hold on
ezplot('(9*y)/5 - x + 5  ')
hold on
ezplot('x + 2*y + 4 ') 
hold on
ezplot('x + 2*y - 4 ') 
hold on
return


figure;
ezplot('- x - 7*y - 4  ')
hold on
ezplot('148*x - 196*y + (x + 7*y)^2 - 684 ')
hold on
ezplot('(9*y)/5 - x + 5  ')
hold on
X = [3.1591,3.1591,3.0591,3.2591,3.0591,3.2591];
Y = [-1.1227,-0.9227,-1.0227,-1.0227,-0.5381,-1.5073];
plot(X,Y);
hold on;

return
figure;
X = [-5,-5,1,0,-5,-5,1];
Y = [5,-4,3,-4,-4,5,3];

plot(X,Y);
hold on;
xlim([-6, 4]); % Set x-axis range from 0 to 6
ylim([-5, 6]); % Set y-axis range from 0 to 6

return
X = [1,0,2,2,1];
Y = [3,-4,0,1,3];
plot(X,Y);
hold on;

return
%d(1)=domain([-5,-4;0,-4;1,3;-5,5],x,y); 
 %           d(2)=domain([0,-4;2,0;2,1;1,3],x,y); 

figure;

hold on

ezplot('0*x + y + 5  ') 
hold on


ezplot('9*y/5-x+5')
hold on
return


figure;
ezplot('x + 0*y + 4 ') 
hold on
f = ezplot('- x/sqrt(2) - y - 7.824 ') 
set(f, 'Color', 'r');

hold on

ezplot('0*x + y + 5  ') 
hold on


ezplot('9*y/5-x+5')
hold on
return



ezplot('x + y + 2')

%ezplot('x + 2*y + 4 <= 0') 
hold on
ezplot('x + 0*y + 1 ') 
hold on
ezplot('0*x + y + 1 ') 
hold on
ezplot('2* x + y + 2')

%ezplot('x + 2*y + 4 <= 0') 
hold on
%ezplot('x + 7*y - 46 <= 0') 
%hold on
return
figure;
ezplot('1-5/2*y-x  ')
hold on
ezplot('- x - 2*y - 4 ')
hold on
ezplot('x - 2*y +44  ')
hold on
return





ezplot('(9*y)/5 - x + 5 ')
hold on
return

x = sym('x')
y = sym('y')
subs(4 - 2*y - x,[x,y],[-20,13])
subs(x-2*y+44,[x,y],[-20,13])
return

figure;
ezplot('8*x+24*y-(x+2*y)^2-32 ')
hold on
ezplot('x + (5*y)/2 - 1 ')
hold on
return



ezplot('4 - 2*y - x ')
 hold on
 ezplot('-x -2* y - 4 ') 
 hold on
ezplot('2*y-x-15 ')
hold on



return

%X = [-5,0,-5,-5];
%Y = [-4,-4,5,-4];
% X = [-5,0,1,-5];
% Y = [5,-4,3,5];
 X = [-5,-5,0,1,-5];
 Y = [5,-4,-4,3,5];

plot(X,Y);
xlim([-6, 2]); % Set x-axis range from 0 to 6
ylim([-5, 6]); % Set y-axis range from 0 to 6
return
ezplot('(9*y)/5 - x + 5 ')
hold on
%ezplot('x + 7*y - 10 ')
%hold on

ezplot('- x - 7*y - 4 ')
hold on
return

ezplot('148*x - 196*y + (x + 7*y)^2 - 684 ')
hold on
ezplot('y/3 - x + 14/3 ')
hold on

s2/3 - s1 + 14/3
return



ezplot('- y - 4 +0*x ')
hold on
ezplot('x + (5*y)/9 + 20/9 ')
hold on
ezplot('x + 3*y - 10 ')
hold on
ezplot('x - (6*y)/7 + 11/7  ')
hold on
ezplot('x - y/7 - 4/7')
hold on
ezplot('- x - 5 +0*y')
hold on
%return
%ezplot('-x + y ')
hold on
return
ezplot('-x -3 + y ')
hold on
ezplot('-x^2 +1 +  y ')
hold on
ezplot('-x^2 -1 + y ')
hold on
return
ezplot('x + y ')
hold on
ezplot('(5*y)/4 - x + 25/4')
hold on
ezplot('- x - (3*y)/2 - 3/2')
hold on
ezplot('27*x - (51*y)/2 + (x + (3*y)/2)^2 - 591/4 ')
hold on
ezplot('(9*y)/2 - 3*x + (x + (3*y)/2)^2 + 9/4')
hold on
return


(9*y)/2 - 3*x + (x + (3*y)/2)^2 + 9/4 <= 0 


ezplot('x+y-1')
hold on
ezplot('-x-y')
hold on
ezplot('y/3 - x - 1/6')
hold on
ezplot('x - (5*y)/4 - 25/4')
hold on
ezplot('(x + y)^2 - 8*x ')

hold on
ezplot('7*x - 5*y - 25')
hold on

return
ezplot('x + y - 1  ')
hold on
ezplot('- x - y  ')

hold on
ezplot('x - y/3 + 1/6 ')
hold on
ezplot('y - x - 1  ')
hold on
ezplot('(x + y)^2 - 8*x  ')
hold on
ezplot('18*y - 18*x - (x + y)^2 + 99  ')


hold on
return



return
%ezplot('-y')
%ezplot("-y")
%hold on
ezplot('log(abs(x))/log(exp(1))')
return
ezplot('x^2*(ln(abs(y)+0.5))'  )
return
%ezplot('exp(-x^2)',[-1,1]  )
%ezplot('sin(x)' ,[-1,1] )
ezplot('abs(x*y)'  )

hold on

%ezplot('x - 1' ) 
%hold on
ezplot('x - y'  )

hold on
ezplot('(2*y)/3 - x - 1') 
hold on
ezplot('- x - (2*y)/3' )
hold on
ezplot('6^(1/2)/4 - (6^(1/2)*x)/2 - y + 7/8'  )
hold on
ezplot('y + (6^(1/2)*x)/2 - 3/8' )
hold on
return
ezplot("x - y")
hold on

ezplot("-3*x + 2*y -3")
hold on
ezplot("x + y")
hold on
ezplot("x + y-2")
hold on
ezplot("x + 3*y/2+3/2")
hold on
ezplot("x + 3*y/2-9/2")
hold on
return
ezplot("3*x + 2*y")
hold on
return
ezplot('s2/2 - s1/2 - (s1 + s2)^2/36 + 11/4')
hold on
ezplot('(27*s1)/44 - (51*s2)/88 + (s1 + (3*s2)/2)^2/44 - 591/176')
hold on
return

(27*s1)/44 - (51*s2)/88 + (s1 + (3*s2)/2)^2/44 - 591/176 
s2/2 - s1/2 - (s1 + s2)^2/36 + 11/4 
(51*s2)/88 - (27*s1)/44 - (s1 + (3*s2)/2)^2/44 + 591/176 <= 0

%ezplot('x-y-2')
%hold on
%ezplot('-x+y')
%hold on
ezplot('8*x-(x+y)^2')
hold on
ezplot("x-y+1")
%ezplot("y-1")
hold on
ezplot("y+x-2")
hold on
ezplot("y+x")
hold on