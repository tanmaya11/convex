x=0:0.1:10;
y1=exp(-x/2);
y2=exp(-x/3);
y3=4;
figure
hold all
plot(x,y1)
plot(x,y2)
%plot(x,y3)
patch([x fliplr(x) ], [y1 fliplr(y2) ], 'g')
hold off