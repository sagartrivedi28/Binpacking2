function Plotdata(Fstruct, Dstruct,ax1,ax2)

Rshift = Fstruct.Rshift;
Ddimx = Fstruct.Ddimx;
Ddimy = Fstruct.Ddimy;
Sline = Fstruct.Sline;
Rmax = Fstruct.Rmax;
Fcor = Fstruct.Fcor;
Fsets = Fstruct.Fsets;

axes(ax1);
cla;
axis(1.1*[-Rmax(1)/2,Rmax(1)/2+1,-Rmax(2)/2-1,Rmax(2)/2+1]);
rectangle('Position',[-Rmax(1)/2,-Rmax(2)/2,Rmax(1),Rmax(2)],'LineWidth',2,'LineStyle','--','EdgeColor','r');
rectangle('Position',[-Rshift(1),-Rshift(2),2*Rshift(1),2*Rshift(2)],'LineWidth',2,'LineStyle','--');
SetsColor = hot(length(Fsets));

for l=1:length(Fsets)
    for i=1:length(Fsets{l})
        rectangle('Position',[Fcor(Fsets{l}(i),1)+(Sline(1)/2)-Rshift(1), Fcor(Fsets{l}(i),2)+(Sline(2)/2)-Rshift(2),Ddimx(Fsets{l}(i))-Sline(1),Ddimy(Fsets{l}(i))-Sline(2)],'LineWidth',1,'FaceColor','c');
        rectangle('Position',[Fcor(Fsets{l}(i),1)-Rshift(1),Fcor(Fsets{l}(i),2)-Rshift(2),Ddimx(Fsets{l}(i)),Ddimy(Fsets{l}(i))],'LineWidth',1,'LineStyle','--','EdgeColor','b');
        text(Fcor(Fsets{l}(i),1)-Rshift(1)+Ddimx(Fsets{l}(i))/2,Fcor(Fsets{l}(i),2)-Rshift(2)+Ddimy(Fsets{l}(i))/2,num2str(Fsets{l}(i)),'FontSize',8,'BackgroundColor',SetsColor(l,:));
    end;
end;
text(0,Rshift(2)+1000,['Detected Minima Cut-Sets: ' num2str(length(Fsets))],'FontSize',12);





xoffset = Dstruct.xoffset;
yoffset = Dstruct.yoffset;

x_in = Dstruct.x_in;
y_in = Dstruct.y_in;
x_on = Dstruct.x_on;
y_on = Dstruct.y_on;
r=Dstruct.r;
rs=Dstruct.rs;

axes(ax2);
cla;
plot(x_in,y_in,'b','Linewidth',2);
hold on;
plot(x_on,y_on,'r--','Linewidth',1);
theta=0:2*pi/360:2*pi;
xc=(r-rs)*cos(theta);
yc=(r-rs)*sin(theta);
plot(xc,yc,'k');   
plot(xoffset,yoffset,'r+');

end
