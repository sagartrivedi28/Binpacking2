function [Fsets,Fwafer] = MPW_DicePacking(Fcor,Ddimx,Ddimy,Dvol,Gpara,ax)


Rmax = Gpara.Rmax;
Sline = Gpara.Sline;
Rno = Gpara.Rno;
N = 1:length(Ddimx);

Cmat = Cgraph(Fcor,Dvol,N);
[Fsets, Fwafer] = CGpartition(Cmat,Dvol,N,Rno);

Rshift = [max(Fcor(:,3))-min(Fcor(:,1)) max(Fcor(:,4))-min(Fcor(:,2))]/2;
axes(ax);
cla;
axis(1.1*[-Rmax(1)/2,Rmax(1)/2+1,-Rmax(2)/2-1,Rmax(2)/2+1]);
rectangle('Position',[-Rmax(1)/2,-Rmax(2)/2,Rmax(1),Rmax(2)],'LineWidth',2,'LineStyle','--','EdgeColor','r');
rectangle('Position',[-Rshift(1),-Rshift(2),2*Rshift(1),2*Rshift(2)],'LineWidth',2,'LineStyle','--');
SetsColor = hot(length(Fsets));

for l=1:length(Fsets)
    for i=1:length(Fsets{l})
        rectangle('Position',[Fcor(Fsets{l}(i),1), Fcor(Fsets{l}(i),2)+(Sline(2)/2),Ddimx(Fsets{l}(i)),Ddimy(Fsets{l}(i))],'LineWidth',1,'FaceColor','c');
        rectangle('Position',[Fcor(Fsets{l}(i),1)-(Sline(1)/2),Fcor(Fsets{l}(i),2)-(Sline(1)/2),Ddimx(Fsets{l}(i))+(Sline(1)/2),Ddimy(Fsets{l}(i))+(Sline(1)/2)],'LineWidth',1,'LineStyle','--','EdgeColor','b');
        text(Fcor(Fsets{l}(i),1)+Ddimx(Fsets{l}(i))/2,Fcor(Fsets{l}(i),2)+Ddimy(Fsets{l}(i))/2,num2str(Fsets{l}(i)),'FontSize',8,'BackgroundColor',SetsColor(l,:));
    end;
end;
text(0,Rshift(2)+1000,['Detected Minima Cut-Sets: ' num2str(length(Fsets))],'FontSize',12);
end



% Conflict Graph Generation Algorithm: Need to consider 300nm adjacent cut
function Cmat = Cgraph(Dcor,Dvol,N)
Cmat = false(N);
for k1=1:N
    for k2=1:N
        if (Dvol(k2)~=0) && (Dvol(k1)~=0)
            if ((Dcor(k1,1) < Dcor(k2,1)) && (Dcor(k2,1) < Dcor(k1,3))) || ((Dcor(k1,1) < Dcor(k2,3)) && (Dcor(k2,3) < Dcor(k1,3))) || ...
               ((Dcor(k1,2) < Dcor(k2,2)) && (Dcor(k2,2) < Dcor(k1,4))) || ((Dcor(k1,2) < Dcor(k2,4)) && (Dcor(k2,4) < Dcor(k1,4)))
                 Cmat(k1,k2)=1;
            end;
        end;
    end;
end;
Cmat = Cmat | Cmat';
end

% Graph Partition Algorithm: Need to consider the SameCut Dice
function [Nsets, Wvector] = CGpartition(Cmat,Dvol,N,Rno)

Tvol = Dvol;
% Vunit = min(Tvol(Tvol~=0));

sets=false(size(Cmat(:,1)));
% Tvol(temp(1))=0;

while sum(Tvol)
    temp = find(Tvol==min(Tvol(Tvol~=0)));
    Sflag = false;
    
    for l=1:size(sets,2) 
            if ~((sum(Cmat(:,temp(1)).*sets(:,l)) || sets(temp(1),l)) || Sflag)
                sets(temp(1),l)=true;
                Sflag=true;
            end;
    end;
    
    if Sflag==false
        sets = [sets false(size(Cmat(:,1)))];
        sets(temp(1),end) = true;
    end;
    
    Tvol(temp(1))=0;
end;

Tvol = Dvol;

while sum(Tvol)
    temp = find(Tvol==max(Tvol));
    
    for l=1:size(sets,2) 
            if ~(sum(Cmat(:,temp(1)).*sets(:,l)) || sets(temp(1),l))
                sets(temp(1),l)=true;                
            end;
    end;  
    
    Tvol(temp(1))=0;
end;
sets(Dvol==0,end)=true;

Nsets={};
for k=1:size(sets,2)
    Nsets{k}=find(sets(:,k))';
end;

Wvector = Wminimize(sets,Dvol,Rno);

end

function Wvector = Wminimize(sets,Dvol,Rno)

Ksets = double(sets);
Bsets = Ksets(1,:);
NDvol = Dvol(1);
l=1;
for r = 2:length(Rno)
    if (Rno(r-1)==Rno(r))
        Bsets(l,:) = Bsets(l,:)+Ksets(r,:);
        NDvol(l) = NDvol(l)+Dvol(r);
    else
        Bsets = [Bsets; Ksets(r,:)];
        NDvol = [NDvol;Dvol(r)];
        l=l+1;
    end;
    
end;

WIvector = zeros(length(sets(1,:)),1);

fun = @(Wvector) sum(Bsets*ceil(Wvector)-NDvol); %fun = @(Wvector) sum(ceil(Wvector));
A=-Bsets;
b=-NDvol;
Aeq=[]; beq=[];
lb = WIvector; ub=Inf*ones(size(WIvector));
Wvector = fmincon(fun,WIvector,A,b,Aeq,beq,lb,ub);

end
