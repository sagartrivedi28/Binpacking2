function Fstruct = DicePackingR7( Did,Ddimx,Ddimy,Dvol,Dreq,Gpara)

IValue = 2000; % in Micrometer
Fstruct.Fcor=0;
Fstruct.Fwafer=0;

Rmax = Gpara.Rmax;
Sline = Gpara.Sline;
Rno = Gpara.Rno;
wsize = Gpara.wsize;

TBflagM = Gpara.TBflag;
TBflag1 = false;
ShiftFlag = false;

if TBflagM    
    ShiftFlag = true;
    CRblock = [Gpara.TBcorner(1)*Gpara.TBcorner(3), Gpara.TBcorner(2)*Gpara.TBcorner(4)];   % width, height, column, row
    CNblock = [Gpara.TBcenter(1)*Gpara.TBcenter(3), Gpara.TBcenter(2)*Gpara.TBcenter(4)];    % width, height, column, row
    
    MinR = max([2*CRblock+2*median(Ddimx),CNblock+2*median(Ddimx), 12000]); 
else
    MinR = 12000;
end;

BladeSplit = Dreq(:,1);
Caround = Dreq(:,2);
SameCut = Dreq(:,3);
Iso1mm = Dreq(:,4);
DiceGroup = Dreq(:,6);

Ddimx = Ddimx+IValue*(Iso1mm~=0);
Ddimy = Ddimy+IValue*(Iso1mm~=0);

%%%%%%%%%%%%% Super Die %%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NDdimy,ind1] = sort(Ddimy,'descend');
NDid = Did(ind1);
NDdimx = Ddimx(ind1);
NDvol = Dvol(ind1);
NDreq = Dreq(ind1,:);

N = length(NDid);
Wcnt = Inf;
NBladeSplit = BladeSplit(ind1);
NDiceGroup = DiceGroup(ind1);
% NIso1mm = Iso1mm(ind1);
% NCaround = Caround(ind1);
NSameCut = SameCut(ind1);

th_NDdimy=1.25;                                    %empty not enough ==> allow shorter y-size dice placement  
th_nRow=1.2;                                       %as previously placed dice's ysize > th_nRow*new_ysize ==> trigger multi-row mode 
[~,ind_mvol]=sort(NDvol,'ascend');          %sort by volume
[~,~,vol_grp]=unique(NDvol(ind_mvol));
ind_vol=zeros(size(ind_mvol));
for i=1:max(vol_grp)                                 %same volume sort by ydim 
ind_vol(vol_grp==i)=sort(ind_mvol(vol_grp==i),'descend');
end



for oldW = Rmax(1):-1000: MinR

i=1; 
T_height=0;
shelf.height = [];
shelf.empty = oldW;
shelf.Dinfo = [];
shelf.level=[];
shelf.nRow=1;       % for multiple row optimization/ Total # rows
shelf.Rid = [];         % Each dice row ID
W = oldW;
TBflag = TBflagM;
ShiftFlag = TBflagM;

height_nRow=0;
Dcor = zeros(N,4);
DPool=true(N,1);
newnRow=false;
SplitFlag = false;
ScutFlag = false;
ShiftUpFlag=false;

if max(NBladeSplit)==2
    DPool2 = DPool & (NBladeSplit==2);
    DPool = ~DPool2;
    SplitFlag=true;
end;

% if max(NDiceGroup)>0
%     
% end;

% if max(NSameCut)>0
%    
%     
% end;


while any(DPool)   
    
%%%%%%%%%%%%%%%%%% TB Placement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if (i<=4 && TBflag &&  shelf(i).empty(1)==W) % Row 1,2 for corner & 3 for center TB placement 
    if ShiftUpFlag
        ShiftUp=300-Sline(2);
        shelf(i-1).height(1) = shelf(i-1).height(1)+ShiftUp;
    else
        ShiftUp=0; 
    end
    shelf(i).level = T_height + ShiftUp;
    shelf(i).nRow = 1;         
    height_empty=0;
    ShiftUpFlag=false;
    if i<=2            
        T_height = T_height + CRblock(2) + ShiftUp;       
        shelf(i).height = CRblock(2);                            
        height_nRow=shelf(i).height;  
        W = oldW-CRblock(1);
        shelf(i).empty = W - CRblock(1);                
        
    elseif i==3        
        T_height = T_height + CNblock(2) + ShiftUp;         
        shelf(i).height = CNblock(2);                            
        height_nRow=shelf(i).height;        
        W = (oldW-CNblock(1))/2;
        shelf(i).empty = W;
                       
    else
        i=i-1;        
        W = oldW;
        shelf(i).empty = shelf(i).empty+(W-CNblock(1))/2;
        TBflag = false;
        TBflag1 = true;
    end;

end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TBflagM && i<=3
    DPoolT = DPool & (NDdimy< th_NDdimy*shelf(i).height);
else
    DPoolT = DPool;
end;


n=find(DPoolT,1,'first');
 if ~any(shelf(i).empty(1) >= NDdimx(n)) 
       n=find(shelf(i).empty(1) >= NDdimx & DPoolT & NDdimy>=height_nRow(1)/th_NDdimy,1,'first');      %  Finding Dices which are high enough & within empty Space.
 end
 
 
 n_multi=[]; 
 
 % IF available valid candidate for placement
 if ~isempty(n) 
     
     if (TBflagM && i<=3) && (shelf(i).height<NDdimy(n))
         T_height = T_height + (NDdimy(n)-shelf(i).height);
         shelf(i).height = NDdimy(n);         
         
     end;
     
      % Single Rows Initialization
     if (shelf(i).empty(1)==W) && (~TBflag)                                           
         
         TBflag1 = false;
         if ShiftUpFlag
             ShiftUp=300-Sline(2);
             shelf(i-1).height(1) = shelf(i-1).height(1)+ShiftUp;
         else
             ShiftUp=0; 
         end         
         shelf(i).height = NDdimy(n);
         shelf(i).level = T_height + ShiftUp;
         shelf(i).nRow = 1;
         height_nRow=shelf(i).height;
         height_empty=0;
         T_height = T_height+NDdimy(n) + ShiftUp;
         ShiftUpFlag=false;                     
         
     else
         nRow=shelf(i).height/NDdimy(n);
         if length(shelf(i).empty)==1 && nRow>=th_nRow        % Multiple Rows Initialization
               shelf(i).nRow=ceil(nRow);
               shelf(i).empty=shelf(i).empty*ones(ceil(nRow),1);
               newnRow=true;
         end
     end

     DPool(n) = false;
     [shelf,Dcor,Dloc(1)] = BulidShelf(shelf,Dcor,W,NDid,NDdimx,NDdimy,i,n,height_nRow,1); %Single Row
         
     shelf(i).Rid = [shelf(i).Rid;1];
     if newnRow
         height_nRow(1)=NDdimy(n); 
     end;
     
 end

 for nr=2:shelf(i).nRow % Multi-Row
    
          if newnRow
             height_empty(nr)=shelf(i).height-sum(height_nRow);
          else
             height_empty(nr)=height_nRow(nr)+(shelf(i).height-sum(height_nRow))*(nr==shelf(i).nRow);
          end
          x_empty=shelf(i).empty(nr)>=NDdimx(ind_vol);
          y_empty=height_empty(nr)>=NDdimy(ind_vol);
          n_multi = ind_vol(find(x_empty & y_empty & DPool(ind_vol),1,'first'));
          if ~isempty(n_multi)
               if newnRow
                   height_nRow(nr)=NDdimy(n_multi);
               end;
              DPool(n_multi) = false;
              [shelf,Dcor,Dloc(nr)] = BulidShelf(shelf,Dcor,W,NDid,NDdimx,NDdimy,i,n_multi,height_nRow,nr);
              shelf(i).Rid = [shelf(i).Rid; nr];
          end
    
 end
      shelf(i).nRow=length(height_nRow);
      newnRow=false;
            

      if isempty(n) && isempty(n_multi) && any(DPool)
          
          if ~TBflag1 && ~TBflag  % Fix the issue
              [Dcor,shelf]=xEqualSpace(Dcor,shelf,W,i);
          end;         
         
            if ~isempty(shelf(i).Dinfo)
                H_i=unique(roundoff(sort(Dcor(shelf(i).Dinfo(:,1),4),'descend')),'stable');        
                H_diff_i=H_i(1)-H_i(min(2,length(H_i)));
                ShiftUpFlag= 0 < H_diff_i && H_diff_i < 300-Sline(2);
             else
                 ShiftUpFlag = 0;
            end;
         
             i=i+1;
             shelf(i).empty=W;
             height_nRow=0;
      end
      
 %%%%%%%%%%%% Split Blade %%%%%%%%%%%%%
      if SplitFlag && ~any(DPool)
          DPool=DPool2;
          SplitFlag=false;
          if ~(isempty(n) && isempty(n_multi))
              i=i+1;
              shelf(i).empty=W;
              height_nRow=0;
          end;           
           T_height = T_height+1000;
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
end;

if ~TBflag1 && ~TBflag   % Fix the issue
    [Dcor,shelf] = xEqualSpace(Dcor,shelf,W,i);
end;

Dcor = EdgeAlignment(Dcor,shelf,W,Dvol,Sline); 

if ShiftFlag
    
    for roll=2:length(shelf)
        [shelf,Dcor] = SwapShelf(shelf,Dcor,roll, min([roll+1,length(shelf)]));        
    end;
    
    ShelfC = max([2,round(length(shelf)/2)]);
    for roll=2:ShelfC
        [shelf,Dcor] = SwapShelf(shelf,Dcor,roll,min([roll+1,ShelfC]));
    end;     
     
end;

%Dcor = ChipAround(Dcor,shelf,Caround);

if ~isempty(wsize)
Rnew = [max(Dcor(:,3))-min(Dcor(:,1)) max(Dcor(:,4))-min(Dcor(:,2))]/10^3;  
Dstruct = DieOffset(Rnew(1),Rnew(2),wsize);
GDPW = sum(Dstruct.map_in(:));
else
GDPW = min(NDvol(NDvol~=0));
end

if T_height< Rmax(2)
Cmat = Cgraph(Dcor,Dvol,N,Sline);
[sets, Wvector] = CGpartition(Cmat,Dvol,N,Rno);
if Wcnt>sum(ceil(Wvector/GDPW))
    Wcnt = sum(ceil(Wvector/GDPW));            
    Fwafer = Wvector;
    Fshelf = shelf;
    Fcor = Dcor;
    FCmat = Cmat;
    Fsets = sets;
    
    % figure; imshow(FCmat,[]); colormap jet;
    Ddimx = Ddimx-IValue*(Iso1mm~=0);
    Ddimy = Ddimy-IValue*(Iso1mm~=0);
    Rshift = [max(Fcor(:,3))-min(Fcor(:,1)) max(Fcor(:,4))-min(Fcor(:,2))]/2;
    Fstruct.Rshift = Rshift;
    Fstruct.Ddimx = Ddimx;
    Fstruct.Ddimy = Ddimy;
    Fstruct.Sline = Sline;
    Fstruct.Rmax = Rmax;
    Fstruct.Fcor = Fcor;
    Fstruct.Fwafer = Fwafer;
    Fstruct.Fsets = Fsets;

end;
end;
shelf=[];

end;

end



% Conflict Graph Generation Algorithm: Need to consider 300nm adjacent cut
function Cmat = Cgraph(Dcor,Dvol,N,Sline)
Cmat = false(N);
Sline_mat=ones(size(Dcor,1),1)*Sline([1,2,1,2]);
Dcor=roundoff(Dcor);
Dcor_m220=roundoff(Dcor-300+Sline_mat);
Dcor_p220=roundoff(Dcor+300-Sline_mat);

for k1=1:N
    for k2=1:N
        if (Dvol(k2)~=0) && (Dvol(k1)~=0)
            if (Dcor_m220(k1,1) < Dcor(k2,1) && Dcor(k2,1) < Dcor_p220(k1,3)...
                && Dcor(k2,1) ~= Dcor(k1,3) && Dcor(k1,1) ~= Dcor(k2,3)...
                && (Dcor(k2,1) ~= Dcor(k1,1) || Dcor(k1,3) ~= Dcor(k2,3)))||...
               (Dcor_m220(k1,1) < Dcor(k2,3) && Dcor(k2,3) < Dcor_p220(k1,3)... 
                && Dcor(k2,1) ~= Dcor(k1,3) && Dcor(k1,1) ~= Dcor(k2,3)...
                && (Dcor(k2,1) ~= Dcor(k1,1) || Dcor(k1,3) ~= Dcor(k2,3)))||...
               (Dcor_m220(k1,2) < Dcor(k2,2) && Dcor(k2,2) < Dcor_p220(k1,4)...
                && Dcor(k2,2) ~= Dcor(k1,4) && Dcor(k1,2) ~= Dcor(k2,4)...
                && (Dcor(k2,2) ~= Dcor(k1,2) || Dcor(k1,4) ~= Dcor(k2,4)))||...
               (Dcor_m220(k1,2) < Dcor(k2,4) && Dcor(k2,4) < Dcor_p220(k1,4)... 
                && Dcor(k2,2) ~= Dcor(k1,4) && Dcor(k1,2) ~= Dcor(k2,4)...
                && (Dcor(k2,2) ~= Dcor(k1,2) || Dcor(k1,4) ~= Dcor(k2,4)))
            
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

function [shelf,Dcor,Dloc] = BulidShelf(shelf,Dcor,W,NDid,NDdimx,NDdimy,i,n,H_nRow,nr)

H_nRow_Room=H_nRow(nr)-NDdimy(n)+(shelf(i).height-sum(H_nRow))*(nr==shelf(i).nRow);
if nr>length(shelf(i).level)
   shelf(i).level(nr)=shelf(i).level(nr-1)+H_nRow(nr-1);
end
Dloc = W-shelf(i).empty(nr);
shelf(i).empty(nr) = shelf(i).empty(nr)-NDdimx(n);
shelf(i).Dinfo = [shelf(i).Dinfo; NDid(n) NDdimx(n) Dloc H_nRow_Room];  % ID, X-dimention, X-location, Height Leftover
Dcor(NDid(n),:) = [Dloc,shelf(i).level(nr),Dloc+NDdimx(n),shelf(i).level(nr)+NDdimy(n)];

end

function [Dcor,shelf]=xEqualSpace(Dcor,shelf,W,i)
                
                Dcor_i = roundoff(Dcor(shelf(i).Dinfo(:,1),:));
               [~,~,level_grp]=unique(Dcor_i(:,2)); % find the xo of multi-row
                InterSect_bw=sum(Dcor_i(:,3)*ones(1,length(Dcor_i(:,1)))==ones(length(Dcor_i(:,3)),1)*(Dcor_i(:,1).'),2)>=2;
                W_shelf=W-min(shelf(i).empty);
                xo_nRow_shelf=Dcor_i(InterSect_bw,3); if isempty(xo_nRow_shelf);xo_nRow_shelf=0;end;
                
                PREnRow=Dcor_i(:,3)<=xo_nRow_shelf;
                level_grp(~PREnRow)=level_grp(~PREnRow)+1*any(PREnRow);

                empty=shelf(i).empty;
                xoffset=zeros(size(empty));
                
                if any(PREnRow) % shift in x from the xo of multi-row
                     xoffset=min(shelf(i).empty)*xo_nRow_shelf/W_shelf; 
                     empty=[xoffset;shelf(i).empty-xoffset];
                     shelf(i).empty=shelf(i).empty-xoffset;                                                 
                     xoffset=[0;xoffset*ones(length(empty)-1,1)];
                end
                
                for s=1:max(level_grp) %equally distribute the empty space in x for shelf i 
                    
                     sublevel=level_grp==s;       
                     Dcor_subset=Dcor(shelf(i).Dinfo(sublevel,1),:);
                     shelf_Dinfo_subset = shelf(i).Dinfo(sublevel,3);
                     [~,ind]=sort(Dcor_subset(:,1));
                     
                     [~,~,H_grp_subset] = unique(Dcor(shelf(i).Dinfo(sublevel,1),[2,4]),'rows'); %Treat adjacent same-height dice as a group
                     Hsame_bw_subset = [false;H_grp_subset(ind(2:end))==H_grp_subset(ind(1:end-1))];
                     delta_list=zeros(length(Hsame_bw_subset),1);
                     for h=2:length(Hsame_bw_subset)
                         if Hsame_bw_subset(h)
                            delta_list(h)=max(delta_list);
                         else
                            delta_list(h)=max(delta_list)+1; 
                         end
                     end
                     
                     if max(delta_list)~=0
                           delta_list=delta_list*empty(s)/(max(delta_list)+1);
                     else
                           delta_list(:)=0;
                     end

                     Dcor_subset(ind,[1,3])=Dcor_subset(ind,[1,3])+delta_list*ones(1,2)+xoffset(s);
                     shelf_Dinfo_subset(ind)=shelf_Dinfo_subset(ind)+delta_list+xoffset(s);
                     
                     Dcor(shelf(i).Dinfo(sublevel,1),:)=Dcor_subset;
                     shelf(i).Dinfo(sublevel,3)=shelf_Dinfo_subset;
                     
                end

end

function Dcor = EdgeAlignment(Dcor,shelf,W,Dvol,Sline)
    
    xshift=(Dvol==0);
    for i=1:length(shelf)
          level(i) = shelf(i).level(1);
    end
    level(end+1)=inf;
    
    while any(~xshift)
       n=find(~xshift,1,'first'); xshift(n)=true;
       ind=find(Dcor(n,2)>=level,1,'last');
       SameShelf = level(ind) <= Dcor(:,2) & Dcor(:,2) < level(ind+1);
       Dcor_shelf=Dcor(SameShelf,:);
       v_min = roundoff([0;Dcor_shelf(:,3)]-Dcor(n,1));
       v_min = max(v_min(v_min<=0));
       v_max = roundoff([Dcor_shelf(:,1);W]-Dcor(n,3));
       v_max = min(v_max(v_max>=0));
       
       Cmat = Cgraph(Dcor,Dvol,size(Dcor,1),Sline);
       C0mat=Cmat(:,n) & ~SameShelf & Dvol~=0;
       C0 = sum(C0mat.*Dvol);
       
       v = roundoff([Dcor(~SameShelf,1)-Dcor(n,1);...
                     Dcor(~SameShelf,3)-Dcor(n,1);...
                     Dcor(~SameShelf,1)-Dcor(n,3);...
                     Dcor(~SameShelf,3)-Dcor(n,3);...
                     Dcor(~SameShelf,3)-Dcor(n,1)+300-Sline(1);...
                     Dcor(~SameShelf,1)-Dcor(n,3)-300+Sline(1)]);

       v = v(v_min<=v & v<=v_max);
       
       Dcor0=Dcor;
       for i=1:length(v)
           Dcor_tmp = Dcor;
           Dcor_tmp(n,[1,3])=Dcor_tmp(n,[1,3])+v(i);
           Cmat = Cgraph(Dcor_tmp,Dvol,size(Dcor,1),Sline);
           Cmat = Cmat(:,n) & Dvol~=0;
           C = sum(Cmat.*Dvol);
           if ~any(~C0mat & Cmat) && C0>C %New Conflict Are Not Allowed
                Dcor0=Dcor_tmp;
                C0=C;
                C0mat=Cmat;
           end
       end
       Dcor=Dcor0;
    end
end




function Dcor = ChipAround(Dcor,shelf,Caround)

for j=1:length(shelf(1).Dinfo(:,1)) 
    if Caround(shelf(1).Dinfo(j,1))  
        [shelf,Dcor]=SwapShelf(shelf,Dcor,1,2);    
        break;
    end;
end;

for j=1:length(shelf(end).Dinfo(:,1)) 
    if Caround(shelf(end).Dinfo(j,1))  
        [shelf,Dcor]=SwapShelf(shelf,Dcor,length(shelf)-1,length(shelf));    
        break;
    end;
end;


for i = 1:length(shelf)
    for j = 1:length(shelf(i).Dinfo(:,1)) 
        
        if Caround(shelf(i).Dinfo(j,1))            
            Rid = shelf(i).Rid(j);
            if j==1
               
                Dgap = Dcor(shelf(i).Dinfo(j+1,1),1)-Dcor(shelf(i).Dinfo(j,1),3);
                Dcor(shelf(i).Dinfo(j+1,1),1) = Dcor(shelf(i).Dinfo(j+1,1),1)-Dgap;
                Dcor(shelf(i).Dinfo(j+1,1),3) = Dcor(shelf(i).Dinfo(j+1,1),3)-Dgap;                               
            else
                inds = find(shelf(i).Rid==Rid);
                Linds = inds(inds<j); if ~isempty(Linds); Linds = Linds(end); end;
                Rinds = inds(inds>j); if ~isempty(Rinds); Rinds = Rinds(1); end;
                
                if ~isempty(Linds)
                    Dgap = Dcor(shelf(i).Dinfo(j,1),1)-Dcor(shelf(i).Dinfo(Linds,1),3);                    
                    Dcor(shelf(i).Dinfo(Linds,1),1) = Dcor(shelf(i).Dinfo(Linds,1),1)+Dgap;
                    Dcor(shelf(i).Dinfo(Linds,1),3) = Dcor(shelf(i).Dinfo(Linds,1),3)+Dgap;                                              
                end;
                
                if ~isempty(Rinds)
                    Dgap = Dcor(shelf(i).Dinfo(Rinds,1),1)-Dcor(shelf(i).Dinfo(j,1),3);                    
                    Dcor(shelf(i).Dinfo(Rinds,1),1) = Dcor(shelf(i).Dinfo(Rinds,1),1)-Dgap;
                    Dcor(shelf(i).Dinfo(Rinds,1),3) = Dcor(shelf(i).Dinfo(Rinds,1),3)-Dgap;                                              
                end;            
                
                
            end
            
            Dcor(shelf(i).Dinfo(j,1),2) = Dcor(shelf(i).Dinfo(j,1),2) + shelf(i).Dinfo(j,4);
            Dcor(shelf(i).Dinfo(j,1),4) = Dcor(shelf(i).Dinfo(j,1),4) + shelf(i).Dinfo(j,4);
        end        
        
    end;
end;
end

function [shelf,Dcor] = SwapShelf(shelf,Dcor,i, j)

if i>j
    temp=i;
    i=j; j=temp;
end;

DiffLevel = shelf(j).level(1)-shelf(i).level(1);
DiffHeight = shelf(j).height(1)-shelf(i).height(1);
s=DiffLevel>0;

for k=1:length(shelf(i).Dinfo(:,1))
    Dcor(shelf(i).Dinfo(k,1),2) = Dcor(shelf(i).Dinfo(k,1),2) +DiffLevel + s*DiffHeight;
    Dcor(shelf(i).Dinfo(k,1),4) = Dcor(shelf(i).Dinfo(k,1),4) +DiffLevel + s*DiffHeight;
end;

for k=1:length(shelf(j).Dinfo(:,1))
    Dcor(shelf(j).Dinfo(k,1),2) = Dcor(shelf(j).Dinfo(k,1),2) -DiffLevel -(~s)*DiffHeight;
    Dcor(shelf(j).Dinfo(k,1),4) = Dcor(shelf(j).Dinfo(k,1),4) -DiffLevel -(~s)*DiffHeight;
end;

for k=1: length(shelf(i).level)
    shelf(i).level(k)=shelf(i).level(k)+DiffLevel+s*DiffHeight;
end;

for k=1: length(shelf(j).level)
    shelf(j).level(k)=shelf(j).level(k)-DiffLevel-(~s)*DiffHeight;
end;

temp = shelf(i);
shelf(i) = shelf(j);
shelf(j)=temp;

end

function value = roundoff(value)
    value=round(value/0.001)*0.001;
end
