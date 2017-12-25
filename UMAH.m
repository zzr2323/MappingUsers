function UMAH_new(dir,FNum,mu,gamma,beta)

% dir£ºpath of your data files
% FNum£ºthe number of the features£¬e.g 3
% mu£¬gamma£¬beta the parameters of the experients

name=[dir,'/Hx.txt'];
load(name); 
Hx=spconvert(Hx);
name=[dir,'/Hy.txt'];
load(name);
Hy=spconvert(Hy);
LNum=677; % the labeled number of users across two networks 
Hf=cell(1,FNum);
for i=1:FNum
    filename=['Hf',num2str(i)];
    name=[dir,'/',filename,'.txt'];
    load(name);
    Hf{i}=spconvert(eval(filename));
end
Dvx=diag(sum(Hx,2));
Dvy=diag(sum(Hy,2));
Dvf=cell(1,FNum);
for i=1:FNum
    Dvf{i}=diag(sum(Hf{i},2));
end
Dex=diag(sum(Hx,1));
Dey=diag(sum(Hy,1));
Def=cell(1,FNum);
for i=1:FNum
    Def{i}=diag(sum(Hf{i},1));
end
Wx=speye(size(Hx,2));
Wy=speye(size(Hy,2));
Wb=cell(1,FNum);
for i=1:FNum
    Wb{i}=speye(size(Hf{i},2));
end
XNum = size(Hx,1);
ExNum=nnz(Hx);
YNum = size(Hy,1);
EyNum=nnz(Hy);
ENum=ExNum+EyNum;
for i=1:FNum
    ENum=ENum+nnz(Hf{i});
end
Dex=diag(diag(Dex).^-1);
Dey=diag(diag(Dey).^-1);
for i=1:FNum
    Def{i}=diag(diag(Def{i}.^-1));
end
Lhx=Dvx-Hx*Wx*Dex*Hx';
Lhy=Dvy-Hy*Wy*Dey*Hy';
Zeros1=zeros(size(Lhx,1),size(Lhy,2));
Zeros2=zeros(size(Lhy,1),size(Lhx,2));
Lxy=[Lhx,Zeros1;Zeros2,Lhy];

Wc=sparse(XNum+YNum,XNum+YNum);
for i=1:LNum
    for j=XNum+1:XNum+LNum
        if i==j-XNum
            Wc(i,j)=1;
        end
    end
end
for i=XNum+1:XNum+LNum
    for j=1:LNum
        if i-XNum==j
            Wc(i,j)=1;
        end
    end
end
Dc=diag(sum(Wc,2));
Lc=Dc-Wc;

Wc2ll=ones(LNum)-eye(LNum);
Wc2xx=zeros(XNum-LNum,XNum-LNum);
Wc2xy=zeros(XNum-LNum,YNum-LNum);
Wc2yy=zeros(YNum-LNum,YNum-LNum);
Wc2=[zeros(LNum),zeros(LNum,XNum-LNum),Wc2ll,zeros(LNum,YNum-LNum);...
    zeros(XNum-LNum,LNum),Wc2xx,zeros(XNum-LNum,LNum),Wc2xy;...
    Wc2ll,zeros(LNum,XNum-LNum),zeros(LNum),zeros(LNum,YNum-LNum);...
    zeros(YNum-LNum,LNum),Wc2xy',zeros(YNum-LNum,LNum),Wc2yy];
Wc2=sparse(Wc2);
Dc2=diag(sum(Wc2,2));
Lc2=Dc2-Wc2+speye(XNum+YNum);

wb=zeros(1,FNum);
for i=1:FNum
    wb(i)=gamma/FNum;
end
Lf=cell(1,FNum);
for i=1:FNum
    Lf{i}=Dvf{i}-Hf{i}*Wb{i}*Def{i}*Hf{i}';
end
count=1;
while(1)
    LF=sparse(XNum+YNum,XNum+YNum);
    for i=1:FNum
        LF=LF+Lf{i}*wb(i)^beta;
    end
    A=(Lxy+LF)/ENum+(mu/LNum)*Lc;
    [V,D]=eig(full(A),full(Lc2));
    Mb=zeros(1,FNum);
    for i=1:FNum
        Mb(i)=trace(V'*Lf{i}*V);
    end
    norm=1/(beta-1);
    wb_update=zeros(1,FNum);
    result='';
    for i=1:FNum
        sum1=0;
        for j=1:FNum
            sum1=sum1+(Mb(i)/Mb(j))^norm;
        end
        wb_update(i)=gamma/sum1;
        result=[result,' ',num2str(wb_update(i))];
    end
    disp([count,result]);
    count=count+1;
    flag=0;
    for i=1:FNum
        if abs(wb_update(i)-wb(i))>0.001
            flag=1;
            break;
        end
    end
    if flag
    	for i=1:FNum
            wb(i)=wb_update(i);
        end
    else
        Name1=[dir,'/V_',num2str(mu),'_',num2str(gamma),'_',num2str(beta),'.mat'];
        Name2=[dir,'/D_',num2str(mu),'_',num2str(gamma),'_',num2str(beta),'.mat'];
        save(Name1,'V','-v7.3');
        save(Name2,'D','-v7.3');
        break;
    end
end
end


