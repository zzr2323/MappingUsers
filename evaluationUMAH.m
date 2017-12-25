function evaluationUMAH(dir,at)

Name1=[dir,'/V_30_3_5.mat'];
load(Name1);
V1=V;
name=[dir,'/Hx.txt'];
load(name); 
Hx=spconvert(Hx);
name=[dir,'/Hy.txt'];
load(name);
Hy=spconvert(Hy);
LNum=677;
name=[dir,'/xgt.txt'];
load(name);
xgt=spconvert(xgt);
evaluationNum=2033;

XNum = size(Hx,1); 
YNum = size(Hy,1);
XunLNum=XNum-LNum;
YunLNum=YNum-LNum;

% name1=[dir,'/d3_lm.txt'];
% fw1=fopen(name1,'w+');
% name2=[dir,'/precision3_lm.txt'];
% fw2=fopen(name2,'w+');
for d=100:100:3000
    S=zeros(XunLNum,YunLNum);
    XY=[V1(LNum+1:XNum,:);V1(XNum+LNum+1:XNum+YNum,:)];
    D=1-pdist(XY(:,1:d),'cosine');
    Z=squareform(D);
    for i=1:XunLNum
        for j=1:YunLNum
             S(i,j)=Z(i,XunLNum+j);
        end    
    end

    [A,Idx]=sort(S,2,'descend'); 

    hit=0;
    for i=1:evaluationNum
        if xgt(i,1)~=0
            if any(Idx(i,:)==xgt(i,1))
                hit=hit+1/find(Idx(i,:)==xgt(i,1));
%             if any(Idx(i,1:at)==xgt(i,1))
%                 hit=hit+(at+1-find(Idx(i,1:at)==xgt(i,1)))/at;
            end        
        end
    end

    precision=hit/evaluationNum;
    disp(d);
    disp(precision);
%     result=[num2str(d),'\n'];
%     fprintf(fw1,result);
%     result=[num2str(precision),'\n'];
%     fprintf(fw2,result);
end
end
















