clear all;clc;

%导入数据
disease_miRNA_interaction=xlsread('miRNAdisease.xlsx');
SSweighted=textread('疾病语义类似性加权矩阵.txt');
SS=textread('疾病语义类似性矩阵.txt');
FSweighted=textread('miRNA功能类似性加权矩阵.txt');
FS=textread('miRNA功能类似性矩阵.txt');

%计算邻接矩阵A （已知相关为1 否则为0）
for x=1:5430
    i=disease_miRNA_interaction(x,1);
    j=disease_miRNA_interaction(x,2);
    A(i,j)=1;
end
AA=A';

%KD (疾病之间的高斯核剖面相似性矩阵)
%计算gamma_disease
sumnorm2_d=0;
for i=1:383
    norm2_d=norm(AA(i,:),2);
    sumnorm2_d=sumnorm2_d+norm2_d;
end
gamma_disease=1/sumnorm2_d;

%计算KD
for i=1:383
    for j=1:383
        delta_d = norm((AA(i,:)-AA(j,:)),2);
        KD(i,j)=exp(-gamma_disease*delta_d);
    end
end
         
%KM（miRNA之间的高斯核剖面相似性矩阵）
%计算gamma_miRNA
sumnorm2_m=0;
for i=1:495
    norm2_m=norm(A(i,:),2);
    sumnorm2_m=sumnorm2_m+norm2_m;
end
gamma_miRNA=1/sumnorm2_m;

%计算KM
for i=1:495
    for j=1:495
        delta_m = norm((A(i,:)-A(j,:)),2);
        KM(i,j)=exp(-gamma_miRNA*delta_m);
    end
end

%计算疾病的综合相似性矩阵 inte_sim_disease
for i=1:383
    for j=1:383
        if SSweighted(i,j)==1
            inte_sim_disease(i,j)=SS(i,j);
        else
            inte_sim_disease(i,j)=KD(i,j);
        end
    end
end

%计算miRNA的综合相似性矩阵 inte_sim_miRNA
for i=1:495
    for j=1:495
        if FSweighted(i,j)==1
            inte_sim_miRNA(i,j)=FS(i,j);
        else
            inte_sim_miRNA(i,j)=KM(i,j);
        end
    end
end

%求初始向量p0
for i=1:383
    t=0;
    xx=[];
    for j=1:495
        if AA(i,j)==1
            t=t+1;     %疾病di有多少个已知相关的miRNA 
            xx=[xx,j]; %与疾病di已知相关的miRNA的编号
            score1(i,j)=1;
        else
            score1(i,j)=0;
        end
    end
    for j1=1:495
        if score1(i,j1)==0
            for j2=xx
                score1(i,j1)=score1(i,j1)+inte_sim_miRNA(j1,j2);
            end
            score1(i,j1)=score1(i,j1)/t;
        end
    end
end

for i=1:495
    t1=0;
    xx1=[];
    for j=1:383
        if A(i,j)==1
            t1=t1+1;     %给定miRNA有多少种已知相关的疾病 
            xx1=[xx1,j]; %与给定的miRNA已知相关的疾病的编号
            score2(i,j)=1;
        else
            score2(i,j)=0;
        end
    end
    for j1=1:383
        if score2(i,j1)==0
            for j2=xx1
                score2(i,j1)=score2(i,j1)+inte_sim_disease(j1,j2);
            end
            score2(i,j1)=score2(i,j1)/t1;
        end
    end
end

%对两种评分进行综合，得到初始向量p0
score11=score1';
for i=1:495
    for j=1:383
        score(i,j)=(score11(i,j)+score2(i,j))/2;
    end
end

%计算W，对miRNA的综合相似性矩阵进行标准化
D=diag(sum(inte_sim_miRNA,2).^(-1/2)); %生成对角阵
W=D.*inte_sim_miRNA.*D;

%进行随机游走

d=input('please input a number between 1 to 383: '); %输入疾病的编号
p0=score(:,d);

%设置重启概率r
r=0.8;
p1=(1-r)*W*p0+r*p0;
p2=(1-r)*W*p1+r*p0;
while norm(p1-p2,2)>1e-6
    p1=p2;
    p2=(1-r)*W*p1+r*p0;
end

%去掉种子结点
for i=1:495
    if A(i,d)==1
        p2(i,1)=0;
    else
        p2(i,1)=p2(i,1);
    end
end

[grade,rank]=sort(p2,'descend');
result=[rank,grade];

