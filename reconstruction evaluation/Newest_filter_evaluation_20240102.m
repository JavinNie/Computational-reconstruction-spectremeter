clear;
clc;
load('200_0.3.mat')
A=T';
% load('Filter.mat')
%%A滤波器矩阵，MM*NN，MM：波长采样数；NN：折射率采样数



[MM,NN]=size(A); %波长采样数MM,折射率采样数NN；
psi=A;%加载滤波器

b=linspace(1,MM,MM);%声明一个和波长采样数等长的序列，不知道FILTER怎么存的，所以没写size(filter,1)
T=0;%延拓,忽略
psi=[rand(T,size(psi,2))/10;psi;rand(T,size(psi,2))/10].';  %抗噪能力对PSI非常敏感 对PHI不敏感（PHI由算法确定） 条件数一次反比噪声容限
aaa=CM_Correlation(psi,NN);
A = [0:1/(MM+2*T-1):1];
A = exp(-(A-0.5).^2/0.0001);%随便给了一个模拟的光谱向量，不关键，主要是为了获取下一步的L
wavlv=ceil(log2(MM))-4;%小波分解级数,里面数字也是波长采样数
[wavA,L]=wavedec(A,wavlv,'db6');%选取小波基为多贝西六阶小波
[LoD,HiD,LoR,HiR] = wfilters('db6');%获取提升算法所需滤波矩阵
n=length(HiR);


%%%%%%%小波提升算法矩阵化%%%%%%
for k=1:wavlv
    i=k+1;
    vec1=ones(1,L(i));
    MAT1=diag(vec1);
    
    MAT2=zeros(L(i)-1,L(i));
    up_sam(k).MAT=MAT1(1,:);
    for j=1:L(i)-1
        up_sam(k).MAT=[up_sam(k).MAT;MAT2(j,:);MAT1(j+1,:)];     
    end
    idwt_l(k).MAT=[];
    idwt_h(k).MAT=[];
    vec_l=[fliplr(LoR) zeros(1,2*L(i)-1-n)];
    vec_h=[fliplr(HiR) zeros(1,2*L(i)-1-n)];
    for j=1:2*L(i)+10
        vec2_l=circshift(vec_l',j-n)';
        vec2_h=circshift(vec_h',j-n)';
        idwt_l(k).MAT=[idwt_l(k).MAT;vec2_l];
        idwt_h(k).MAT=[idwt_h(k).MAT;vec2_h];
    end
    m=2*L(i)+n-L(i+1)-2;
    P(k).MAT=[zeros(L(i+1),floor(m/2)) eye(L(i+1)) zeros(L(i+1),ceil(m/2))];
    WL(k).MAT=P(k).MAT*idwt_l(k).MAT*up_sam(k).MAT;
    WH(k).MAT=P(k).MAT*idwt_h(k).MAT*up_sam(k).MAT;
%     for j=1:L(i+1)
%         if(j<L(i+1))
%             vec2=[zeros(1,j-1) LoR zeros(1,2*L(i)-j-11)];
%  
%         else
%             vec2=[LoR(12) zeros(1,2*L(i)-12) LoR(1:11)];
%         end
%         idwt(k).MAT=[idwt(k).MAT;vec2];
%    end
    if(k==1)
        phi_l=WL(k).MAT;
        phi_h=WH(k).MAT;
    else
        phi_l=WL(k).MAT*phi_l;
        phi_h=WL(k).MAT*phi_h;
    end
end
%%
% %小波分解四级是这样，见前文的wavlv变量，要是五次的话就是phi_h2=WL(5).MAT*WL(4).MAT*WL(3).MAT*WH(2).MAT依此类推，主要是懒得改了
% phi_h2=WL(4).MAT*WL(3).MAT*WH(2).MAT;
% phi_h3=WL(4).MAT*WH(3).MAT;
% phi_h4=WH(4).MAT;
% %编码矩阵，觉得高频不重要或者没什么高频分量的话，可以像这样只把低频部分对应的编码矩阵放进去，这样平滑一些
% phi=[phi_l phi_h phi_h2]; % phi_h3 phi_h4  
%% 自动按阶数完成编码矩阵计算
N_order=3;%% 手动修改，根据需要选择保留前几阶小波分量；
if(N_order>wavlv)
    N_order=wavlv;
end
phi=[phi_l phi_h];
for ii=1:N_order-2
    phi_temp=1;
    for jj=wavlv:-1:2
        phi_temp=phi_temp*WL(jj).MAT;
    end
    phi=[phi,phi_temp];
end


dictionary=psi*phi;%字典 编码矩阵和测量矩阵的积

for i=1:length(phi(1,:))
    vec(i)=1;%max(abs(dictionary(:,i)));
    dictionary(:,i) =dictionary(:,i)/vec(i);  %字典做了一次归一化，这里归一化向量直接是全1可忽略,没什么用
end


%%%生成测试波形%%%

%波形原型%
N=1;  %数量
for i=1:N
    p=randi(8);%光谱峰数，可以给个定值
    a=zeros(1,MM); 
    rand('seed',16875); %想让每次产生波形都一样，那就像这样固定随机数种子
    for j=1:p
        alpha=rand()*50;%振幅
        lamda=rand()*MM;%波长
        sigma=rand()*10+5;%峰宽
        a=a+alpha.*exp(-(lamda-b).^2./sigma^2);
    end
    testset(i,:)=a;
    testset(i,:) = (testset(i,:)/max(testset(i,:)));%归一化
end
testset = [testset(:,1).*ones(N,T) testset testset(:,MM).*ones(N,T)]';

%%分辨率测试波形，峰宽1nm-41nm%%
N=20;
for i=1:N
    p=1;
    a=zeros(1,MM); 
    rand('seed',16875);
    for j=1:p
        alpha=1;
        lamda=125;
        sigma=2*i-1;%
        a=a+alpha.*exp(-(lamda-b).^2./sigma^2);
    end
    testset_1(i,:)=a;
    testset_1(i,:) = (testset_1(i,:)/max(testset_1(i,:)));
end
testset_1 = [testset_1(:,1).*ones(N,T) testset_1 testset_1(:,MM).*ones(N,T)]';

%%动态范围测试波形，在波形原型基础上放缩幅值%%
N=41;
for i=1:N
    dynamic=i-21;
    p=randi(8);%1
    a=zeros(1,MM); 
    rand('seed',16875);
    for j=1:p
        alpha=rand()*50;%1
        lamda=rand()*MM;%125
        sigma=rand()*10+5;%10
        a=a+alpha.*exp(-(lamda-b).^2./sigma^2);
    end
    testset_2(i,:)=a;
    testset_2(i,:) = 10^(dynamic/10)*(testset_2(i,:)/max(testset_2(i,:)));
end
testset_2 = [testset_2(:,1).*ones(N,T) testset_2 testset_2(:,MM).*ones(N,T)]';
L= 1/2/max(max(eig(dictionary'*dictionary)));

%%%抗噪性测试波形，原型基础上加不同噪声
for i=1:1:20
    SNR=10*(i-1);
    y(:,i)=awgn(psi*testset,SNR);
    [x0(:,i),~,~ ] = omp(y(:,i), dictionary,50,1e-4);
    %[x0(:,i),~] = FISTA(dictionary,y(:,i),1,L,200);
end
x0 = vec'.*x0;
test_r0=phi*x0;
test_r0(test_r0<0)=0;



%%%%%求解%%%%%

%这里懒得写通过参数一键选择求解方式的代码了，提供了OMP和FISTA两种，前者直接求解抗噪性差，后者迭代算法重建效果差.根据需要自行注释
%上面抗噪性部分把波形产生和求解放一起了，其实应该分开


for i=1:1:20
    y(:,i)=psi*testset_1(:,i);
    [x1(:,i),~,~ ] = omp(y(:,i), dictionary,50,1e-4);
    %[x1(:,i),~] = FISTA(dictionary,y(:,i),1e-5/L,L,200);
end
x1 = vec'.*x1;
test_r1=phi*x1;
test_r1(test_r1<0)=0;

for k=20:-1:-20
    i=k+21;
    y(:,i)=psi*testset_2(:,i);
    [x2(:,i),~,~ ] = omp(y(:,i), dictionary,50,1e-4);
    %[x2(:,i),~] = FISTA(dictionary,y(:,i),1e-5/L,L,200);
end
x2 = vec'.*x2;
test_r2=phi*x2;
test_r2(test_r2<0)=0;

for i=1:20
    error0(i)=sum((testset(T+1:T+MM,1)-test_r0(T+1:T+MM,i)).^2)/MM;
end
for i=1:20
    error1(i)=sum((testset_1(T+1:T+MM,i)-test_r1(T+1:T+MM,i)).^2)/MM;
end
for i=1:41
    error2(i)=10^(4.2-0.2*i)*sum((testset_2(T+1:T+MM,i)-test_r2(T+1:T+MM,i)).^2)/MM;
end
figure();
x=linspace(-20,20,41);
plot(x,log10(error2),'b','LineWidth',2);xlabel('动态范围dB'); ylabel('归一化logMSE');

figure();
x=linspace(0,190,20);
plot(x,log10(error0),'b','LineWidth',2);xlabel('SNR/dB'); ylabel('logMSE');

figure();
x=linspace(1,41,20);
plot(x,log10(error1),'b','LineWidth',2);xlabel('分辨率/nm'); ylabel('logMSE');

figure();
plot(testset(T+1:T+MM,end),'b','LineWidth',2);hold on;xlabel('波长/nm');
plot(test_r0(T+1:T+MM,end),'r','LineWidth',2);legend('GroundTruth','Rebuild');

function [x, r, indexes] = omp(y, A, k, tol)
% OMP算法求解Ax=y中的稀疏解x
% 输入参数：
% y：观测向量，大小为m x 1
% A：测量矩阵，大小为m x n
% k：稀疏度
% tol：收敛阈值
% 输出参数：
% x：稀疏解，大小为n x 1
% r：残差，大小为m x 1
% indexes：选中的系数下标，大小为k x 1
[m, n] = size(A);
r = y;
x = zeros(n, 1);
indexes = zeros(k, 1);
for i = 1:k
    % 计算投影系数
    proj = abs(A' * r);
    [~, index] = max(proj);
    indexes(i) = index;
    % 求解最小二乘问题
    x(indexes(1:i)) = A(:, indexes(1:i)) \ y;
    % 计算残差
    r = y - A(:, indexes(1:i)) * x(indexes(1:i));
    % 判断是否达到收敛条件
    if norm(r) < tol
        break;
    end
end
end

function [x, obj_val] = FISTA(A, b, lambda, L, max_iter)
% 基于FISTA算法的L1正则化问题求解
% 输入参数：
% A：线性变换矩阵
% b：观测向量
% lambda：正则化参数
% L：Lipschitz常数
% max_iter：最大迭代次数
% 输出参数：
% x：L1最小化问题的解
% obj_val：目标函数值
% 初始化
[m, n] = size(A);
x = zeros(n, 1);
y = x;
xp = x ;
t = 1;
obj_val = zeros(1, max_iter);
% 迭代
for iter = 1:max_iter
    
    % 计算梯度
    grad = A' * (A * y - b);
    % 更新d
    d = y - L * grad;
    % 软阈值处理
    x = soft_threshold(d, lambda * L);
    
    % 计算步长
    tp = (1 + sqrt(1 + 4 * t^2)) / 2;
    
    % 更新y
    y = x + (t - 1) / tp * (x - xp);
    xp = x;
    t = tp;
    % 计算目标函数值
    xn=[x(1:25);x(26:50);x(51:90)];
    obj_val(iter) = norm(A * x - b)^2 / 2 + lambda * norm(xn, 1);
end
end
function [y] = soft_threshold(x, lambda)
% 软阈值函数
y = sign(x) .* max(abs(x) - lambda, 0);
end


function ave_corr=CM_Correlation(CM,filter_num)
%CM : each Colum represents a filter
%filter_num: number of filters
[l,c]=size(CM);
if c==filter_num
%     pass
elseif l==filter_num
    CM=CM';
else
    error('[error]CM_Correlation: size of CM is wrong')
end

Correlation=corr(CM);
eye_ind=eye(filter_num);
Correlation(logical(eye_ind))=0;
max_corr=max(Correlation);
ave_corr=mean(max_corr);
% max_corr=max(max_corr)
end