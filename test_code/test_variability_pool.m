clear all
close all



dtype='Xsel'
iexp_type=1;
exp_type='AN'
ises=15

[X,mX,stdX,oris,cons] = sortdata(dtype,ises,exp_type);

C = 100;
O = [0 30];

I1 = find(oris==O(1) & cons==C);
I2 = find(oris==O(2) & cons==C);

X1 = [X{1,I1}; X{2,I1}]; 
X2 = [X{1,I2}; X{2,I2}];
D = [X1; X2];
evts = [ones(size(X1,1),1); -1*ones(size(X1,1),1)];

% X = bsxfun(@rdivide, X, sqrt(sum(X.^2,1)));

Ncv=10;
NC = size(D,2);
Ws = zeros(NC+1,Ncv+1);
inxs_cv = get_Kfoldcvinxs(evts,Ncv,false);

lam1=0.1;

for icv = 1 : Ncv+1
    if icv == Ncv+1
        inxsm = 1:size(D,1);
    else
        inxsm = inxs_cv.inxs_test{icv};
    end
    trdat = D(inxsm,:);
    trlab = evts(inxsm);
    conf.W=[];conf.btrain = true;
    if icv==Ncv+1
        conf.lambda1=lam1*Ncv;
    else
        conf.lambda1=lam1;
    end
    conf.lambda2=0;
    [W, out_tr] = sl_smlr(trdat', trlab,conf);
    Ws(:,icv)=W(:,1);
end

figure; hold on;
plot(Ws(:,1:Ncv));
plot(Ws(:,end),'k', 'LineWidth',2); 
mW = mean(Ws,2);

W = Ws(:,end);
y= D*W(1:end-1)+W(end);

cl=[5 17 43    8 12  29 35 36  45]%1
y1 =D(:,cl)*W(cl);
figure; plot(y1)


D1 = D(:,cl);
yp0 = bsxfun(@plus,D1, [0:length(cl)-1]);
yp = bsxfun(@times,D1, W(cl)');

m = (mean(X2,1)-mean(X1,1));
s1 = sqrt(var(X1,0,1));
s2 = sqrt(var(X2,0,1));
dp = (mean(X2,1) - mean(X1,1))./sqrt(var(X1,0,1)+var(X2,0,1));
[I J]=sort(W(1:end-1),'ascend');
nm = (m-min(m))/(max(m)-min(m));
ns1 = (s1-min(s1))/(max(s1)-min(s1));
ns2 = (s2-min(s2))/(max(s2)-min(s2));
ndp = (dp-min(dp))/(max(dp)-min(dp));
nI  = (I-min(I))/(max(I)-min(I));

figure; plot([nI ndp(J)' nm(J)' ns1(J)' ns2(J)'])


% figure; plot(yp0)
figure; 
plot([sum(yp(:,:),2) yp(:,1:3)-3])





