xdata = [0:1:359];
x0 = [0 0 1 20 1];
fun = @(x)get_vonmises(x,xdata);

N = 10000;
R=zeros(length(xdata),N);
C0=linspace(0,180,N);
for i=1:N
    x =x0;
    x(1) = C0(i);
    R(:,i) = fun(x);    
end
 

figure; imagesc(R)
%--------- I have to load a from plot_orisc_crses.m
[ahat,ci] = gamfit(a);
x = [0:0.1:20];
Y = gampdf(x,ahat(1),ahat(2));
P = Y/sum(Y);
figure; plot(x,P);
hold on; 
m=hist(a,x);
plot(x,m/sum(m))

n = round(P*N);

n20=n;
n40=n;

R40 =R;
R20=R;
rem40 = randperm(N);
rem20 = randperm(N);
for i = 1 : length(n)
    s = x(i);
    n1 = n40(i);    
    if n1>0
        R40(:,rem40(1:n1))= s*R40(:,rem40(1:n1));
        rem40 = rem40(n1+1:end);
    end
    n2 = n20(i);    
    if n2>0
        R20(:,rem20(1:n2))= s*R20(:,rem20(1:n2));
        rem20 = rem20(n2+1:end);
    end
end
    

 figure; plot(R40(90,:));
 bin = linspace(0,360,31);
 R40x = zeros(size(R,1),30);
 R100x = zeros(size(R,1),30);
 for i = 1:30
     inx = C0>=bin(i) & C0<bin(i+1);
     R40x(:,i)=mean(R40(:,inx),2);     
     R100x(:,i)=mean(R(:,inx),2);     
 end
%  
figure; hold on;
plot(bin(1:end-1)+0.5*(bin(2)-bin(1)),R100x(90,:)/max(R100x(90,:)),'r.-','linewidth',2);
plot(bin(1:end-1)+0.5*(bin(2)-bin(1)),R40x(90,:)/max(R100x(90,:)),'b.-','linewidth',2)
set(gca,'fontsize',20)
xlim([0 180])

figure; hold on;
  figure; plot(C0,R40(90 ,:)) 
set(gca,'fontsize',20)
xlim([0 180])
    
  figure; plot(R40([90 120],:)')   
 figure; plot(R20([90 120],:)')   
  figure; plot([R40(90 ,:)' R20(90 ,:)'])   
 


