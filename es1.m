clear all 
close all 
clc

load data1a.mat

N = length(u);

Ne = 5000; 
Nv = N-Ne;

u = u - mean(u);
y = y - mean(y);

ue = u(1:Ne);
ye = y(1:Ne);
uv = u(Ne+1:end);
yv = y(Ne+1:end);

DATAe = [ye ue];
DATAv = [yv uv];

N0 = 15;

%% ARX
for nk = 1:3,
    for k = 1:10,
        na = k; nb = k;
        n = na;
        n_arx(k) = n;
        model_arx = arx(DATAe,[na,nb,nk]);
        %figure,resid(DATAe,model_arx,'CORR',30)
        yh_arx = compare(DATAv, model_arx, inf);
        MSE = norm(yv(N0+1:end)-yh_arx(N0+1:end))^2/(Nv-N0);
        RMSE_arx(k,nk)= sqrt(MSE);
        AIC_arx(k,nk) = n*2/(Nv-N0)+log(MSE);
    end
    close all 
end

RMSE_arx,AIC_arx

figure,plot(n_arx,RMSE_arx,'o-'),grid on, title('trade off RMSE-ARX')
legend('nk = 1','nk = 2','nk = 3')
figure,plot(n_arx,AIC_arx,'o-'),grid on, title('trade off AIC-ARX')
legend('nk = 1','nk = 2','nk = 3')

%% ARMAX
for nk = 1:3,
    for k = 1:10,
        na = k; nb = k; nc = k;
        n = na;
        n_armax(k) = n;
        model_armax = armax(DATAe,[na,nb,nc,nk]);
        %figure,resid(DATAe,model_armax,'CORR',30)
        yh_armax = compare(DATAv, model_armax, inf);
        MSE = norm(yv(N0+1:end)-yh_armax(N0+1:end))^2/(Nv-N0);
        RMSE_armax(k,nk)= sqrt(MSE);
        AIC_armax(k,nk) = n*2/(Nv-N0)+log(MSE);
    end
    close all 
end

RMSE_armax,AIC_armax

figure,plot(n_armax,RMSE_armax,'o-'),grid on, title('trade off RMSE-ARMAX')
legend('nk = 1','nk = 2','nk = 3')
figure,plot(n_armax,AIC_armax,'o-'),grid on, title('trade off AIC-ARMAX')
legend('nk = 1','nk = 2','nk = 3')

%% OE
for nk = 1:3,
    for k = 1:10,
        nb = k; nf = k;
        n = nf;
        n_oe(k) = n;
        model_oe = oe(DATAe,[nb,nf,nk]);
        %figure,resid(DATAe,model_oe,'CORR',30)
        yh_oe = compare(DATAv, model_oe, inf);
        MSE = norm(yv(N0+1:end)-yh_oe(N0+1:end))^2/(Nv-N0);
        RMSE_oe(k,nk)= sqrt(MSE);
        AIC_oe(k,nk) = n*2/(Nv-N0)+log(MSE);
    end
    close all 
end

RMSE_oe,AIC_oe

figure,plot(n_oe,RMSE_oe,'o-'),grid on, title('trade off RMSE-OE')
legend('nk = 1','nk = 2','nk = 3')
figure,plot(n_oe,AIC_oe,'o-'),grid on, title('trade off AIC-OE')
legend('nk = 1','nk = 2','nk = 3')

%% best model

best_arx = arx(DATAe,[9,9,1]);
best_armax = armax(DATAe,[3,3,3,1]);
best_oe = oe(DATAe,[4,4,1]);
best_model = best_armax;
present(best_model)
[A,B,C,D] = polydata(best_model)
F2 = tf(B,A,1,'variable','z^-1')
F1 = tf(C,A,1,'variable','z^-1')

%% ARX(2,3,2)

clear all 
close all 
clc

load data1a.mat

N = length(u);

u = u - mean(u);
y = y - mean(y);

na = 2;
nb = 3;
nk = 2;

n = na+nb; 

t_min = max([na+1,nb+nk]);
Y = y(t_min:N);
PHI = [];

for t = t_min:N,
    PHI = [PHI;-y(t-1) -y(t-2) u(t-2) u(t-3) u(t-4)];
end
theta_LS = PHI\Y

alpha = 100;
theta_RLS3 = zeros(n,N);
V = alpha*eye(n);
phi = [];

for t = t_min:N,
    phi = [-y(t-1) -y(t-2) u(t-2) u(t-3) u(t-4)]';
    beta = 1+phi'*V*phi;
    V = V-inv(beta)*V*phi*phi'*V;
    K = V*phi;
    eps = y(t)-phi'*theta_RLS3(:,t-1);
    theta_RLS3(:,t)=theta_RLS3(:,t-1)+K*eps;
end
theta_RLS3_f = theta_RLS3(:,end)

figure,plot(1:N,theta_RLS3)
for i = 1:size(theta_LS),
    refline(0,theta_LS(i))
end

% EUI2
sigma = sqrt(diag(inv(PHI'*PHI))); 
eps = 1;
theta_min = theta_LS-eps*sigma;
theta_max = theta_LS+eps*sigma;
EUI2 = [theta_min theta_max]

%PUI 
alpha = norm(Y-PHI*theta_LS);
if alpha^2>eps^2,
    disp('PUI cannot be computed')
end

%EUIinf

A = inv(PHI'*PHI)*PHI';
eps = 0.05; 

for k = 1:length(theta_LS),
    theta_Min(k,1) = A(k,:)*(Y-eps*sign(A(k,:))');
    theta_Max(k,1) = A(k,:)*(Y+eps*sign(A(k,:))');
end
 EUIinf = [theta_Min theta_Max]