clear;
clc;
close all;

load('pon_data_20201210.mat','B');

[row ,column] = size(B);
B_cut = zeros(1,1);
p = zeros(1,1);
for i = 1: row/2
    for j = 1:column
        if B(i*2,j) == 0
            break
        end
        B_cut(1,j) = B(i*2-1,j);  
        B_cut(2 ,j) = B(i*2,j);
%         [p(i,1:5),S(i)] = polyfit(B_cut(2,:),B_cut(1,:),4); %change
%         RR(i) = 1 - (S(i).normr/norm(B_cut(1,:) - mean(B_cut(1,:))))^2
    end
    [p(i,1:5),S(i)] = polyfit(B_cut(2,:),B_cut(1,:),4); %change
    RR(i) = 1 - (S(i).normr/norm(B_cut(1,:) - mean(B_cut(1,:))))^2
    
    
end
% figure;
% plot(B_cut(2,:),B_cut(1,:),'.');
% figure;

for i = 46: 46
    for j = 1:column
        if B(i*2,j) == 0
            break
        end
        B_cut(1,j) = B(i*2-1,j);  
        B_cut(2 ,j) = B(i*2,j);
        
    end

    
    % change
    modelFun = @(b,x) b(1) + b(2)*x.^1 +b(3)*x.^2+b(4)*x.^3+b(5)*x.^4 ;
%     modelFun = @(b,x) b(1) + b(2)*x.^1 +b(3)*x.^2;
    start = [1 ;10 ;10 ;10 ;10 ];
    nlm = fitnlm(B_cut(2,:), B_cut(1,:),modelFun,start);
%     xx = linspace(-30,30)';
%     line(xx,predict(nlm,xx),'linestyle','--','color','k')
    R_squ_get(i) = nlm.Rsquared.Adjusted;
    R_b_get(i,1:5) = nlm.Coefficients.Estimate;

    
end
figure;
plot(B_cut(2,:),B_cut(1,:),'.')
figure;

hold on;
plot(RR);
title('normal R')
% hold on;
figure;
plot(R_squ_get);
title('ML R')

figure;
x2 = B(92,1:end);
y2 = B(91,1:end);  %change odd


for j = 1:column
        if B(i*2,j) == 0
            break
        end
        y(:,j) = y2(:,j);  
        x(:,j) = x2(:,j);
        
end

plot(x,y,'.');
x_lin = -30:0.1:30;
[y2,delta] = polyval(p(46,:),x_lin,S(46));  % change

hold on;
plot(x_lin,y2,'r');
plot(x_lin,y2+2*delta,'m--',x_lin,y2-2*delta,'m--');




modelFun = @(b,x) b(1) + b(2)*x.^1 +b(3)*x.^2+b(4)*x.^3+b(5)*x.^4 ;

start = [1 ;10 ;10 ;10 ;10 ];
% opts = statset('Display','iter')
% nlm = fitnlm(x,y,modelFun,start,'Options',opts);
 nlm = fitnlm(x,y,modelFun,start);
xx = linspace(-30,30)';
hold on;
line(xx,predict(nlm,xx),'linestyle','--','color','k')
R_squ_get = nlm.Rsquared.Ordinary;


[ypred,ypredci] = predict(nlm,xx,'Simultaneous',true);
plot(x,y,'k.',xx,ypred,'b-',xx,ypredci,'r:')
xlabel('x') 
ylabel('y')
% ylim([-150 350])
legend({'Data','Weighted fit','95% Confidence Limits'}, ...
    'location','SouthEast')
figure();
r = nlm.Residuals.Standardized;
plot(x,r,'b^')
xlabel('x')
ylabel('Standardized Residuals')



%%
% 
% [ypred,ypredci] = predict(nlm,xx,'Simultaneous',true);
% plot(x,y,'k.',xx,ypred,'b-',xx,ypredci,'r:')
% xlabel('x') 
% ylabel('y')
% % ylim([-150 350])
% legend({'Data','Weighted fit','95% Confidence Limits'}, ...
%     'location','SouthEast')
% figure();
% r = nlm.Residuals.Standardized;
% plot(x,r,'b^')
% xlabel('x')
% ylabel('Standardized Residuals')
% 
% figure;
% plot(RR);
% 
% R_squ_get = nlm.Rsquared.Ordinary;
% p_value_get = nlm.Coefficients.pValue;
% p_value_2 = coefTest(nlm);
%         
