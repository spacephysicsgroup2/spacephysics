clear;
clc;
close all;

event1 = cell2mat(struct2cell(load('pon_data_20140112.mat','B')));
event2 = cell2mat(struct2cell(load('pon_data_20140219.mat','B')));
event3 = cell2mat(struct2cell(load('pon_data_20140325.mat','B')));
event4 = cell2mat(struct2cell(load('pon_data_20201210.mat','B')));
event5 = cell2mat(struct2cell(load('pon_data_20211223.mat','B')));

event1_t = cell2mat(struct2cell(load('pon_data_20140112.mat','B_teacher')));
event2_t = cell2mat(struct2cell(load('pon_data_20140219.mat','B_teacher')));
event3_t = cell2mat(struct2cell(load('pon_data_20140325.mat','B_teacher')));
event4_t = cell2mat(struct2cell(load('pon_data_20201210.mat','B_teacher')));
event5_t = cell2mat(struct2cell(load('pon_data_20211223.mat','B_teacher')));


theory1 = cell2mat(struct2cell(load('pon_data_theory_20201210.mat','B')));
theory2 = cell2mat(struct2cell(load('pon_data_theory_20211223.mat','B')));

Bz1 = cell2mat(struct2cell(load('pon_data_20140112.mat','Bz_pres')));

[R_b_1,R_sq_1] = y(event1);
[R_b_2,R_sq_2] = y(event2);
[R_b_3,R_sq_3] = y(event3);
[R_b_4,R_sq_4] = y(event4);
[R_b_5,R_sq_5] = y(event5);

[R_b_4_t,R_sq_4_t] = y(event4_t);
[R_b_5_t,R_sq_5_t] = y(event5_t);


[R_b_theory1,R_sq_theory1] = y(theory1);
[R_b_theory2,R_sq_theory2] = y(theory2);

[row ,column] = size(R_b_1);
out1 = err_2(R_b_theory1,R_b_4_t);
out1_2 = err_2(R_b_4_t,R_b_4);
out2 = err_2(R_b_theory2,R_b_5_t);%teacher 經驗
out2_2 = err_2(R_b_5_t,R_b_5); %teacher 經驗

figure;
plot(out1);
hold on;
plot(out1_2);
title('20201210');
xlabel('time');
ylabel('sum residual sum of squares');
legend('Teacher-Theory','Teacher-Experience');

figure;
plot(out2);
hold on;
plot(out2_2);
title('20211223');
xlabel('time');
ylabel('sum residual sum of squares');
legend('Teacher-Theory','Teacher-Experience');

figure;

new = [Bz1 R_b_1];


x = 1:100;
plot(R_sq_5);
hold on;
plot(R_sq_4);
plot(R_sq_3);
plot(R_sq_2);
plot(R_sq_1);
legend;





function [R_b_get,R_squ_get] = y(B)

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

for i = 1: row/2
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
end


function out = err_2(R_b_real,R_b_th)

[row ,column] = size(R_b_real);
late = 1;
for i = 1:row
    for x = -20:1:20
        y = polyval(R_b_real(i,:),x);
        y2 = polyval(R_b_th(i,:),x);
        err(late) = (y-y2)^2;
        late = late + 1;
    end
    out(i) = sum(err(:));
    late = 1;
end
end




