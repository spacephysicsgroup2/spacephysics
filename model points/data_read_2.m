clear;
clc;
close all;

T = readtable('20201210_theory.csv');
% T2 = readtable('Shue_Model_20211223.csv');
A = table2array(T);
% A2 = table2array(T2);

[row column] = size(A);
j = 1;
B = zeros(1,1);
Bz_pres = zeros(1,1);

for i = 1 : row-1
    
    if A(i,1) == A(i+1,1)
        
        B(2*A(i,1) -1,j) =  A(i,2);
        B(2*A(i,1) ,j) =  A(i,3);
        j = j+1;
        fprintf('%f\n',j);
    else 
        B(2*A(i,1) -1,j) =  A(i,2);
        B(2*A(i,1) ,j) =  A(i,3);
%         Bz_pres(A(i,1),1) = A(i,4);
%         Bz_pres(A(i,1),2) = A(i,5);
        
        j = 1;
        
    end
end
% Bz_pres(A(end,1),1) = A(end,4);
% Bz_pres(A(end,1),2) = A(end,5);


%teacher
% [row2 column2] = size(A2);
% j = 1;
% B_teacher = zeros(1,1);


% for i = 1 : row2-1
%     
%     if A2(i,1) == A2(i+1,1)
%         
%         B_teacher(2*A2(i,1) -1,j) =  A2(i,2);
%         B_teacher(2*A2(i,1) ,j) =  A2(i,3);
%         j = j+1;
%         fprintf('%f\n',j);
%     else 
%         B_teacher(2*A2(i,1) -1,j) =  A2(i,2);
%         B_teacher(2*A2(i,1) ,j) =  A2(i,3);
%        
%         j = 1;
%         
%     end
% end
% 
% 



C  = B(2,:,:);
save('pon_data_theory_20201210.mat','B');



