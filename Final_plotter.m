function Final_plotter(darg)
close all;
%clear;
clc;
pwd

pathh=strcat('SVRG_BB/Results_2022/',darg,'/');
find_par=1; %0 means accuracy and 1 means cost

reg = 1e-3;
fnt=20; %Font
lgft=23; %
msize=10; %Marker size
lsize=2; % Legend size 
isbar=1;
cost=[1e-2,0.4];
%cost = ['auto'];
%time=[0,2000];
time = 'auto';
%epoch=[0,20];
epoch='auto';
fname1={'SVRG','SVRG-BB','SVRG-2BB','SVRG-2D','SVRG-2BBS-ed','SVRG-2BBS-ec','SVRG-2BBS-e1'};%
method = {'svrg','svrg_bb','svrg_2bb','svrg_2d','svrg_2bbs_eta_decay','svrg_2bbs_eta_constant','svrg_2bbs_eta_one'};
Act_cost=[]; % Actual cost
Act_dv=[]; %Actual daviation 
best = zeros(7,14);
d = [];
for i = 1:7
    best(i,:) = other_best(strcat(pathh,char(method(i))), find_par, reg);
    Act_cost = [Act_cost best(i,6)];
    Act_dv = [Act_dv best(i,7)];
    d{i} = load(strcat(pathh,char(method(i)),sprintf('_%.1e_R_%.1e.mat',best(i,11),best(i,12))));
    ep = length(d{1}.S1.epoch);
end

f_opt = min(Act_cost);
f_opt1 = f_opt;
epsilon = 1e-5;
ind = find(Act_cost==f_opt);
fprintf('Method = %s',char(method(ind)));
g = [];
marker = ['o','<','*','s','p','v','^'];
% G1 = figure;%nexttile;
% G2 = figure;%nexttile;
% G3 = figure;%nexttile;

figure;
for j = 1:7
    
    %variance and CPU time
    g{j}.variance = d{j}.S1.variance(1:end);
    g{j}.time = d{j}.S1.otime(1:end);
    g{j}.epoch = d{j}.S1.epoch(1:end);
    g1 = plot(g{j}.time,g{j}.variance,'Marker',marker(j),'MarkerSize',msize,'LineWidth',lsize); hold on;
    if j==7
    %ylim(cost);
    title('Variance/Epoch');
    Gr1 = gca;
    set(gca,'Fontsize',fnt);
    Gr1.YScale = 'log';
    legend(fname1,'Location','NorthOutside','Orientation','horizontal','Box','on');
    hold off
    end
end    
disp('\n First');      

figure;
for j = 1:7
    %mean and standart daviation of cost
    g{j}.cost_mean = mean(d{j}.S1.ocost(1:end,:),2);
    g{j}.cost_std = isbar*std(d{j}.S1.ocost(1:end,:),[],2);
    %g2 = errorbar(G2,g{j}.epoch,g{j}.cost_mean,g{j}.cost_std,'Marker',marker(j),'MarkerSize',msize,'LineWidth',lsize); hold(G2,'on');
    g2 = plot(g{j}.time,(g{j}.cost_mean - f_opt + epsilon)/(1+f_opt),'Marker',marker(j),'MarkerSize',msize,'LineWidth',lsize); hold('on');
    
    if j==7
    %ylim(cost);
    title('Cost/Time');
    Gr1 = gca;
    set(gca,'Fontsize',fnt);
    Gr1.YScale = 'log';
    legend(fname1,'Location','NorthOutside','Orientation','horizontal','Box','on');
    hold off
    end
end  
disp('\n second');      
figure;
for j = 1:7
    %mean and standart daviation of optimal cost
    g{j}.opt_mean = abs(mean(d{j}.S1.opt_gap(1:end,:),2));
    g{j}.opt_std = abs(isbar*std(d{j}.S1.opt_gap(1:end,:),[],2));
    g3 = errorbar(g{j}.time,g{j}.opt_mean,g{j}.opt_std,'Marker',marker(j),'MarkerSize',msize,'LineWidth',lsize); hold('on');
    
    if j==7
    %ylim(cost);
    title('Optcost/Epoch');
    Gr1 = gca;
    set(gca,'Fontsize',fnt);
    Gr1.YScale = 'log';
    legend(fname1,'Location','NorthOutside','Orientation','horizontal','Box','on');
    hold off
    end
    
    
    %mean and standart daviation of validation cost
%     g{j}.vcost_mean = mean(d{j}.S1.vcost(1:end,:),2);
%     g{j}.vcost_std = std(d{j}.S1.vcost(1:end,:),[],2);
        
    
    %mean and standart daviation of Train ACC
%     g{j}.train_ac_mean = mean(d{j}.S1.train_ac(1:end,:),2);
%     g{j}.train_ac_std = std(d{j}.S1.train_ac(1:end,:),[],2);
    
    %mean and standart daviation of Train ACC
%     g{j}.val_ac_mean = mean(d{j}.S1.val_ac(1:end,:),2);
%     g{j}.val_ac_std = std(d{j}.S1.val_ac(1:end,:),[],2);
    
    %Lacotte style
    %g{j}.Lacotte_style = (g{j}.cost_mean - f_opt + epsilon)/(1+f_opt);
    

end

% for
% 
% e=1;
%  c = 1e-0;
% c1 = (mean(d.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c2 = (mean(d1.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c3 = (mean(d2.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c4 = (mean(d4.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c5 = (mean(d5.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c6 = (mean(d6.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% 
% 
% c1 = (mean(d.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c2 = (mean(d1.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c3 = (mean(d2.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c4 = (mean(d4.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c5 = (mean(d5.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% c6 = (mean(d6.S1.ocost(1:end,:),2) - f_opt + epsilon)/(1+f_opt1)*c;
% 
% % s1 = isbar*std(d.S1.ocost(1:end,:),[],2)/(1+f_opt1);
% % s2 = isbar*std(d1.S1.ocost(1:end,:),[],2)/(1+f_opt1);
% % s3 = isbar*std(d2.S1.ocost(1:end,:),[],2)/(1+f_opt1);
% % s4 = isbar*std(d4.S1.ocost(1:end,:),[],2)/(1+f_opt1);
% % % s5 = isbar*std(d5.S1.ocost(1:end,:),[],2)/(1+f_opt1);
% % disp('svrg');
% % [d.S1.ocost c1]
% % disp('svrg-bb');
% % [d1.S1.ocost c2]
% % disp('svrg-2bb');
% % [d2.S1.ocost c3]
% % disp('svrg-dh');
% % [d4.S1.ocost c4]
% % disp('svrg-new step');
% %[d5.S1.ocost c5]
%  figure; hold on;
% plot(d.S1.epoch(1:endd,:),(c1),'r--*','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d1.S1.epoch(1:endd,:),(c2),'b--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d2.S1.epoch(1:endd,:),(c3),'g-s','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d3.S1.epoch(2:endd,:),(d3.S1.variance(2:endd,:)),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d4.S1.epoch(1:endd,:),(c4),'m--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d5.S1.epoch(1:endd,:),(c5(1:endd)),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d6.S1.epoch(1:endd,:),(c6(1:endd)),'k--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% xlabel('Epoch','Fontsize',lgft)
% ylabel('Original cost','Fontsize',lgft)
% ylim(cost);
% xlim(epoch);
% G = gca;
% set(gca,'Fontsize',fnt);
% G.YScale = 'log';
% lgh=legend(fname1,'location','northeast','NumColumns',5);
% % %saveas (gcf, 'Variance_Epoch' , 'epsc' )
% 
% 
% e=1;
% 
% figure;hold on;
% %title("Optimal gap Vs. Epoch");
% errorbar(d.S1.epoch(e:endd,:),((mean(abs(d.S1.ocost(e:endd,:)),2))-f_opt+epsilon)/(1+f_opt),isbar*(std(d.S1.ocost(e:endd,:),[],2))/(1+f_opt1),'r--*','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(d1.S1.epoch(e:endd,:),((mean(abs(d1.S1.ocost(e:endd,:)),2))-f_opt+epsilon)/(1+f_opt),isbar*(std(d1.S1.ocost(e:endd,:),[],2))/(1+f_opt1),'b--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(d2.S1.epoch(e:endd,:),((mean(abs(d2.S1.ocost(e:endd,:)),2))-f_opt+epsilon)/(1+f_opt),isbar*(std(d2.S1.ocost(e:endd,:),[],2))/(1+f_opt1),'g-s','MarkerSize',msize,'LineWidth',lsize);hold on;
% %errorbar(d3.S1.epoch(e:endd,:),abs((mean(abs(d3.S1.ocost(e:endd,:)),2))-f_opt+epsilon)/(1+f_opt),isbar*(std(d3.S1.ocost(e:endd,:),[],2)),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(d4.S1.epoch(e:endd,:),((mean(abs(d4.S1.ocost(e:endd,:)),2))-f_opt+epsilon)/(1+f_opt),isbar*(std(d4.S1.ocost(e:endd,:),[],2))/(1+f_opt1),'m--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(d5.S1.epoch(e:endd,:),((mean(abs(d5.S1.ocost(e:endd,:)),2))-f_opt+epsilon)/(1+f_opt),isbar*(std(d5.S1.ocost(e:endd,:),[],2))/(1+f_opt1),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% xlabel('Epoch','Fontsize',lgft)
% ylabel('relative error','Fontsize',lgft)
% ylim(cost);
% xlim(epoch);
% legend(fname1,'location','northeast','NumColumns',1);
% G = gca;
% G.YScale = 'log';
% set(gca,'Fontsize',fnt);
% %legendmarkeradjust(8,12)
% saveas (gcf, 'opt_gap_Epoch' , 'epsc' )
% 
% 
% figure;hold on;
% %title("Optimal gap Vs. time");
% e=1;
% errorbar(d.S1.otime(e:endd,:),((mean(abs(d.S1.ocost(e:endd,:)),2))-f_opt),isbar*(std(d.S1.ocost(e:endd,:),[],2))/(1+f_opt),'r--*','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(d1.S1.otime(e:endd,:),((mean(abs(d1.S1.ocost(e:endd,:)),2))-f_opt),isbar*(std(d1.S1.ocost(e:endd,:),[],2))/(1+f_opt),'b--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(d2.S1.otime(e:endd,:),((mean(abs(d2.S1.ocost(e:endd,:)),2))-f_opt),isbar*(std(d2.S1.ocost(e:endd,:),[],2))/(1+f_opt),'g-s','MarkerSize',msize,'LineWidth',lsize);hold on;
% %errorbar(d3.S1.otime(e:endd,:),abs((mean(abs(d3.S1.ocost(e:endd,:)),2))-f_opt),isbar*(std(d3.S1.ocost(e:endd,:),[],2))/(1+f_opt),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(d4.S1.otime(e:endd,:),((mean(abs(d4.S1.ocost(e:endd,:)),2))-f_opt),isbar*(std(d4.S1.ocost(e:endd,:),[],2))/(1+f_opt),'m--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% errorbar(mean(d5.S1.otime(e:endd,:),2),((mean(abs(d5.S1.ocost(e:endd,:)),2))-f_opt),(isbar*std(d5.S1.ocost(e:endd,:),[],2))/(1+f_opt),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% xlabel('time','Fontsize',lgft)
% ylabel('relative error','Fontsize',lgft)
% ylim(cost);
% xlim(time);
% % xticks(20:40:340);
% G = gca;
% G.YScale = 'log';
% set(gca,'Fontsize',fnt);
% legend(fname1,'location','northeast');
% %legendmarkeradjust(8,12)
% saveas (gcf, 'opt_gap_time' , 'epsc' )
% 
% %endd=20;
% fig=figure;hold on;
% %title("Variance Vs. Epoch");
% plot(d.S1.epoch(2:endd,:),(d.S1.variance(2:endd,:)),'r--*','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d1.S1.epoch(2:endd,:),(d1.S1.variance(2:endd,:)),'b--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d2.S1.epoch(2:endd,:),(d2.S1.variance(2:endd,:)),'g-s','MarkerSize',msize,'LineWidth',lsize);hold on;
% %plot(d3.S1.epoch(2:endd,:),(d3.S1.variance(2:endd,:)),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d4.S1.epoch(2:endd,:),(d4.S1.variance(2:endd,:)),'m--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d5.S1.epoch(2:endd,:),mean(d5.S1.variance(2:endd,:),2),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% xlabel('Epoch','Fontsize',lgft)
% ylabel('Variance','Fontsize',lgft)
% ylim(cost);
% xlim(epoch);
% G = gca;
% set(gca,'Fontsize',fnt);
% G.YScale = 'log';
% lgh=legend(fname1,'location','northeast','NumColumns',5);
% saveas (gcf, 'Variance_Epoch' , 'epsc' )
% %saveLegendToImage(fig, lgh, 'legend3', 'epsc');
% 
% 
% figure;hold on;
% %title("Variance Vs. Time");
% plot(d.S1.otime(2:endd,:),(d.S1.variance(2:endd,:)),'r--*','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d1.S1.otime(2:endd,:),(d1.S1.variance(2:endd,:)),'b--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d2.S1.otime(2:endd,:),(d2.S1.variance(2:endd,:)),'g-s','MarkerSize',msize,'LineWidth',lsize);hold on;
% %plot(d3.S1.otime(2:endd,:),(d3.S1.variance(2:endd,:)),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(d4.S1.otime(2:endd,:),(d4.S1.variance(2:endd,:)),'m--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% plot(mean(d5.S1.otime(2:endd,:),2),mean(d5.S1.variance(2:endd,:),2),'k:>','MarkerSize',msize,'LineWidth',lsize);hold on;
% xlabel('Time','Fontsize',lgft)
% ylabel('Variance','Fontsize',lgft)
% ylim(cost);
% xlim(time);
% % xticks(20:40:340);
% legend(fname1,'location','northeast');
% G = gca;
% G.YScale = 'log';
% set(gca,'Fontsize',fnt);
% saveas (gcf, 'Variance_Time' , 'epsc' )


%  figure;hold on;
% % title("Train Acc Vs. Epoch");
%  errorbar(d.S1.epoch(e:endd,:),mean(d.S1.train_ac(e:endd,:),2),isbar*std(d.S1.train_ac(e:endd,:),[],2),'r-*','MarkerSize',msize,'LineWidth',lsize);hold on;
%  errorbar(d.S1.epoch(2:endd,:),mean(d1.S1.train_ac(2:endd,:),2),isbar*std(d1.S1.train_ac(2:endd,:),[],2),'b--*','MarkerSize',msize,'LineWidth',lsize);hold on;
%  errorbar(d.S1.epoch(e:endd,:),mean(d2.S1.train_ac(e:endd,:),2),isbar*std(d2.S1.train_ac(e:endd,:),[],2),'g-s','MarkerSize',msize,'LineWidth',lsize);hold on;
%  errorbar(d.S1.epoch(2:endd,:),mean(d3.S1.train_ac(2:endd,:),2),isbar*std(d3.S1.train_ac(2:endd,:),[],2),'k--*','MarkerSize',msize,'LineWidth',lsize);hold on;
%  errorbar(d.S1.epoch(e:endd,:),mean(d4.S1.train_ac(e:endd,:),2),isbar*std(d4.S1.train_ac(e:endd,:),[],2),'m--o','MarkerSize',msize,'LineWidth',lsize);hold on;
%  errorbar(d.S1.epoch(e:endd,:),mean(d5.S1.train_ac(e:endd,:),2),isbar*std(d5.S1.train_ac(e:endd,:),[],2),'y--o','MarkerSize',msize,'LineWidth',lsize);hold on;
% 
%  xlabel('Epoch','FontSize',lgft)
%  ylabel('Train acc','FontSize',lgft)
%  ylim(acc);
%  xlim(epoch);
%  %legend(fname,'location','southeast');
% % G = gca;
% % G.YScale = 'log';
%  set(gca,'Fontsize',fnt);
% % saveas (gcf, 'Train_Acc_Epoch' , 'epsc' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure;hold on;
% % title("Val Acc Vs. Epoch");
%  errorbar(d.S1.epoch(e:endd,:),mean(d.S1.val_ac(e:endd,:),2),isbar*std(d.S1.val_ac(e:endd,:),[],2),'r-*','MarkerSize',msize,'LineWidth',lsize); hold on;
%  errorbar(d.S1.epoch(2:endd,:),mean(d1.S1.val_ac(2:endd,:),2),isbar*std(d1.S1.val_ac(2:endd,:),[],2),'b-*','MarkerSize',msize,'LineWidth',lsize);
%  errorbar(d.S1.epoch(e:endd,:),mean(d2.S1.val_ac(e:endd,:),2),isbar*std(d2.S1.val_ac(e:endd,:),[],2),'g-s','MarkerSize',msize,'LineWidth',lsize); hold on;
%  errorbar(d.S1.epoch(2:endd,:),mean(d3.S1.val_ac(2:endd,:),2),isbar*std(d3.S1.val_ac(2:endd,:),[],2),'k-*','MarkerSize',msize,'LineWidth',lsize);
%  errorbar(d.S1.epoch(e:endd,:),mean(d4.S1.val_ac(e:endd,:),2),isbar*std(d4.S1.val_ac(e:endd,:),[],2),'m--o','MarkerSize',msize,'LineWidth',lsize); hold on;
%  errorbar(d.S1.epoch(e:endd,:),mean(d5.S1.val_ac(e:endd,:),2),isbar*std(d5.S1.val_ac(e:endd,:),[],2),'k:>','MarkerSize',msize,'LineWidth',lsize); hold on;
%  xlabel('Epoch','FontSize',lgft)
%  ylabel('Test acc','FontSize',lgft)
%   ylim(accc);
%  xlim(epoch);
%  %legend(fname,'location','southeast');
% % %G = gca;
% % %G.YScale = 'log';
%  set(gca,'Fontsize',fnt);
% % saveas (gcf, 'Val_Acc_Epoch' , 'epsc' )

end