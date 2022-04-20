close all; clear all;
load 'Divergence_matrix.mat';

fig1=figure('Position',[100 1500 800 500],'Color',[1 1 1]); %[left bottom width height];   
sub = subplot(1,1,1);
h_a = 0.001;
h_b = 0.0005;
rgb = imagesc(0.2:h_a:1,  0:h_b:0.14, abs(Masyvas_diverg'));   
set(gca,'YDir','normal');
colormap(parula(5));
colorbar('Ticks', [0.1 0.3 0.5 0.7 0.9], 'TickLabels', {'(v)', '(iv)', '(iii)', '(ii)','(i)'});

pos1 = get(sub(1),'Position'); set(sub(1),'Position',[pos1(1)-0.01 pos1(2)-0.01 pos1(3) pos1(4)]); 

annotation(fig1,'textbox',[0.89 0.03 0.055 0.055],'String',{'\it a'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
 
annotation(fig1,'textbox',[0.01 0.96 0.2 0.055],'String',{'\it b'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
set(gca,'FontName','Times New Roman','FontSize',28); 
