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

%- Markers for parameter of other illustrated examples -
aa=0.8; 
bb=0.03; % Fig. 3 (type (i))
hold on; plot (aa, bb,'*','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[1 .4 .4]);
bb=0.005;% Fig. 6  (type (ii))
hold on; plot (aa, bb,'d','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[1 .4 .4]);
bb=0.015;% Fig. 9  (type (iii))
hold on; plot (aa, bb,'s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[1 .4 .4]);
bb=0.035;% Fig. 12  (type (iv))
hold on; plot (aa, bb,'v','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[1 .4 .4]);
bb=0.07; % Fig. 14  (type (v))
hold on; plot (aa, bb,'x','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[1 .4 .4]);

%- Wada boundaries -
r=5;
[kiek1,kiek2]=size(Masyvas_diverg);
M=Masyvas_diverg*4;
Mask2(1:kiek1,1:kiek2)=M;
L=5;
spalvu_sk=r;

for iii = 1:kiek1-L+1
    for jjj = 1:kiek2-L+1
        istart = iii; iend = iii + L - 1;
        jstart = jjj; jend = jjj + L - 1;
        mat = M(istart:iend,jstart:jend);
        colors = 0;
        for i = 1:spalvu_sk
            k = i - 1;
            prob = nnz(~(mat-k))/(L*L);
            if prob>0, colors = colors + 1; end
        end
        if colors<3 
            Mask2(istart+(L-1)/2,jstart+(L-1)/2)=5;
        end
    end
end

fig2=figure('Position',[1000 1500 800 500],'Color',[1 1 1]); %[left bottom width height]; 
sub=subplot(1,1,1); 
h_a = 0.001;
h_b = 0.0005;
imagesc(0.2:h_a:1,  0:h_b:0.14, abs(Mask2((1+(L-1)/2):(end-(L-1)/2),(1+(L-1)/2):(end-(L-1)/2))'));   
set(gca,'YDir','normal');

cmap2 = [parula(5); 1 1 1];
colormap(cmap2);
colorbar('Ticks', [0.1 0.3 0.5 0.7 0.9]*4, 'TickLabels', {'(v)', '(iv)', '(iii)', '(ii)','(i)'});

pos1 = get(sub(1),'Position'); set(sub(1),'Position',[pos1(1)-0.01 pos1(2)-0.01 pos1(3) pos1(4)]); %[left buttom ]
 
annotation(fig2,'textbox',[0.89 0.03 0.055 0.055],'String',{'\it a'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
 
annotation(fig2,'textbox',[0.01 0.96 0.2 0.055],'String',{'\it b'},'FontName',...
            'Times New Roman','FontSize',40,'FitBoxToText','off','LineStyle','none');
set(gca,'FontName','Times New Roman','FontSize',28); % axis on; 

