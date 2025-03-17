mat_files=["parametric_analysis_cs_bs_v1.mat",... % #1
    "parametric_analysis_cs_ba_v1.mat",... % #2 - used for figures in main manuscript
    "parametric_analysis_cs_ba_v2.mat",... % #3
    "parametric_analysis_ca_bs_v2.mat",... % #4
    "parametric_analysis_ca_ba_v1.mat",... % #5
    "sensitivity_analysis_cs_ba_v1_5yr_avg",... % #6
    "cs_ba_rewrite_v1"]; % #7 Vector of data files to be studied
file_number=2;
dir_loc='mat-data/';
file_name=strcat(dir_loc,mat_files(file_number));
load(file_name); % Loads your file from the array of filenames above

%% Load ColorFiles
colorfile = matfile("colormaps/colormap_auto_benefits2.mat"); 
colorize = colorfile.custommap;
colorizerev = colorfile.custommaprev;
bencolorfile = matfile("colormaps/colormap_auto_benefits.mat"); 
bencolorize = bencolorfile.benmap;
bencolorizerev = bencolorfile.benmaprev;
xscolorfile = matfile("colormaps/colormap_auto_xs.mat"); 
xscolorize = xscolorfile.xsmap;
xscolorizerev = xscolorfile.xsmaprev;
lrcolorfile = matfile("colormaps/auto_colormap_lr.mat"); 
lrcolorize = lrcolorfile.lrcustommap;
lrcolorizerev = lrcolorfile.lrcustommaprev;
eqcolorfile = matfile("colormaps/colormap_auto_equal.mat"); 
eq_color = eqcolorfile.equal_cmap;
eq_color_rev = eqcolorfile.equal_cmap_rev;
opeqcolorfile = matfile("colormaps/colormap_auto_op_eq.mat"); 
opeq_color = opeqcolorfile.op_equal_cmap;
opeq_colorrev = opeqcolorfile.op_equal_cmap_rev;
com1_colorfile=matfile("colormaps/colormap_auto_com1.mat"); 
com1_color= com1_colorfile.com1_cmap;
com1_colorrev= com1_colorfile.com1_cmap_rev;
com2_colorfile=matfile("colormaps/colormap_auto_com2.mat"); 
com2_color= com2_colorfile.com2_cmap;
com2_colorrev= com2_colorfile.com2_cmap_rev;

%% General Plot Variables
% font sizes
fig_font='Arial';
axis_font=16;
title_font=16;
sg_font=18;
gca_font=12;
MarkSize=40;
legendboxes=16;

% labels
xlabel_c2vec='Community 2 Property Values ($1M)';
ylabel_c1vec='Community 1 Property Values ($1M)';

%to make transect lines
xtransect=[C2_vec(1) C2_vec(end)];
ytransect=[C1_vec(1) C1_vec(end)];
xtransect_e6=[C2_vec(1)/1e6 C2_vec(end)/1e6];
ytransect_e6=[C1_vec(1)/1e6 C1_vec(end)/1e6];

%% Figure 1: Behaviors

figure (101) %figure 5a-b in paper
sgtitle('Behaviors','FontSize',sg_font,'fontname',fig_font)
subplot(1,2,1)
pcolor(C1_vec/1e6,C2_vec/1e6,Beh_cord_pa)
colormap(colorizerev)
hold on
p7 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#A42517','MarkerEdgeColor','none')
p6 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#C82D1D','MarkerEdgeColor','none')
p5 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#E46E62','MarkerEdgeColor','none')
p4 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#F9E3E1','MarkerEdgeColor','none')
p3 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#D1D9E8','MarkerEdgeColor','none')
p2 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#7DAAD8','MarkerEdgeColor','none')
p1 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#0459B4','MarkerEdgeColor','none')

legend([p1, p2,p3,p4,p5,p6,p7],{'Seaward Growth','seaward growth / hold the line',...
    'Hold the Line','Mixed','Slow Retreat','Mixed Retreat',...
    'Full Retreat'},Location='southoutside',FontSize=legendboxes,Box='off')
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
title('Behaviors: Coordination','FontSize',title_font,'fontname',fig_font)
clim([0 9])
% c1=colorbar;

subplot(1,2,2)
pcolor(C1_vec/1e6,C2_vec/1e6,Beh_cons_pa)
colormap(colorizerev)
shading flat
hold on
p7 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#A42517','MarkerEdgeColor','none')
p6 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#C82D1D','MarkerEdgeColor','none')
p5 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#E46E62','MarkerEdgeColor','none')
p4 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#F9E3E1','MarkerEdgeColor','none')
p3 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#D1D9E8','MarkerEdgeColor','none')
p2 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#7DAAD8','MarkerEdgeColor','none')
p1 = plot(NaN,NaN,'Marker','s','LineStyle','None','MarkerSize',MarkSize,'MarkerFaceColor','#0459B4','MarkerEdgeColor','none')

legend([p1, p2,p3,p4,p5,p6,p7],{'Seaward Growth','seaward growth / hold the line',...
    'Hold the Line','Mixed','Slow Retreat','Mixed Retreat',...
    'Full Retreat'},Location='southoutside',FontSize=legendboxes,Box='off')
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
title('Behaviors: Non-Coordination','FontSize',title_font,'fontname',fig_font)
clim([0 9])
% c2=colorbar;

%% Figure 2: Nourishment Volumes Relative to Other Community 
max_cord=max(max(abs(V1_nrsh_cord_pa-V2_nrsh_cord_pa)));
max_cons=max(max(abs(V1_nrsh_cons_pa-V2_nrsh_cons_pa)));
diff_max=max(max_cord/1e7, max_cons/1e7);

figure (102) %figure 5a-b in paper
sgtitle({'Difference in Nourishment Volumes per Million Cubic Meters', 'Community 1 - Community 2'},'FontSize',sg_font,'fontname',fig_font)
subplot(1,2,1)
pcolor(C1_vec/1e6,C2_vec/1e6,(V1_nrsh_cord_pa-V2_nrsh_cord_pa)/1e7)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
colormap(opeq_colorrev)
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
title('Coordination','FontSize',title_font,'fontname',fig_font)
% clim([-2.5e7 2.5e7])
clim([-1*diff_max diff_max])
colorbar('southoutside')

subplot(1,2,2)
pcolor(C1_vec/1e6,C2_vec/1e6,(V1_nrsh_cons_pa-V2_nrsh_cons_pa)/1e7)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
colormap(opeq_colorrev)
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',gca_font)
title('Non-Coordination','FontSize',title_font,'fontname',fig_font)
% clim([-2.5e7 2.5e7])
clim([-1*diff_max diff_max])
colorbar('southoutside')

%% Figure 3: Nourishment Volumes: Coordination
max_cord_c1=max(max(abs(V1_nrsh_cord_pa)));
max_cord_c2=max(max(abs(V1_nrsh_cons_pa)));
max_cord=max([max_cord,max_cons]);


figure (103) %figure 5a-b in paper
sgtitle('Nourishment Volumes: Coordination')
subplot(1,2,1)
pcolor(C1_vec/1e6,C2_vec/1e6,V1_nrsh_cord_pa)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
% colormap(opeq_colorrev)
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',sg_font)
title('Community 1','FontSize',title_font,'fontname',fig_font)
clim([0 max_cord])
c1=colorbar;

subplot(1,2,2)
pcolor(C1_vec/1e6,C2_vec/1e6,V2_nrsh_cord_pa)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
% colormap(opeq_color)
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',sg_font)
title('Community 2','FontSize',title_font,'fontname',fig_font)
clim([0 max_cord])
c2=colorbar;

colormap(subplot(1,2,1),com1_colorrev)
colormap(subplot(1,2,2),com2_colorrev)


%% Figure 4: Nourishment Volumes - Conservative Non-Coordination
max_cons_c1=max(max(abs(V1_nrsh_cons_pa)));
max_cons_c2=max(max(abs(V2_nrsh_cons_pa)));
max_cons=max([max_cons_c1,max_cons_c2]);

figure (104) %figure 5a-b in paper
sgtitle('Nourishment Volumes: Conservative Non-Coordination')
subplot(1,2,1)
pcolor(C1_vec/1e6,C2_vec/1e6,V1_nrsh_cons_pa)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
% colormap(opeq_colorrev)
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',sg_font)
title('Community 1','FontSize',title_font,'fontname',fig_font)
clim([0 max_cons])
c1=colorbar;

subplot(1,2,2)
pcolor(C1_vec/1e6,C2_vec/1e6,V2_nrsh_cons_pa)
hold on
plot(xtransect_e6,ytransect_e6,'-k',LineWidth=2)
% colormap(opeq_colorrev)
shading flat
pbaspect([1 1 1])
xlabel(xlabel_c2vec,'FontSize',axis_font)
ylabel(ylabel_c1vec,'FontSize',axis_font)
set(gca,'FontSize',sg_font)
title('Community 2','FontSize',title_font,'fontname',fig_font)
clim([0 max_cons])
c2=colorbar;

colormap(subplot(1,2,1),com1_colorrev)
colormap(subplot(1,2,2),com2_colorrev)


