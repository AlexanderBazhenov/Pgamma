% 2023-01-20
% 2021-09-09
% Pgamma data
% Physics Letters 1992
clear all
close all
% dirroot = 'D:\Data\ST\2021\T\'
dirroot = 'D:\Data\ST\2023\T\'
dirOld =  'd:\Data\ST\2022\T\'
% 2024-03-4
dirroot = 'D:\ST\2024\T\'
dirOld =  'd:\ST\2022\T\'
%
dirroot ='e:\Users\Public\Documents\ST\2023\T\'
dirOld =  'e:\Users\Public\Documents\ST\2022\T\'
% 2022-06-16
% ki
cd(dirroot), pwd

addpath(dirroot)
addpath(dirOld)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phys Lett 1992
Pgamma1992Data;
PgammaComptonData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OxfordBlue = [0, .33, .71]
RoyalMail = 4.58*[0.218, .032, 0.042]
Pantone = 3*[0.128, 0.140, 0.036]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLUMN DATA PRESENTATION
% Scattering Diagram
figure
h=errorbar (1:length(PgammaPh), PgammaPh, PgammaPhStat,".b");
set (h, 'linewidth', 2);
set(h, 'color', OxfordBlue )
## set(h, 'color', RoyalMail )
## set(h, 'color', Pantone )
xlim( [1-1 length(PgammaPh)+1] )

title_str=strcat('\it P_{\gamma} 1992')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
title('')
xlabel('');
ylabel('');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
xticks([1:15])
yticks([-15 -10 -5 0 5])
xticklabels([1:15])
yticklabels([-15 -10 -5 0 5])
figure_name_out=strcat('PgammaPhNoTickLabels', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAND DATA PRESENTATION
% 2022-07-11
figure
nums = 1:numel(PgammaPh);
vals =PgammaPh;
err = PgammaPhStat;
h = errorbar (vals, nums, err, err,">.");
ylim([1-1 length(PgammaPh)+1])
set (h, 'linewidth', 2);

title_str=strcat('\it P_{\gamma} 1992')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
title('')
xlabel('\it \delta, \times 10^{-5}');
ylabel('\it Sample number');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
yticks([1:15])
xticks([-15 -10 -5 0 5])
yticklabels([1:15])
xticklabels([-15 -10 -5 0 5])
figure_name_out=strcat('PgammaPhNoTickLabelsR90', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pkg load interval
PgammaStd = midrad(PgammaPh, PgammaPhStat)

JKPh =jaccardKRSet(PgammaStd)
[oskorbin_center, k] = estimate_uncertainty_center(PgammaStd)

nums = 1:length(PgammaStd);
vals = PgammaPh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OSKORBIN REGULARIZATION
figure
errorbar (vals, nums, k*rad(PgammaStd), k*rad(PgammaStd),">.k");
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(PgammaStd), rad(PgammaStd),">.r");
set(h,"linewidth",2)
set(h,"markersize",12)
title ("Interval measurements of Pgamma");
xlabel("Pgamma")
ylabel("Measurement number")
plot(oskorbin_center*[1 1], ylim, 'm--')
set(gca, 'fontsize', 14);
box('off')
figure_name_out=strcat('PgammaPhOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PgammaStd
max(inf(PgammaStd))
min(sup(PgammaStd))


nums = 1:length(PgammaStd);
vals = PgammaPh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KAUCHER ESTIMATION OF DATA
figure
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(PgammaStd), rad(PgammaStd),">.k");
set(h,"linewidth",2)
set(h,"markersize",12)
title ("Interval measurements of Pgamma");
xlabel("Pgamma")
ylabel("Measurement number")
plot(max(inf(PgammaStd))*[1 1], ylim, 'b--')
plot(min(sup(PgammaStd))*[1 1], ylim, 'b--')
plot(PgammaKxc*[1 1], ylim, 'm--')
title ("Interval measurements of Pgamma");
xlabel("Pgamma")
ylabel("Measurement number")
figure_name_out=strcat('PgammaPhMin', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2021-11-23
% interval mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERVAL MODE
X=PgammaStd
[mode, mu_array, max_mu, mode_ind, c_array, C]= modeIR4(X);

%%%%%%%%%%%%%%%%%%%%%%%% C ARRAY
figure
errorbar (1:length(PgammaPh), PgammaPh, PgammaPhStat,".b");
xlim( [1-1 length(PgammaPh)+1] )
hold on
for ii=1:length(C)
  xx=1:length(X);
  yy=C(ii)*ones(length(X),1);
 plot(xx,yy, '--k')
end
title_str=strcat('\it C array')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it Run number');
ylabel('\it Data');
figure_name_out=strcat('PgammaPhCarray', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%% MODE
figure
errorbar (1:length(PgammaPh), PgammaPh, PgammaPhStat,".b");
xlim( [1-1 length(PgammaPh)+1] )
hold on
mode_array=[inf(mode), sup(mode)];
for ii=1:length(mode_array)
  xx=1:length(X);
  yy=mode_array(ii)*ones(length(X),1);
 plot(xx,yy, '-r')
end
title_str=strcat('\it Mode')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it Run number');
ylabel('\it Data');
figure_name_out=strcat('PgammaPhMode', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODE HIST
% 2021-11-25
figure
stairs(C(1:length(C)-1), mu_array,'-b', 'linewidth', 2)
hold on
% 1st last
xx=[C(1) C(1)]
yy=[0 mu_array(1)]
h = plot(xx,yy, "-b", 'linewidth', 2)
xx=[C(end) C(end)]
yy=[0 mu_array(end)]
plot(xx,yy, "-b", 'linewidth', 2)
        xx=[C(end-1) C(end)]
        yy=[mu_array(end) mu_array(end)]
        plot(xx,yy, '-b', 'linewidth', 2)
ylim([0 max_mu+1])
% mu_max
for ii=1:length(mode_ind)
  cnow=c_array(mode_ind(ii));
        xx=[inf(cnow) sup(cnow)]
        yy=[mu_array(mode_ind(ii)) mu_array(mode_ind(ii))]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_array(mode_ind(ii))]
        plot(xx,yy, '--r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_array(mode_ind(ii)+1)]
        plot(xx,yy, '--r', "linewidth", 1)
end
% 2022-03-23
box('off')
set(gca, 'fontsize', 14);
title_str=strcat('\it ModeFig')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it data');
ylabel('\it \mu');
figure_name_out=strcat('PgammaPhModeFig', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA AND MODE TOGETHER
nums = 1:length(PgammaStd);
vals = PgammaPh;
figure
ylim([0 nums(end)+1])
h = errorbar (vals, nums, rad(PgammaStd), rad(PgammaStd),">.b");
set(h,"linewidth",1)
set(h,"markersize",12)
title ("Interval measurements of Pgamma");
xlabel("Pgamma")
ylabel("Measurement number")
ylim([0 length(PgammaStd)+1])
figure_name_out=strcat('PgammaPhRotated', '.png')
print('-dpng', '-r300', figure_name_out), pwd

figure
subplot(2,1,1)
ylim([0 nums(end)+1])
h = errorbar (vals, nums, rad(PgammaStd), rad(PgammaStd),">.b");
set(h,"linewidth",1)
set(h,"markersize",12)
title ("Interval data and mode");
%xlabel("data")
ylabel('\it Measurement number')
ylim([0 length(PgammaStd)+1])
%
subplot(2,1,2)
stairs(C(1:length(C)-1), mu_array,'-b')
hold on
% 1st last
xx=[C(1) C(1)]
yy=[0 mu_array(1)]
plot(xx,yy, "-b")
xx=[C(end) C(end)]
yy=[0 mu_array(end)]
plot(xx,yy, "-b")
        xx=[C(end-1) C(end)]
        yy=[mu_array(end) mu_array(end)]
        plot(xx,yy, '-b')
ylim([0 max_mu+1])
% mu_max
for ii=1:length(mode_ind)
  cnow=c_array(mode_ind(ii));
        xx=[inf(cnow) sup(cnow)]
        yy=[mu_array(mode_ind(ii)) mu_array(mode_ind(ii))]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_array(mode_ind(ii))]
        plot(xx,yy, ':r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_array(mode_ind(ii)+1)]
        plot(xx,yy, ':r', "linewidth", 1)
end
title_str='' %strcat('\it ModeFig')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it data');
ylabel('\it \mu');
figure_name_out=strcat('PgammaPhDataModeFig', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% /DATA AND MODE TOGETHER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[oskorbin_center, k] = estimate_uncertainty_center(PgammaStd)

nums = 1:length(PgammaStd);
vals = PgammaPh;

figure
errorbar (vals, nums, k*rad(PgammaStd), k*rad(PgammaStd),">.k");
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(PgammaStd), rad(PgammaStd),">.r");
set(h,"linewidth",2)
set(h,"markersize",12)
title ("Interval measurements of Pgamma");
xlabel("Pgamma")
ylabel("Measurement number")
plot(oskorbin_center*[1 1], ylim, 'm--')
figure_name_out=strcat('PgammaPhOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd

PgammaStdCover = midrad(PgammaPh, k*PgammaPhStat)
figure
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(PgammaStdCover), rad(PgammaStdCover),">.k");
set(h,"linewidth",1)
set(h,"markersize",12)
title ("Interval measurements of Pgamma");
xlabel("data")
ylabel("Measurement number")
plot(oskorbin_center*[1 1], ylim, 'm--')
figure_name_out=strcat('PgammaPhOskorbinBlack', '.png')
print('-dpng', '-r300', figure_name_out), pwd

XCover=PgammaStdCover

% interval mode
% 2022-09-09
[modec, mu_arrayc, max_muc, mode_indc, c_arrayc, Cc]= modeIR4(XCover);

figure
stairs(Cc(1:length(Cc)-1), mu_arrayc,'-b')
xlim([-19 9])
hold on
% 1st last
xx=[Cc(1) Cc(1)]
yy=[0 mu_arrayc(1)]
plot(xx,yy, "-b")
xx=[Cc(end) Cc(end)]
yy=[0 mu_arrayc(end)]
plot(xx,yy, "-b")
        xx=[Cc(end-1) Cc(end)]
        yy=[mu_arrayc(end) mu_arrayc(end)]
        plot(xx,yy, '-b')
ylim([0 max_muc+1])
% mu_max
for ii=1:length(mode_indc)
  cnow=c_arrayc(mode_indc(ii));
        xx=[inf(cnow) sup(cnow)]
        yy=[mu_arrayc(mode_indc(ii)) mu_arrayc(mode_indc(ii))]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_arrayc(mode_indc(ii))]
        plot(xx,yy, ':r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_arrayc(mode_indc(ii)+1)]
        plot(xx,yy, ':r', "linewidth", 1)
end
title_str='' %strcat('\it ModeFig')
title_str='Interval data regularizated mode';
title(title_str, 'FontSize', 12, 'Fontweight', 'normal')
xlabel('\it data');
ylabel('\it \mu');
figure_name_out=strcat('PgammaPhOskorbinModeFig', '.png')
print('-dpng', '-r300', figure_name_out), pwd

figure
subplot(2,1,1)
ylim([0 nums(end)+1])
h = errorbar (vals, nums, rad(XCover), rad(XCover),">.b");
set(h,"linewidth",1)
set(h,"markersize",12)
title_str='"Interval data regularizated Oskorbin center';
title(title_str, 'FontSize', 12, 'Fontweight', 'normal')
%xlabel("data")
ylabel('\it Measurement number')
ylim([0 length(PgammaStd)+1])
xlim([-19 9])
%
subplot(2,1,2)
stairs(Cc(1:length(C)-1), mu_arrayc,'-b')
xlim([-19 9])
hold on
% 1st last
xx=[Cc(1) Cc(1)]
yy=[0 mu_arrayc(1)]
plot(xx,yy, "-b")
xx=[Cc(end) Cc(end)]
yy=[0 mu_arrayc(end)]
plot(xx,yy, "-b")
        xx=[Cc(end-1) Cc(end)]
        yy=[mu_arrayc(end) mu_arrayc(end)]
        plot(xx,yy, '-b')
ylim([0 max_muc+1])
% mu_max
for ii=1:length(mode_indc)
  cnow=c_arrayc(mode_indc(ii));
        xx=[inf(cnow) sup(cnow)]
        yy=[mu_arrayc(mode_indc(ii)) mu_arrayc(mode_indc(ii))]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_arrayc(mode_ind(ii))]
        plot(xx,yy, ':r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_arrayc(mode_ind(ii)+1)]
        plot(xx,yy, ':r', "linewidth", 1)
end
title_str='' %strcat('\it ModeFig')
title_str='Interval data regularizated mode';
title(title_str, 'FontSize', 12, 'Fontweight', 'normal')
xlabel('\it data');
ylabel('\it \mu');
figure_name_out=strcat('PgammaPhDataModeFig', '.png')
print('-dpng', '-r300', figure_name_out), pwd

[mode, mu_array, max_mu, mode_ind, c_array, C]= modeIR(PgammaStd);

figure
stairs(C(1:length(C)-1), mu_array,'-k')
hold on
% 1st last
xx=[C(1) C(1)]
yy=[0 mu_array(1)]
plot(xx,yy, "-k")
xx=[C(end) C(end)]
yy=[0 mu_array(end)]
plot(xx,yy, "-k")
        xx=[C(end-1) C(end)]
        yy=[mu_array(end) mu_array(end)]
        plot(xx,yy, '-k')
ylim([0 max_mu+1])
% mu_max
for ii=1:length(mode_ind)
  cnow=c_array(mode_ind(ii));
        xx=[inf(cnow) sup(cnow)]
        yy=[mu_array(mode_ind(ii)) mu_array(mode_ind(ii))]
        plot(xx,yy, '-r', "linewidth", 3)
##        xx=[inf(cnow) sup(cnow)]
##        yy=[0 0]
##        plot(xx,yy, '-r', "linewidth", 3)
##        xx=[inf(cnow) inf(cnow)]
##        yy=[0 mu_array(mode_ind(ii))]
##        plot(xx,yy, ':r', "linewidth", 1)
##        xx=[sup(cnow) sup(cnow)]
##        yy=[0 mu_array(mode_ind(ii)+1)]
##        plot(xx,yy, ':r', "linewidth", 1)
end
xlim([-19 9])
% Oskorbin
stairs(Cc(1:length(C)-1), mu_arrayc,'-b')
xlim([-19 9])
hold on
% 1st last
xx=[Cc(1) Cc(1)]
yy=[0 mu_arrayc(1)]
plot(xx,yy, "-b")
xx=[Cc(end) Cc(end)]
yy=[0 mu_arrayc(end)]
plot(xx,yy, "-b")
        xx=[Cc(end-1) Cc(end)]
        yy=[mu_arrayc(end) mu_arrayc(end)]
        plot(xx,yy, '-b')
ylim([0 max_muc+1])
% mu_max
for ii=1:length(mode_indc)
  cnow=c_arrayc(mode_indc(ii));
        xx=[inf(cnow) sup(cnow)]
        yy=[mu_arrayc(mode_indc(ii))-1 mu_arrayc(mode_indc(ii))]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_arrayc(mode_ind(ii))]
        plot(xx,yy, ':r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_arrayc(mode_ind(ii)+1)]
        plot(xx,yy, ':r', "linewidth", 1)
end

title_str='Mode - Ini and Regularized Data' %strcat('\it ModeFig')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it data');
ylabel('\it \mu');
figure_name_out=strcat('PgammaPhDataModeIniOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2021-12-06
% ADD OLD DATA (1981)
PgammaPhOld=-12
PgammaPhStatOld=4
PgammaPhTot=[PgammaPh; PgammaPhOld]
PgammaPhStatTot=[PgammaPhStat; PgammaPhStatOld]

figure
errorbar (1:length(PgammaPhTot)-1, PgammaPhTot(1:end-1), PgammaPhStatTot(1:end-1),".b");
hold on
errorbar (length(PgammaPhTot), PgammaPhTot(end), PgammaPhStatTot(end),".r");
xlim( [1-1 length(PgammaPhTot)+1] )
title_str=strcat('\it P_{\gamma} all')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it Run number');
ylabel('\it Data');
set(gca, 'fontsize', 14);
figure_name_out=strcat('PgammaPhAll', '.png')
print('-dpng', '-r300', figure_name_out), pwd

pkg load interval
PgammaStdTot = midrad(PgammaPhTot, PgammaPhStatTot)
[oskorbin_center_Tot, k_Tot] = estimate_uncertainty_center(PgammaStdTot)

nums = 1:length(PgammaStdTot);
vals = PgammaPhTot;

figure
errorbar (vals, nums, k_Tot*rad(PgammaStdTot), k_Tot*rad(PgammaStdTot),">.k");
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(PgammaStdTot), rad(PgammaStdTot),">.r");
set(h,"linewidth",2)
set(h,"markersize",12)
title ("Interval measurements of Pgamma");
xlabel("Data")
ylabel("Measurement number")
set(gca, 'fontsize', 14);
plot(oskorbin_center_Tot*[1 1], ylim, 'm--')
figure_name_out=strcat('PgammaPhTotOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interval mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PgammaStdTotCover = midrad(PgammaPhTot, k_Tot*PgammaPhStatTot)
X=PgammaStdTotCover

% interval mode
[modeTot, mu_arrayTot, max_muTot, mode_indTot, c_arrayTot, CTot]= modeIR(X);

% C array
figure
errorbar (1:length(PgammaPhTot), PgammaPhTot, k_Tot*PgammaPhStatTot,".b");
xlim( [1-1 length(PgammaPhTot)+1] )
hold on
for ii=1:length(CTot)
  xx=1:length(X);
  yy=CTot(ii)*ones(length(X),1);
 plot(xx,yy, '--k')
end
title_str=strcat('\it C array')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('\it Run number');
ylabel('\it Data');
set(gca, 'fontsize', 12);
figure_name_out=strcat('PgammaPhTotOskorbinCarray', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% mode
figure
stairs(CTot(1:length(CTot)-1), mu_arrayTot,'-b')
xlim([min(CTot)-1 max(CTot)+1])
hold on
% 1st last
xx=[CTot(1) CTot(1)]
yy=[0 mu_arrayTot(1)]
plot(xx,yy, "-b")
xx=[CTot(end) CTot(end)]
yy=[0 mu_arrayTot(end)]
plot(xx,yy, "-b")
        xx=[CTot(end-1) CTot(end)]
        yy=[mu_arrayTot(end) mu_arrayTot(end)]
        plot(xx,yy, '-b')
ylim([0 max_muTot+1])
% mu_max
for ii=1:length(mode_indTot)
  cnow=c_arrayTot(mode_indTot(ii));
        xx=[inf(cnow) sup(cnow)]
        yy=[mu_arrayTot(mode_indTot(ii)) mu_arrayTot(mode_indTot(ii))]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_arrayTot(mode_indTot(ii))]
        plot(xx,yy, ':r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_arrayTot(mode_indTot(ii)+1)]
        plot(xx,yy, ':r', "linewidth", 1)
end
title_str='' %strcat('\it ModeFig')
title_str='Interval data regularizated mode';
title(title_str, 'FontSize', 12, 'Fontweight', 'normal')
xlabel('\it data');
ylabel('\it \mu');
set(gca, 'fontsize', 12);
figure_name_out=strcat('PgammaPhTotOskorbinModeFig', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% /add old data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variability estimates
% 2022-01-20

PgammaStd
max(inf(PgammaStd))
min(sup(PgammaStd))
IPgammaK=ki(max(inf(PgammaStd)), min(sup(PgammaStd)))
PgammaKxc=( max(inf(PgammaStd)) + min(sup(PgammaStd)) ) /2
IPgammaUni=ki(min(inf(PgammaStd)), max(sup(PgammaStd)))
% I_Uni
mid(IPgammaUni)
% K
K=sum(PgammaStd)/length(PgammaStd)
Pgammac=mid(K)
rad(PgammaStd)

radI=(max(sup(PgammaStd))-min(inf(PgammaStd)))/2
for k=1:length(PgammaStd)
  Delta(k)= distIR (PgammaStd(k), infsup(Pgammac, Pgammac));
end
norm(Delta, 1)
norm(Delta, 2)
norm(Delta, Inf)

max(inf(PgammaStd))-min(sup(PgammaStd))
max(sup(PgammaStd))-min(inf(PgammaStd))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covering Subset
K=sum(PgammaStd)/length(PgammaStd)
Pgammac=mid(K)
Outliers=[]
for ii=1:length(PgammaStd)
  if sup(PgammaStd(ii)) < mid(K)
    Outliers=[Outliers, ii]
  endif
  if inf(PgammaStd(ii)) > mid(K)
    Outliers=[Outliers, ii]
  endif
end
  AllData=1: length(PgammaStd)
 Insiders= setdiff(AllData, Outliers)

 supI = min(sup(PgammaStd(Insiders)))
infI = max(inf(PgammaStd(Insiders)))
I=infsup(infI, supI)
midI=mid(I)
radI=rad(I)

 figure
ylim([0 nums(end)+1])
hold on
h1 = errorbar (vals(Insiders), Insiders, rad(PgammaStd(Insiders)), rad(PgammaStd(Insiders)),">.k");
h2=errorbar (vals(Outliers), Outliers, rad(PgammaStd(Outliers)), rad(PgammaStd(Outliers)),">.r");
set(h1,"linewidth",2)
set(h1,"markersize",12)
set(h2,"linewidth",2)
set(h2,"markersize",12)
plot( supI*[1 1], ylim, 'b--')
plot(infI*[1 1], ylim, 'b--')
%plot(Pgammac*[1 1], ylim, 'm--')
%title ("Interval measurements of Pgamma");
xlabel("Data")
ylabel("Measurement number")
set(gca, 'fontsize', 14);
figure_name_out=strcat('PgammaPhInsidersOultilers', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% 2022-03-05
figure
hold on
%h1 = errorbar (vals(Insiders), Insiders, rad(PgammaStd(Insiders)), rad(PgammaStd(Insiders)),">.k");
h=errorbar (Insiders, PgammaPh(Insiders), PgammaPhStat(Insiders),".b");
set (h, 'linewidth', 2);
xlim( [1-1 length(PgammaPh)+1] )
%set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
xticks([1:15])
yticks([-14 -12 -10 -8 -6 -4 -2 0])
plot( xlim, supI*[1 1], 'k--')
plot( xlim, infI*[1 1], 'k--')
figure_name_out=strcat('PgammaPhCoverNoTickLabels', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-03-22 Example Cover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CoverInd
ExampleCoverInd= [ 1, 8, 15]
supI = min(sup(PgammaStd(ExampleCoverInd)))
infI = max(inf(PgammaStd(ExampleCoverInd)))
I=infsup(infI, supI)
midI=mid(I)
radI=rad(I)

figure
hold on
%h1 = errorbar (vals(Insiders), Insiders, rad(PgammaStd(Insiders)), rad(PgammaStd(Insiders)),">.k");
h=errorbar (ExampleCoverInd, PgammaPh(ExampleCoverInd), PgammaPhStat(ExampleCoverInd),".b");
% 2022-03-23
set (h, 'linewidth', 2);
%set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 1);
box('off')
xticklabels([])
yticklabels([])
xticks(ExampleCoverInd)
yticks([infI midI supI])
% space for I
xlim( [-2 length(PgammaPh)+1] )
% plot I
plot( xlim, supI*[1 1], 'k--')
plot( xlim, midI*[1 1], 'r--')
plot( xlim, infI*[1 1], 'k--')
% information set I
%
xshift = -1
plot( [0 0]+xshift, [infI supI], 'r-', 'linewidth', 2)
% horizontal
plot( [-0.2 0.2]+xshift, supI*[1 1], 'r-', 'linewidth', 2)
plot( [-0.2 0.2]+xshift, infI*[1 1], 'r-', 'linewidth', 2)
plot( 0+xshift, (midI+infI)/2, 'sk' )
plot( 0+xshift, midI, 'ok' )
figure_name_out=strcat('ExampleCover', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-03-23 Example Cover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NonCoverInd
ExampleNonCoverInd= [ 2, 7, 15]
supJ = max(sup(PgammaStd(ExampleNonCoverInd)))
infJ = min(inf(PgammaStd(ExampleNonCoverInd)))
J=infsup(infJ, supJ)
midJ=mid(J)
radJ=rad(J)

figure
hold on
%h1 = errorbar (vals(Insiders), Insiders, rad(PgammaStd(Insiders)), rad(PgammaStd(Insiders)),">.k");
h=errorbar (ExampleNonCoverInd, PgammaPh(ExampleNonCoverInd), PgammaPhStat(ExampleNonCoverInd),".b");
% 2022-03-23
set (h, 'linewidth', 2);
%set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 1);
box('off')
xticklabels([])
yticklabels([])
xticks(ExampleNonCoverInd)
yticks([infJ midJ supJ])
% space for I
xlim( [-2 length(PgammaPh)+1] )
% plot I
plot( xlim, supJ*[1 1], 'k--')
plot( xlim, midJ*[1 1], 'r--')
plot( xlim, infJ*[1 1], 'k--')
% information set I
%
xshift = -1
plot( [0 0]+xshift, [infJ supJ], 'r-', 'linewidth', 2)
% horizontal
plot( [-0.2 0.2]+xshift, supJ*[1 1], 'r-', 'linewidth', 2)
plot( [-0.2 0.2]+xshift, infJ*[1 1], 'r-', 'linewidth', 2)
plot( 0+xshift, (midJ+infJ)/2, 'sk' )
plot( 0+xshift, midJ, 'ok' )
figure_name_out=strcat('ExampleNonCover', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Datanow=PgammaStd(Insiders)

for k=1:length(Datanow)
  datanow = Datanow(k);
  DeltaInsiders(k)= distIR (datanow, infsup(midI, midI));
end
norm(DeltaInsiders, 1)
norm(DeltaInsiders, 2)
norm(DeltaInsiders, Inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     ADD COMPTON DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PgammaComptonData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PGAMMA - COMPTON LINEAR MODEL
x= PgammaComptonStd
y= PgammaStd(1:length(PgammaComptonStd))
X = [ x.^0 x ];

infA = inf(X)
supA = sup(X)
infb = inf(y)
supb = sup(y)

[tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,0)


nums = 1:length(x);
errx = rad(x)
erry = rad(y)

figure
h = errorbar (mid(x), mid(y), errx, errx, erry, erry, "~>")
set (h, "linestyle", "none")
axis('equal')
xlabel("Compton")
ylabel("Photo Peak")
set(gca, 'fontsize', 14);
figure_name_out=strcat('PgammaComptonPhotoPeak', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTS PHOTOPEAK + COMPTON
% 2023-02-11
k21= 1.57
x= k21*PgammaComptonStd
y= PgammaStd(1:length(PgammaComptonStd))
X = [ x.^0 x ];

infA = inf(X)
supA = sup(X)
infb = inf(y)
supb = sup(y)

[tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,0)


nums = 1:length(x);
errx = rad(x)
erry = rad(y)

figure
h = errorbar (mid(x), mid(y), errx, errx, erry, erry, "~>")
set (h, "linestyle", "none")
axis('equal')
xlabel("Compton")
ylabel("Photo Peak")
set(gca, 'fontsize', 14);
figure_name_out=strcat('PgammaComptonPhotoPeakk21', '.png')
print('-dpng', '-r300', figure_name_out), pwd




xp = [ min(inf(x)); max(sup(x)) ]
Xp = [ xp.^0 xp ];
yp=Xp*argmax

%  Pgamma - Compton linear zero
x= PgammaComptonStd
y= PgammaStd(1:length(PgammaComptonStd))
% 2022-04-12
X = [x ];
infA = inf(X)
supA = sup(X)
infb = inf(y)
supb = sup(y)
[tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,0)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-03-30 PHOTOPEAK + COMPTON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2023-01-21
% PAIR DATA SCATTERING
PgammaPairPhCompton
%
min(inf(PgammaStd))
max(sup(PgammaStd))
max(inf(PgammaStd))
min(sup(PgammaStd))
jaccardKRSet(PgammaStd)
%
min(inf(PgammaComptonStd))
max(sup(PgammaComptonStd))
max(inf(PgammaComptonStd))
min(sup(PgammaComptonStd))
jaccardKRSet(PgammaComptonStd)

% 2023-02-11
%%%%%%%%%%%%%%%%%%%%  JK PgammaStd - PgammaComptonStd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JK_PhComp_array = jaccard2(PgammaStd(1:numel(PgammaComptonStd)), PgammaComptonStd)
k21 = 1.57
JK_PhComp_arrayk21 = jaccard2(PgammaStd(1:numel(PgammaComptonStd)), k21*PgammaComptonStd)
%
filename = "JK_PhComp.txt";
fid = fopen (filename, "w");
for step = 1: numel(PgammaComptonStd)
input1 = JK_PhComp_array(step)
input2 = JK_PhComp_arrayk21(step)
XYstr = LateXTableStr (step, input1, input2)
  fprintf(fid, XYstr);
end
fclose (fid);
%%%%%%%%%%%%%%%%%%%%  JK PgammaStd - PgammaComptonStd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[modeP, mu_arrayP, max_muP, mode_indP, c_arrayP, CP]= modeIR4(PgammaStd)
[modeC, mu_arrayC, max_muC, mode_indC, c_arrayC, CC]= modeIR4(PgammaComptonStd)

% 2022-12-20
% Ini Data
k21=1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCATTERING DIAGRAM CORRECTED
figure
hold on
h1=errorbar (1:length(PgammaPh), PgammaPh, PgammaPhStat,".b");
set (h1, 'linewidth', 2);
h2=errorbar ([1:length(PgammaComptonStd)]+0.5, k21*PgammaCompton, k21*PgammaComptonStat,".r");
set (h2, 'linewidth', 2);
xlim( [1-1 length(PgammaPh)+1] )

title_str=strcat('\it P_{\gamma} 1992')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
title('')
xlabel('');
ylabel('');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
yticks([-15 -10 -5 0])
yticklabels([-15 -10 -5 0])
xticks([1 12 15])
xticklabels([1 12 15])
figure_name_out=strcat('PgammaPhComptonIniNoTickLabels', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Olga Skrekel
% Pb Background
BgIntensity = [ 3.3336 2.186 1.618]
BgIntensityR = BgIntensity /BgIntensity(3)
%k21= BgIntensity(1)/BgIntensity(2)

% 2022-04-08
% paper equation
k21= (1 +  BgIntensity(1)/BgIntensity(3) )/ (1 +  BgIntensity(2)/BgIntensity(3) )
k21=  BgIntensity(1)/BgIntensity(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCATTERING DIAGRAM CORRECTED
figure
hold on
h1=errorbar (1:length(PgammaPh), PgammaPh, PgammaPhStat,".b");
set (h1, 'linewidth', 2);
h2=errorbar ([1:length(PgammaComptonStd)]+0.5, k21*PgammaCompton, k21*PgammaComptonStat,".r");
set (h2, 'linewidth', 2);
xlim( [1-1 length(PgammaPh)+1] )

title_str=strcat('\it P_{\gamma} 1992')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
title('')
xlabel('');
ylabel('');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
yticks([-15 -10 -5 0])
yticklabels([-15 -10 -5 0])
xticks([1 12 15])
xticklabels([1 12 15])
figure_name_out=strcat('PgammaPhComptonCorrectedNoTickLabels', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INTERVAL MODE PHOTOPEAK + COMPTON
PgammaTotini = [ PgammaStd; PgammaComptonStd]
JKTotini  =jaccardKRSet(PgammaTotini)
modePini

% k21=1.575
PgammaTot = [ PgammaStd; k21*PgammaComptonStd]
JKTot  =jaccardKRSet(PgammaTot)

[mode, mu_array, max_mu, mode_ind, c_array, C]= modeIR4(PgammaTot);
% 2023-01-24
PgammaModeXY
% /INTERVAL MODE PHOTOPEAK + COMPTON

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODE HIST
% 2021-11-25
figure
%stairs(C(1:length(C)-1), mu_array,'-b', "linewidth", 2)
stairs(C(1:length(C)), [mu_array mu_array(end)],'-b', "linewidth", 2)
hold on
% 1st last
xx=[C(1) C(1)]
yy=[0 mu_array(1)]
plot(xx,yy, "-b", "linewidth", 2)
xx=[C(end) C(end)]
yy=[0 mu_array(end)]
plot(xx,yy, "-b",  "linewidth", 2)
        xx=[C(end-1) C(end)]
        yy=[mu_array(end) mu_array(end)]
        plot(xx,yy, '-b',  "linewidth", 2)
ylim([0 20])
xlim([min(C)-0.5 max(C)+0.5])
% mu_max
for ii=1:length(mode_ind)
  cnow=c_array(mode_ind(ii));
##        xx=[inf(cnow) sup(cnow)]
##        yy=[mu_array(mode_ind(ii)) mu_array(mode_ind(ii))]
##        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_array(mode_ind(ii))]
        plot(xx,yy, ':r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_array(mode_ind(ii)+1)]
        plot(xx,yy, ':r', "linewidth", 1)
end
title_str=strcat('')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('');
ylabel('');
set(gca, 'linewidth', 2);
set(gca, 'FontSize', 14);
box('off')
xticklabels([])
yticklabels([])
xticks([-15 -10 -5 0 5])
yticks([ 5 10 15 20])
xticklabels([-15 -10 -5 0 5])
yticklabels([5 10 15 20])
figure_name_out=strcat('PgammaPhModeTotFig', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OSKORBIN_CENTER PHOTOPEAK + COMPTON

[oskorbin_center_Total, k_Total] = estimate_uncertainty_center(PgammaTot )

nums = 1:length(PgammaTot);
vals = mid(PgammaTot);

figure
errorbar (vals, nums,   k_Total*rad(PgammaTot), k_Total*rad(PgammaTot),">.k");
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(PgammaTot), rad(PgammaTot),">.r");
set(h,"linewidth",2)
set(h,"markersize",12)
set(gca, 'linewidth', 2);
set(gca, 'FontSize', 14);
box('off')
plot(oskorbin_center_Total*[1 1], ylim, 'm--')
figure_name_out=strcat('PgammaPhTotalOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 2022-04-08
% Indexes for optimization
k21=BgIntensity(1)/BgIntensity(2)
kcorrCompt=1.575
PgammaTot = [ PgammaStd; k21*PgammaComptonStd];
JKTot  =jaccardKRSet(PgammaTot)
[mode, mu_array, max_mu, mode_ind, c_array, C]= modeIR4(PgammaTot);
max_mu
[oskorbin_center_k, k_now] = estimate_uncertainty_center(PgammaTot);
k_now

JK_array = []
max_mu_array = []
k_array = []
kcorrCompt_array=1:0.025:2
for jj = 1:length(kcorrCompt_array)
  kc_now=kcorrCompt_array(jj)
  X = [ PgammaStd; kc_now*PgammaComptonStd];
  [JK_now, max_mu_now, k_now] = SetMeasures(X);
  JK_array = [JK_array JK_now];
  max_mu_array = [max_mu_array max_mu_now];
  k_array = [k_array k_now];
end
% save optimization arrays
save PgammaArrays.mat JK_array max_mu_array k_array kcorrCompt_array
% load PgammaArrays
% good optimization interval
max_mu_ind = find (max_mu_array >= max(max_mu_array))
Oskorbin_ind = find (k_array <= min(k_array)+0.0001)
max_mu_Oskorbin_ind = intersect(max_mu_ind, Oskorbin_ind)
ind_out = setdiff(1:length(kcorrCompt_array), max_mu_Oskorbin_ind )
% out of good
[indmax iinmax2] = max(diff(ind_out))
ind_out_L = ind_out(1:iinmax2)
ind_out_R = setdiff( ind_out, ind_out_L )
%
kcorrCompt_array(max_mu_Oskorbin_ind )
%
[max_JK_array, max_JK_array_ind] = max(JK_array)
[min_JK_array, min_JK_array_ind] = min(JK_array)
k_array_opt =kcorrCompt_array(max_JK_array_ind)
max_mu_array_opt=max_mu_array(max_JK_array_ind)

% 2022-07-08
figure
hold on
h1 = plot(kcorrCompt_array, JK_array, '-r', "linewidth", 2)
xx = [k_array_opt k_array_opt]
yy = [ min_JK_array max_JK_array ]
h2 = plot(xx, yy, '--k', "linewidth", 2)
box('off')
set(gca, 'linewidth', 2);
set(gca, 'FontSize', 14);
ylim([min_JK_array -0.195])
xticks([1 1.2 1.4 1.6 1.8 2.0])
yticks([-0.26 -0.24 -0.22 -0.20])
xticklabels([1 1.2 1.4 1.6 1.8 2.0])
yticklabels([-0.26 -0.24 -0.22 -0.20 ])
xlabel('\it k_{12}')
ylabel('\it JK')
figure_name_out=strcat('PgammaPhTotalJK', '.png')
print('-dpng', '-r300', figure_name_out), pwd


figure
subplot(2,1,1)
plot(kcorrCompt_array, JK_array, '-k', "linewidth", 2)
box('off')
set(gca, 'linewidth', 2);
set(gca, 'FontSize', 14);
ylim([-0.28 -0.195])
xticks([1 1.2 1.4 1.6 1.8 2.0])
yticks([-0.26 -0.24 -0.22 -0.20])
xticklabels([1 1.2 1.4 1.6 1.8 2.0])
yticklabels([-0.26 -0.24 -0.22 -0.20 ])
subplot(2,1,2)
plot(kcorrCompt_array, k_array, '-k', "linewidth", 2)
box('off')
set(gca, 'linewidth', 2);
set(gca, 'FontSize', 14);
box('off')
xticklabels([])
yticklabels([])
xticks([1 1.2 1.4 1.6 1.8 2.0])
yticks([1.7 1.8 1.9 2.0 ])
xticklabels([1 1.2 1.4 1.6 1.8 2.0])
yticklabels([1.7 1.8 1.9 2.0 ])
figure_name_out=strcat('PgammaPhTotalJKOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP HERE !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-10-10
% 2022-10-11
%
% PlotJKOskorbinMu

%%%%%%%%%%%%%%%%%%%%%%%%  /PLOT 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2022-04-09
k21=BgIntensity(1)/BgIntensity(2)
X = [ PgammaStd; k21*PgammaComptonStd];
[JK_now, max_mu_now, k_now] = SetMeasures(X);

max(inf(X))
min(sup(X))

figure
stairs(C(1:length(C)), [mu_array mu_array(end)],'-b', "linewidth", 2)
hold on
% 1st last
xx=[C(1) C(1)]
yy=[0 mu_array(1)]
plot(xx,yy, "-b", "linewidth", 2)
xx=[C(end) C(end)]
yy=[0 mu_array(end)]
plot(xx,yy, "-b",  "linewidth", 2)
##        xx=[C(end-1) C(end)]
##        yy=[mu_array(end) mu_array(end)]
##        plot(xx,yy, '-b', , "linewidth", 2)
ylim([0 20])
xlim([min(C)-0.5 max(C)+0.5])
% mu_max
for ii=1:length(mode_ind)
  cnow=c_array(mode_ind(ii));
##        xx=[inf(cnow) sup(cnow)]
##        yy=[mu_array(mode_ind(ii)) mu_array(mode_ind(ii))]
##        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) sup(cnow)]
        yy=[0 0]
        plot(xx,yy, '-r', "linewidth", 3)
        xx=[inf(cnow) inf(cnow)]
        yy=[0 mu_array(mode_ind(ii))]
        plot(xx,yy, '--r', "linewidth", 1)
        xx=[sup(cnow) sup(cnow)]
        yy=[0 mu_array(mode_ind(ii)+1)]
        plot(xx,yy, '--r', "linewidth", 1)
end
title_str=strcat('')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
xlabel('');
ylabel('');
set(gca, 'linewidth', 1);
set(gca, 'FontSize', 14);
box('off')
xticklabels([])
yticklabels([])
xticks([-15 -10 -5 0 5])
yticks([ 5 10 15 20])
xticklabels([-15 -10 -5 0 5])
yticklabels([5 10 15 20])
figure_name_out=strcat('PgammaPhTotalMuArray', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% Scattering Diagram
figure
hold on
h1=errorbar (1:length(PgammaPh), PgammaPh, PgammaPhStat,".b");
set (h1, 'linewidth', 2);
h2=errorbar ([1:length(PgammaComptonStd)]+0.5, k21*PgammaCompton, k21*PgammaComptonStat,".r");
set (h2, 'linewidth', 2);
xlim( [1-1 length(PgammaPh)+1] )
title('')
xlabel('');
ylabel('');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
yticks([-15 -10 -5 0])
yticklabels([-15 -10 -5 0])
xticks([1 12 15])
xticklabels([1 12 15])
figure_name_out=strcat('PgammaPhComptonCorrectedk21NoTickLabels', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd

vals=mid(X);
nums=1:length(vals);
figure
errorbar (vals, nums, k_now*rad(X), k_now*rad(X),">.k");
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(X), rad(X),">.r");
set(h,"linewidth",2)
set(h,"markersize",12)
set(gca, 'fontsize', 14);
plot(oskorbin_center_k*[1 1], ylim, 'm--')
figure_name_out=strcat('PgammaPhTotOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd

OnesX=ones(length(X),1)
PgammaX= mid (sum(X./rad(X)) / sum(OnesX./rad(X)))
StdPgammaX=std(mid(X))/sqrt(length(X)-1)

% 2022-04-09
x= k21*PgammaComptonStd
y= PgammaStd(1:length(PgammaComptonStd))
X = [ x.^0 x ];
infA = inf(X)
supA = sup(X)
infb = inf(y)
supb = sup(y)
[tolmax,argmax,envs,ccode] = tolsolvty(infA,supA,infb,supb,0)

figure
[xsort, xind]=sort(mid(x))
ysort = mid(y(xind))
%[ysort, yind]=sort(mid(y))
plot(xsort, ysort, 'sb')

dirwork = 'D:\Data\ST\2022\T\octave-interval'
dirwork ='e:\Users\Public\Documents\ST\2022\T\octave-interval'
cd(dirwork)
pwd
addpath(genpath('./m'))
% cd(dirroot)
% cd octave-interval/
% rmpath(genpath('./m'))
% abstract error
epsilon= 5*ones(length(x),1)
% 2022-04-12
epsilon= rad(y)
% manual
Kwid=2.32
epsilon= Kwid*rad(y(xind))

lb = [0];
irp_steam2 = ir_problem(xsort, ysort, epsilon, lb);
## ������� ������������ ������ ��������� ������ y = beta2 * x
b_int2 = ir_outer(irp_steam2)          # ������������ ������
b_rad2 = 0.5 * (b_int2(2) - b_int2(1)) # ������ ������������ ������
b_mid2 = mean(b_int2)                  # ����� ��������� ��� �������� ������ ��������� beta2
## ������ ��������� beta2 ������� ���������� ���������
b_lsm2 = xsort \ ysort


## ������� ���������� ������������ ��� ������ y = beta2 * x
xlimits = [-13 1];
ir_plotmodelset(irp_steam2,xlimits)
grid off
hold on
ir_scatter(irp_steam2,'bo')
ir_scatter(ir_problem(xp,yp2mid,yp2rad),'r.')
title('')
xlabel('');
ylabel('');
%set (h, "linestyle", "none")
%axis('equal')
xlim([-15 5])
ylim([-25 10])
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
yticks([-15 -10 -5 0])
yticklabels([-15 -10 -5 0])
xticks([-15 -10 -5 0])
xticklabels([-15 -10 -5 0])
figure_name_out=strcat('PgammaPhEqCompton','ExrtaWid',num2str(Kwid), '.png')
print('-dpng', '-r300', figure_name_out), pwd




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-04-19
% BG
% Indexes for optimization
k23=0.12
k23=0.55
PgammaCorrBG=PgammaStd*(1+k23)-k23*PgammaBGStd;
JKTot  =jaccardKRSet(PgammaCorrBG)
[mode, mu_array, max_mu, mode_ind, c_array, C]= modeIR2(PgammaCorrBG);
max_mu
[oskorbin_center_k, k_now] = estimate_uncertainty_center(PgammaCorrBG);
k_now

JK_array = []
max_mu_array = []
k_array = []
k23_array=0:0.02:0.7
for k23_now=k23_array
  X=PgammaStd*(1+k23_now)-k23_now*PgammaBGStd;
  [JK_now, max_mu_now, k_now] = SetMeasures(X);
  JK_array = [JK_array JK_now];
  max_mu_array = [max_mu_array max_mu_now];
  k_array = [k_array k_now];
end

figure
subplot(2,1,1)
plot(k23_array, JK_array, '-k', "linewidth", 2)
box('off')
set(gca, 'linewidth', 2);
set(gca, 'FontSize', 14);
ylim([-0.28 -0.195])
xticks([1 1.2 1.4 1.6 1.8 2.0])
yticks([-0.26 -0.24 -0.22 -0.20])
xticklabels([1 1.2 1.4 1.6 1.8 2.0])
yticklabels([-0.26 -0.24 -0.22 -0.20 ])
subplot(2,1,2)
plot(k23_array, k_array, '-k', "linewidth", 2)
box('off')
set(gca, 'linewidth', 2);
set(gca, 'FontSize', 14);
box('off')
xticklabels([])
yticklabels([])
xticks([1 1.2 1.4 1.6 1.8 2.0])
yticks([1.7 1.8 1.9 2.0 ])
xticklabels([1 1.2 1.4 1.6 1.8 2.0])
yticklabels([1.7 1.8 1.9 2.0 ])
figure_name_out=strcat('PgammaPhBGcorrJKOskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% Scattering Diagram
figure
hold on
h1=errorbar (1:length(PgammaPh), PgammaPh, PgammaPhStat,".b");
set (h1, 'linewidth', 2);
h2=errorbar ([1:length(PgammaPh)]+0.5, mid(PgammaCorrBG), rad(PgammaCorrBG),".r");
set (h2, 'linewidth', 2);
%
xlim( [1-1 length(PgammaPh)+1] )
ylim( [min(inf(PgammaCorrBG)) max(sup(PgammaCorrBG))] )
plot(xlim, max(inf(PgammaCorrBG))*[1 1], 'b--')
plot(xlim, min(sup(PgammaCorrBG))*[1 1], 'b--')
title('')
xlabel('');
ylabel('');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
yticks([-15 -10 -5 0])
yticklabels([-15 -10 -5 0])
xticks([1 12 15])
xticklabels([1 12 15])
figure_name_out=strcat('PgammaPhBGCorrectedk23NoTickLabels', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd

figure
plot(C(1:end-1),mu_array)


PgammaStd
max(inf(PgammaCorrBG))
min(sup(PgammaCorrBG))
IPgammaK=ki(max(inf(PgammaStd)), min(sup(PgammaStd)))
PgammaKxc=( max(inf(PgammaStd)) + min(sup(PgammaStd)) ) /2
IPgammaUni=ki(min(inf(PgammaStd)), max(sup(PgammaStd)))
% I_Uni
mid(IPgammaUni)
% K
K=sum(PgammaStd)/length(PgammaStd)
mid(K)


nums = 1:length(PgammaStd);
vals = PgammaPh;

figure
ylim([0 nums(end)+1])
hold on
h = errorbar (vals, nums, rad(PgammaStd), rad(PgammaStd),">.k");
set(h,"linewidth",2)
set(h,"markersize",12)
title ("Interval measurements of Pgamma");
xlabel("Pgamma")
ylabel("Measurement number")
plot(max(inf(PgammaStd))*[1 1], ylim, 'b--')
plot(min(sup(PgammaStd))*[1 1], ylim, 'b--')
plot(PgammaKxc*[1 1], ylim, 'm--')
title ("Interval measurements of Pgamma");
xlabel("Pgamma")
ylabel("Measurement number")
figure_name_out=strcat('PgammaPhMin', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




PgammaPhComptonBg = [PgammaPh; k21*PgammaCompton]
PgammaPhComptonStatBg = [PgammaPhStat; k21*PgammaComptonStat]
OnesPgammaPhComptonBg = ones(length(PgammaPhComptonBg),1)


PgammaPhComptonMeanW = sum(PgammaPhComptonBg./PgammaPhComptonStatBg)/sum(OnesPgammaPhComptonBg./PgammaPhComptonStatBg)

w  = OnesPgammaPhComptonBg./PgammaPhComptonStatBg

S = std(PgammaPhComptonBg)
sqrt(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-06-16
% ki

addpath(genpath('./kinterval-0.0.1'))

Pgammaki = ki(inf(PgammaStd), sup(PgammaStd))
PBGgammaki = ki(inf(PgammaBGStd), sup(PgammaBGStd))
BGRatio = 0.05
Pgammakinew = Pgammaki + BGRatio*opp(PBGgammaki)

Pgammanew = infsup(inf(Pgammakinew), sup(Pgammakinew))

jaccardKRSet(Pgammanew )
jaccardKRSet(PgammaStd)
[mode, mu_array, max_mu, mode_ind, c_array, C]= modeIR3(Pgammanew );

X = Pgammanew
[JK, max_mu, k_now] = SetMeasures(X)

num = 1:length(PgammaStd)

% Scattering Diagram
figure
h1=errorbar (num+0.5, mid(Pgammanew), rad(Pgammanew),".r");
hold on
h2=errorbar (num, mid(PgammaStd), rad(PgammaStd),".b");

set (h1, 'linewidth', 2);
xlim( [1-1 length(PgammaStd)+1] )

title_str=strcat('\it P_{\gamma} 1992')
title_str=strcat('RBG=', num2str(BGRatio), ' JK=', num2str(jaccardKRSet(Pgammanew )))

title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
title('')
xlabel('');
ylabel('');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
xticks([1:15])
yticks([-15 -10 -5 0 5])
xticklabels([1:15])
yticklabels([-15 -10 -5 0 5])
figure_name_out=strcat('PgammaKR',' RBG=', num2str(BGRatio), '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd


JK_array = []
max_mu_array = []
k_array = []
k23_array=0:0.02:0.2
for k23_now=k23_array
  Pgammaki = ki(inf(PgammaStd), sup(PgammaStd));
PBGgammaki = ki(inf(PgammaBGStd), sup(PgammaBGStd));
Pgammakinew = Pgammaki + k23_now*opp(PBGgammaki);
X = infsup(inf(Pgammakinew), sup(Pgammakinew));

%  X=PgammaStd*(1+k23_now)-k23_now*PgammaBGStd;
  [JK_now, max_mu_now, k_now] = SetMeasures(X);
  JK_array = [JK_array JK_now];
  max_mu_array = [max_mu_array max_mu_now];
  k_array = [k_array k_now];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-05

X = PgammaStd

[mode, mu_array, max_mu, mode_ind, c_array, C, multi]= modeIR4(X);
retval = ModePlot (mode, mu_array, max_mu, mode_ind, c_array, C, multi)
ylim([0 8.5])
xlim([-15 5])
box('off')
xlabel('\it \delta, \times 10^{-5}');
figure_name_out=strcat('PgammaMuArrayIni', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c in X
[CX, c_array, C]= cinX(X)

figure
spy(CX')
xlabel('c-array')
ylabel('X set')
set(gca, 'Fontsize', 14)
xlim([0.5 length(c_array)+0.5])
ylim([0.5 length(X)+0.5])
grid on
title('C in X array')
figure_name_out=strcat('C in X array', '.png')
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mu_array
mu_array = sum(CX, 2)
% mode ind
[mode_mult, ~] = max(mu_array)
% mode
mode = c_array(mode_ind)
%
max_clique_array=[]
max_clique_union=[]

for ii=1:length(mode_ind)
max_clique_now = find( CX(mode_ind(ii),:) > 0);
max_clique_array = [max_clique_array; max_clique_now]
max_clique_union =  union(max_clique_union, max_clique_now)
end
outside_data = setdiff( 1:length(X), max_clique_union)
CXmode  =CX(mode_ind,:)
sumCXmode= sum(CXmode, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grandX
[Y, I, ~] = unique(CX, "rows", "stable")
XonC_array = sum(Y, 1)

[X_mult, X_ind] = max(XonC_array)
grandX = X(X_ind)

m =  length(c_array)
n = length(X)
JK_array = JK_CX(X)
dist_array = dist_CX(X)

% Scattering Diagram
if length(mode_ind) == 1
  retval = CliquePlot (X, mode, max_clique_union)
end
%retval = MultiCliquePlot(X, mode, max_clique_array)
% horizontal
retval = MultiCliquePlot(X, CX, mode, mode_ind, max_clique_array, max_clique_union)
% plot legend
  xlim( [1-1 length(X)+3] )
  legendElement = infsup(1,2)
  Xlegend = [legendElement-8  legendElement-6  legendElement-4  legendElement-2 legendElement ]
  for jj=2:length(Xlegend )
  h = errorbar (length(X)+2, mid(Xlegend(jj)), rad(Xlegend(jj)),".r" );
  set (h, 'linewidth', -0.5+ round(jj/2) );
    text(length(X)+3, mid(Xlegend(jj)), num2str(jj-1), 'fontsize', 14)
end
jj=1
  h = errorbar (length(X)+2, mid(Xlegend(jj)), rad(Xlegend(jj)),".b" );
  set (h, 'linewidth', 0.5 );
  text(length(X)+3, mid(Xlegend(jj)), num2str(0), 'fontsize', 14)
 text(length(X)+2, 4, '\it XC(X)', 'fontsize', 14)
% / plot legend
xticks([1:15] )
xticklabels([1:15])
yticks([-15 -10 -5 0 5])
yticklabels([-15 -10 -5 0 5])
figure_name_out=strcat('PgammaPhMaxCliqueMulti', '.png')
%

% vertical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-16
ScatteringClique
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 1);
xlim([-15.5 10])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure_name_out=strcat('PgammaPhMaxCliqueMulti90', '.png')
figure_name_out=strcat('PgammaPhComptonMaxCliqueMulti90', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-07
% Combined Data

k21 = 1.575
X = [ PgammaStd; k21*PgammaComptonStd];

% Scattering Diagram


% 2022-07-12
[mode, mu_array, max_mu, mode_ind, c_array, C, multi]= modeIR4(X);
retval = ModePlot (mode, mu_array, max_mu, mode_ind, c_array, C, multi)
box('off')
xticklabels([])
yticklabels([])
yticks([0 5 10 15 18])
xticks([-15 -10 -5 0 5])
yticklabels([0 5 10 15 18])
xticklabels([-15 -10 -5 0 5])
xxlim = xlim
xx= [ xxlim(1)  xxlim(2) ]
yy = [max_mu max_mu]
h = plot(xx, yy, '--k')
figure_name_out=strcat('PgammaPhComptonMode', '.png')
print('-dpng', '-r300', figure_name_out), pwd

[CX, c_array, C]= cinX(X)
size(CX)
%
% 2023-01-23
figure
pcolor(CX')
set(gca, 'fontsize', 14)
ax = gca
figure_name_out=strcat('PgammaPhComptonCXmatrix', '.png')
print('-dpng', '-r300', figure_name_out), pwd

figure
spy(CX)
set(gca, 'fontsize', 14)
ax = gca
figure_name_out=strcat('PgammaPhComptonCXmatrix', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%
max_clique_array=[]
max_clique_union=[]

for ii=1:length(mode_ind)
max_clique_now = find( CX(mode_ind(ii),:) > 0);
max_clique_array = [max_clique_array; max_clique_now]
max_clique_union =  union(max_clique_union, max_clique_now)
end
outside_data = setdiff( 1:length(X), max_clique_union)
CXmode  =CX(mode_ind,:)
sumCXmode= sum(CXmode, 1)

if length(mode_ind) == 1
  retval = CliquePlot (X, mode, max_clique_union)
end
set(gca, 'fontsize', 14)
figure_name_out=strcat('PgammaPhComptonMaxClique', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023-01-25
PhNum =1:15
CompNum=16:27
%
figure
hold on
%
h1 = errorbar (PhNum , mid(X(PhNum )), rad(X(PhNum )),".b");
set (h1, 'linewidth', 2);
h2 = errorbar (CompNum-numel(PhNum)+0.2, mid(X(CompNum)), rad(X(CompNum)),".r");
set (h2, 'linewidth', 2);
xlim( [1-1 length(PhNum)+1] )
set(gca, 'fontsize', 14)
figure_name_out=strcat('PgammaPhComptonCorrPairs', '.png')
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2022-07-13
% Compatibility graphics_toolkit

part=[1 2 3 4]
PgammaStdpart = PgammaStd(part)
PgammaPhStatpart=PgammaPhStat(part)
PgammaStdpart'

% Subset
xlimarray=[0.1 3.9]

figure
subplot(1, 2, 1)
h=errorbar (PgammaPh(part), 1:length(PgammaStdpart),  PgammaPhStat(part),">.b");
set(h,"linewidth",2)
set(h,"markersize",14)
ylim( [min(part)-0.5 max(part)+0.5 ])
xlim( [-14 6])
xlabel('\it \delta, \times 10^{-5}');
ylabel('\it Sample number');
set(gca,"linewidth",1)
set(gca,"Fontsize",14)

subplot(1, 2, 2)
hold on
% vertices
xx = [ -4 -9 -4 2   ]
yy = [ 1 2.5 4 2.5  ]
plot(xx, yy, 'sk')
text(xx(1), yy(1)-0.25, '1', 'Fontsize',14, 'HorizontalAlignment', 'center')
text(xx(3), yy(3)+0.25, '3', 'Fontsize',14, 'HorizontalAlignment', 'center')
text(xx(2)-1, yy(2), '2', 'Fontsize',14, 'HorizontalAlignment', 'center')
text(xx(4)+1, yy(4), '4', 'Fontsize',14, 'HorizontalAlignment', 'center')
%
box('off')
axis('off')
ylim( [min(part)-0.5 max(part)+0.5 ])
xlim( [-10 2])
% plot chords
Xp = [ xx(1) xx(2) ], Yp = [ yy(1) yy(2) ], plot( Xp, Yp, '-b'  )
Xp = [ xx(1) xx(3) ], Yp = [ yy(1) yy(3) ], plot( Xp, Yp, '-b'  )
Xp = [ xx(1) xx(4) ], Yp = [ yy(1) yy(4) ], plot( Xp, Yp, '-b'  )
Xp = [ xx(2) xx(3) ], Yp = [ yy(2) yy(3) ], plot( Xp, Yp, '-b'  )
Xp = [ xx(2) xx(4) ], Yp = [ yy(2) yy(4) ], plot( Xp, Yp, '-b'  )


figure_name_out=strcat('GraphPgammaPh1234', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2022-07-19
% Jaccard 2 sets
X = PgammaStd(1:12)
Y = PgammaComptonStd
XY = [X' Y']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
% Metrics
[RelRes, HighLev, Vee_array, Wedge_array, JK_array] = MetricsXY (X, Y);
% /Metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RAW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(4, 1, 1)
hold on
h1 =errorbar (1:length(X), mid(X), rad(X),".b");
xlim([0.2 12.8])
nums = 1:length(X)
nums = nums + 0.2
h2 =errorbar (nums, mid(Y), rad(Y),".r");
xx = [1-0.5 12+0.5]
yy = [minXY minXY]
plot(xx, yy, '--k')
yy = [maxXY maxXY]
plot(xx, yy, '--k')
xticks([1:12])
box off
ylabel('Data')
%
subplot(4, 1, 2)
plot(JK_array)
ylabel('JK Photopeak-Compton')
xlim([0.2 12.8])
xticks([1:12])
box off
%
subplot(4, 1, 3)
hold on
plot(abs(Vee_array))
plot(abs(Wedge_array))
xlim([0.2 12.8])
xticks([1:12])
box off
ylabel('\vee, \wedge Photopeak-Compton')
%
subplot(4, 1, 4)
plot(sumCXmode(1:12))
xlim([0.2 12.8])
ylabel('Sample in Clique')
xlabel('Data sample')
xticks([1:12])
box off
%
figure_name_out=strcat('PhotopeakComptonMetrics', '.png')
print('-dpng', '-r300', figure_name_out), pwd


figure
subplot(3, 1, 1)
hold on
h1 =errorbar (1:length(X), mid(X), rad(X),".b");
xlim([0.2 12.8])
nums = 1:length(X)
nums = nums + 0.2
h2 =errorbar (nums, mid(Y), rad(Y),".r");
xx = [1-0.5 12+0.5]
yy = [minXY minXY]
plot(xx, yy, '--k')
yy = [maxXY maxXY]
plot(xx, yy, '--k')
xticks([1:12])
box off
ylabel('Data')
%
subplot(3, 1, 2)
plot(RelRes)
ylabel('Rel Res')
xlim([0.2 12.8])
xticks([1:12])
box off
%
subplot(3, 1, 3)
plot(-HighLev)
ylabel('High Leverage')
xlim([0.2 12.8])
xticks([1:12])
box off
%
figure_name_out=strcat('PhotopeakComptonMetrics2', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-20
% scattering diagram

X = PgammaStd
Y = PgammaComptonStd
[minX, maxX] = wedgeSet(X);
midX = ( maxX + minX ) / 2;
radX = ( maxX - minX ) / 2;
% Metrics
for ii=1:length(X)
  datanow = X(ii);
  RelResX(ii) = ( mid(datanow ) -  midX ) / rad(datanow );
end
for ii=1:length(X)
  datanow = X(ii);
  HighLevX(ii) = radX  / rad(datanow );
end
[minY, maxY] = wedgeSet(Y);
midY = ( maxY + minY ) / 2;
radY = ( maxY- minY ) / 2;
for ii=1:length(Y)
  datanow = Y(ii);
  RelResY(ii) = ( mid(datanow ) -  midY ) / rad(datanow );
end
for ii=1:length(Y)
  datanow = Y(ii);
  HighLevY(ii) = radY  / rad(datanow );
end
figure
% text(0.55, 0, 'Internal', 'Fontsize', 14);
% text(1.3, 0, 'External', 'Fontsize', 14);
figure_name_out=strcat('PhotopeakStatusDiag', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% Initial data
X = PgammaStd
Y = PgammaComptonStd

XY = [X' Y']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
% Metrics
%[RelRes, HighLev, Vee_array, Wedge_array, JK_array] = MetricsXY (X, Y);
figure
hold on
xxlim = [0.5 1.5]
% plot ouliers area
       vertices = [0 -1 ; xxlim(2) (-xxlim(2)-1); 0 (-xxlim(2)-1); 0 -1 ];
        PP = patch (vertices(:,1), vertices(:,2), 0.95*[1 1 1]);
        set(PP, 'edgecolor', [ 0 0 0]); set(PP, 'facecolor', 0.95*[1 1 0]);
        vertices = [0 1 ; xxlim(2) (xxlim(2)+1); 0 (xxlim(2)+1); 0 1];
        PP = patch (vertices(:,1), vertices(:,2), 0.95*[1 1 1]);
        set(PP, 'edgecolor', [ 0 0 0]); set(PP, 'facecolor', 0.95*[1 1 0]);
% /plot ouliers area
% text(0.6, 2.2, 'Ouliers', 'Fontsize', 14);
% text(0.6, -2.2, 'Ouliers', 'Fontsize', 14);
% plot external area
       vertices = [0 -1 ; 1 0; 0 1; 0 -1 ];
        PP = patch (vertices(:,1), vertices(:,2), [0 0 0]);
        set(PP, 'edgecolor', [ 0 0 0]); set(PP, 'facecolor', 0.95*[0 1 0]);
% /plot external area
% text(0.55, 0, 'Internal', 'Fontsize', 14);
% text(1.2, 0, 'External', 'Fontsize', 14);
% Samples
for ii = 1: length(X)
h = plot(-HighLevX(ii), RelResX(ii), 'ob');
text(-HighLevX(ii)+0.02, RelResX(ii)-0.04, num2str(ii), 'Fontsize', 14);
end
for ii = 1:length(Y)
h = plot(-HighLevY(ii), RelResY(ii), 'or');
text(-HighLevY(ii)+0.02, RelResY(ii)-0.04, num2str(ii), 'Fontsize', 14);
end
% /Samples

% Lines
xx = [1 1]
yy = [-2.5 2.5]
plot(xx, yy, '--k')
xx = [0 1]
yy = [1 0]
plot(xx, yy, '--k')
xx = [0 1]
yy = [-1 0]
plot(xx, yy, '--k')
xx = [0 1.5]
yy = [1 2.5]
plot(xx, yy, '--k')
xx = [0 1.5]
yy = [-1 -2.5]
plot(xx, yy, '--k')
% /Lines

xlabel('high leverage')
ylabel('relative residual')
set(gca, 'Fontsize', 14)
set(gca, 'Linewidth', 1)
ylim( [- 2.5 2.5] )
xlim([0.3 1.4])
xx = [0.3 0.3]
yy = [-2.5 2.5]
plot(xx, yy, '-k')
box('off')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k21 = 1.57
X = PgammaStd
Y = k21*PgammaComptonStd
XY = [X' Y']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
% Metrics
[RelResXY, HighLevXY] = InfluencePlotXYtot (X, Y)
% plot ouliers area
% plot external area
figure
hold on
PLOT_Influence_Diagram_Empty
% Samples
for ii = 1: length(XY)
h = plot(-HighLevXY(ii), RelResXY(ii), 'ob');
text(-HighLevXY(ii)+0.02, RelResXY(ii)-0.04, num2str(ii), 'Fontsize', 14);
end
% Lines
xlim([0.25 1.4])
ylim( [- 2.5 2.5] )
xx = [0.25 0.25]
yy = [-2.5 2.5]
plot(xx, yy, '-k')

figure_name_out=strcat('PhotopeakComptonStatusDiag', '.png')
figure_name_out=strcat('PhotopeakK21ComptonStatusDiag', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = PgammaStd
Y = k21*PgammaComptonStd
PartOskorbin = 0.9
kOskorbin = 1.75 * PartOskorbin
XX = midrad(mid(X), kOskorbin*rad(X))
YY = midrad(mid(Y), kOskorbin*rad(Y))
XY = [XX' YY']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
[RelResXY, HighLevXY] = InfluencePlotXYtot (XX, YY)
%
PLOT_Influence_Diagram_Empty
%
for ii = 1: length(XY)
h = plot(-HighLevXY(ii), RelResXY(ii), 'ob');
%text(-HighLevXY(ii)+0.005, RelResXY(ii)-0.02, num2str(ii), 'Fontsize', 14);
end
%
title_str = strcat('k=', num2str(PartOskorbin ))
title(title_str, 'Fontsize', 14)
%
xlim([0.025 0.2])
ylim( [- 1.5 1.5] )
xx = [0.025 0.025]
yy = [- 1.5 1.5]
plot(xx, yy, '-k')
box('off')
figure_name_out=strcat('PhotopeakK21ComptonStatusDiag', 'PartOskorbin=',num2str(PartOskorbin ), '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-07-23
for PartOskorbin = 0.6 : 0.1 : 2
% PartOskorbin = 1 % 1.2 % 1.1 %1.0 %0.9
kOskorbin = 1.75 * PartOskorbin
XX = midrad(mid(X), kOskorbin*rad(X))
YY = midrad(mid(Y), kOskorbin*rad(Y))
XY = [XX' YY']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
[RelResXY, HighLevXY] = InfluencePlotXYtot (XX, YY)
%
PLOT_Influence_Diagram_Empty
%
for ii = 1: length(XY)
h = plot( abs(HighLevXY(ii)), RelResXY(ii), 'ob');
%text(-HighLevXY(ii)+0.005, RelResXY(ii)-0.02, num2str(ii), 'Fontsize', 14);
end
%
title_str = strcat('k=', num2str(PartOskorbin,'%.1f'))
title(title_str, 'Fontsize', 14)
%
box('off')
xlim([0 1.4])
ylim( [- 2.5 2.5] )
%
xlabel(' | high leverage | ')
ylabel('relative residual')
set(gca, 'Fontsize', 14)
set(gca, 'Linewidth', 1)
%
figure_name_out=strcat('PhotopeakK21ComptonStatusDiag', 'PartOskorbin=',num2str(PartOskorbin,'%.1f'), '.png')
print('-dpng', '-r300', figure_name_out), pwd
end
xlim([0.025 0.2])
ylim( [- 1.5 1.5] )
xx = [0.025 0.025]
yy = [- 1.5 1.5]
plot(xx, yy, '-k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load video
 PartOskorbin = 0.6 : 0.1 : 2
  nframes = length(PartOskorbin);

 fn = fullfile (dirroot, "GreenSpace.mp4");
 w = VideoWriter (fn);
 open (w);


 for jj = 1:nframes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PartOskorbin = 1 % 1.2 % 1.1 %1.0 %0.9
kOskorbin = 1.75 * PartOskorbin(jj)
XX = midrad(mid(X), kOskorbin*rad(X))
YY = midrad(mid(Y), kOskorbin*rad(Y))
XY = [XX' YY']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
[RelResXY, HighLevXY] = InfluencePlotXYtot (XX, YY)
%
figure
hold on
PLOT_Influence_Diagram_Empty
%
for ii = 1: length(XY)
h = plot( abs(HighLevXY(ii)), RelResXY(ii), 'ob');
%text(-HighLevXY(ii)+0.005, RelResXY(ii)-0.02, num2str(ii), 'Fontsize', 14);
end
%
title_str = strcat('k=', num2str(PartOskorbin(jj),'%.1f'))
title(title_str, 'Fontsize', 10)
%
box('off')
xlim([0 1.4])
ylim( [- 2.5 2.5] )
%
xlabel(' | high leverage | ')
ylabel('relative residual')
set(gca, 'Fontsize', 10)
set(gca, 'Linewidth', 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   writeVideo (w, getframe (gcf));
 close all
 endfor
 close (w)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 2022-08-29
 % rad Ipsilon < 0

 k21 = 1.57
X = PgammaStd
Y = k21*PgammaComptonStd
XY = [X' Y']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
% Metrics
[RelResXY, HighLevXY] = InfluencePlotXYtot (X, Y)
% plot ouliers area
% plot external area
xxlim = [-1.5 1.5]

figure
hold on
PLOT_Influence_Diagram_Empty_Double
% Samples
for ii = 1: length(XY)
h = plot(HighLevXY(ii), RelResXY(ii), 'ob');
text(HighLevXY(ii)+0.02, RelResXY(ii)-0.04, num2str(ii), 'Fontsize', 14);
end
% Lines
xlim([-1.5 1.4])
ylim( [- 2.5 2.5] )
figure_name_out=strcat('PhotopeakComptonStatusDiag1', '.png')
figure_name_out=strcat('PhotopeakK21ComptonStatusDiag1', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%% video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 PartOskorbin = 0.6 : 0.1 : 2
  nframes = length(PartOskorbin);

 fn = fullfile (dirroot, "GreenSpace.mp4");
 w = VideoWriter (fn);
 open (w);


 for jj = 1:nframes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PartOskorbin = 1 % 1.2 % 1.1 %1.0 %0.9
kOskorbin = 1.75 * PartOskorbin(jj)
XX = midrad(mid(X), kOskorbin*rad(X))
YY = midrad(mid(Y), kOskorbin*rad(Y))
XY = [XX' YY']
[minXY, maxXY] = wedgeSet(XY)
midXY = ( maxXY + minXY ) / 2
radXY = ( maxXY - minXY ) / 2
[RelResXY, HighLevXY] = InfluencePlotXYtot (XX, YY)
%
figure
hold on
PLOT_Influence_Diagram_Empty_Double
%
for ii = 1: length(XY)
h = plot( HighLevXY(ii), RelResXY(ii), 'ob');
%text(-HighLevXY(ii)+0.005, RelResXY(ii)-0.02, num2str(ii), 'Fontsize', 14);
end
%
title_str = strcat('k=', num2str(PartOskorbin(jj),'%.1f'))
title(title_str, 'Fontsize', 10)
%
box('off')
xlim([-1.5 1.4])
ylim( [- 2.5 2.5] )
%
xlabel(' high leverage  ')
ylabel('relative residual')
set(gca, 'Fontsize', 10)
set(gca, 'Linewidth', 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   writeVideo (w, getframe (gcf));
 close all
 endfor
 close (w)
%%%%%%%%%%%%%%%% /video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2022-09-08

Colors_special
figure
hold on
% vertices
xx = [ -4 -9 -4 2   ]
yy = [ 1 2.5 4 2.5  ]
h = plot(xx, yy, 'sk')
set(h, 'color', OxfordBlue), set(h, 'Markersize', 10), set(h, 'linewidth', 2)
text(xx(1), yy(1)-0.5, '1', 'Fontsize',14, 'HorizontalAlignment', 'center')
text(xx(3), yy(3)+0.5, '3', 'Fontsize',14, 'HorizontalAlignment', 'center')
text(xx(2)-0.5, yy(2), '2', 'Fontsize',14, 'HorizontalAlignment', 'center')
text(xx(4)+0.5, yy(4), '4', 'Fontsize',14, 'HorizontalAlignment', 'center')
%
box('off')
axis('off')
xlim( [-10 2])
ylim([ 0 5])
% plot chords
Xp = [ xx(1) xx(2) ], Yp = [ yy(1) yy(2) ], h=plot( Xp, Yp, '-b'  ),
set(h, 'color', OxfordBlue), set(h, 'linewidth', 2)
Xp = [ xx(1) xx(3) ], Yp = [ yy(1) yy(3) ], h=plot( Xp, Yp, '-b'  )
set(h, 'color', OxfordBlue), set(h, 'linewidth', 2)
%Xp = [ xx(1) xx(4) ], Yp = [ yy(1) yy(4) ], plot( Xp, Yp, '-b'  )
Xp = [ xx(2) xx(3) ], Yp = [ yy(2) yy(3) ], h=plot( Xp, Yp, '-b'  )
set(h, 'color', OxfordBlue), set(h, 'linewidth', 2)
%Xp = [ xx(2) xx(4) ], Yp = [ yy(2) yy(4) ], plot( Xp, Yp, '-b'  )


figure_name_out=strcat('Graph1234', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scattering diagram
% 2022-09-08
figure
nums = 1:numel(PgammaPh);
vals =PgammaPh;
err = PgammaPhStat;
h = errorbar (vals, nums, err, err,">.");
ylim([1-1 length(PgammaPh)+1])
set (h, 'linewidth', 2);
set (h, 'color', OxfordBlue);
%
X = PgammaStd;
hold on
for jj=1:numel(X)
xx = [ inf(X(jj)) inf(X(jj)) ];
yy = [ 0 jj ];
h = plot(xx, yy, '--k');
xx = [ sup(X(jj)) sup(X(jj)) ];
h = plot(xx, yy, '--k');
end
%
title_str=strcat('\it P_{\gamma} 1992')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
title('')
xlabel('\it \delta, \times 10^{-5}');
ylabel('\it Sample number');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 2);
box('off')
xticklabels([])
yticklabels([])
yticks([1:15])
xticks([-15 -10 -5 0 5])
yticklabels([1:15])
xticklabels([-15 -10 -5 0 5])
figure_name_out=strcat('PgammaPhEdgesR90', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% /Scattering diagram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Scattering diagram Total
% 2022-09-09
figure
nums = 1:numel(X);
vals =mid(X);
err = rad(X);
h = errorbar (vals, nums, err, err,">.");
ylim([1-1 length(X)+1])
set (h, 'linewidth', 2);
set (h, 'color', OxfordBlue);
%
hold on
for jj=1:numel(X)
xx = [ inf(X(jj)) inf(X(jj)) ];
yy = [ 0 jj ];
h = plot(xx, yy, '--k');
set(h, 'linewidth', 0.5);
xx = [ sup(X(jj)) sup(X(jj)) ];
h = plot(xx, yy, '--k');
set(h, 'linewidth', 0.5);
end
%
title_str=strcat('\it P_{\gamma} 1992')
title(title_str, 'FontSize', 14, 'Fontweight', 'normal')
title('')
xlabel('\it \delta, \times 10^{-5}');
ylabel('\it Sample number');
set (h, "linestyle", "none")
%axis('equal')
set(gca, 'fontsize', 14);
set(gca, 'linewidth', 1);
box('off')
xticklabels([])
yticklabels([])
yticks([1:numel(X)])
xticks([-15 -10 -5 0 5])
yticklabels([1:numel(X)])
xticklabels([-15 -10 -5 0 5])
figure_name_out=strcat('PgammaPhComptonEdgesR90', '.png')
%figure_name_out=strcat('PgammaPh', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% /Scattering diagram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-09-18
% refine hist
X = PgammaStd'
minC12 = min( inf(X) ), maxC12 = max( sup(X) ), midC12 = ( min( inf(X) )+ max( sup(X) ) )/2,
clear C,
C = [  infsup(minC12, midC12), infsup(midC12, maxC12)  ]

Cnow = C
HistVals = histInt (X, Cnow)
ylabel('\it \mu')
xlabel('\it \delta, \times 10^{-5}');
figure_name_out=strcat('PgammaPhHistIntHalf', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% refine
while (length(Cnow) < length(X)) & (max(HistVals)>0)
[Cedges, Cint ] = Xrefined2 (Cnow);
Cnow = Cint;
HistVals = histInt (X, Cnow);
HistVals
figure_name_out=strcat('PgammaPhHistBins=', num2str(length(Cint)), '.png')
print('-dpng', '-r300', figure_name_out), pwd
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-09-19
X = PgammaStd'
% engrosse hist
[Cedges, Cint ] = Xengrossed2 (X)
Cnow = Cint
HistVals = histInt (X, Cnow)
while (length(Cnow) > 2)
[Cedges, Cint ] = Xengrossed2 (Cnow)
Cnow = Cint
HistVals = histInt (X, Cint)
figure_name_out=strcat('PgammaPhHistEngrossBins=', num2str(length(Cint)), '.png')
print('-dpng', '-r300', figure_name_out), pwd
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
