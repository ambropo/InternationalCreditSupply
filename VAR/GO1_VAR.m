%% VAR
%==========================================================================
% Estimates a panel VAR using specifed data
%==========================================================================
clear all; clear session; close all; clc
warning off all

%% 1. USER DEFINED SETTINGS
%==========================================================================
% Load settings
outdir = 'OUT/';
frequency  = 'Q';
filename = {'ALL_Log'};
eval( ['load ' outdir 'Data' frequency '_' filename{1}]);

% Save settings
writeout   = 1;
highquality = 1;

% Choose endogenous settings
shock_name = 'LEV';
Hvnames_long = {'Leverage';'Cross-border Credit';'Consumption';'House Price';'Real Exch. Rate';'Current Account'};
Hvnames      = {'LEV';'KF';'CONS';'RHP';'RXUSD';'CA';};

% VAR settings
Hnvar = length(Hvnames);
shock = find(strcmp(shock_name,Hvnames));
user_VARfo = '1985Q1'; 
user_VARlo = '2012Q4'; 
XoXlag = 0;


%% 2. VAR ESTIMATION
%==========================================================================
% Load IV data (for robustness)
[IVdata, IVtext] = xlsread('DataQ_57C.xls','IV');
IVnames = IVtext(1,2:end);
for ii=1:length(IVnames)
    IV.(IVnames{ii}) = IVdata(:,ii);
end

% VAR specification
nlag = 1; det = 3; nsteps = 40; ident='oir'; VARtol=.999; NOBStol=40;
VARfo = find(strcmp(user_VARfo,dates)); VARlo = find(strcmp(user_VARlo,dates)); 

% Initialize
IR_all = nan(nsteps,Hnvar,ncountry);
IR_all_norescale = nan(nsteps,Hnvar,ncountry);
VD_all = nan(nsteps,Hnvar,ncountry);
eigen   = nan(ncountry,1);
exclude = nan(ncountry,1);
exclude_text = cell(ncountry,1);
RESID = nan(VARlo-VARfo+1,ncountry,Hnvar);

for jj=1:ncountry

    % Create matrices of variables for the VAR
    ENDO.(cnames{jj}) = nan(VARlo-VARfo+1,Hnvar);
    if exist('EXOnvar','var'); EXOG.cnames{jj} = nan(VARlo-VARfo+1,EXOnvar); end
    
    % Make ENDO.(cnames{jj})
    for ii=1:Hnvar
        ENDO.(cnames{jj})(:,ii) = DATA.(Hvnames{ii})(VARfo:VARlo,jj);
    end
    
    % Make EXO
    if exist('EXOnvar','var')
        for ii=1:EXOnvar
            EXOG.cnames{jj}(:,ii) = IV.(EXOvnames{ii})(VARfo:VARlo);
        end
    end
    
    % Define a common sample between all series 
    if ~exist('EXOnvar','var')
        [TEMP, aux_fo, aux_lo] = CommonSample(ENDO.(cnames{jj}));
        ENDO.(cnames{jj}) = TEMP(:,1:Hnvar);
    else
        [TEMP, aux_fo, aux_lo] = CommonSample([ENDO.(cnames{jj}) EXOG.cnames{jj}]);
        ENDO.(cnames{jj}) = TEMP(:,1:Hnvar);
        EXOG.cnames{jj} = TEMP(:,Hnvar+1:end);
    end
    clear TEMP
    
    % Estimate the VAR
    if length(ENDO.(cnames{jj}))<NOBStol || strcmp('USA',cnames(jj))
        VAR.(cnames{jj}) = NaN;
        IR.(cnames{jj}) = NaN;
        VD.(cnames{jj}) = NaN;
        eigen(jj,1) = NaN;
        exclude(jj) = 1;
        exclude_text(jj) = {'nobs'};
    else
        % Estimate VAR
        if ~exist('EXOnvar','var')
            [VAR.(cnames{jj}), VARopt.(cnames{jj})] = VARmodel(ENDO.(cnames{jj}),nlag,det);
        else
            [VAR.(cnames{jj}), VARopt.(cnames{jj})] = VARmodel(ENDO.(cnames{jj}),nlag,det,EXOG.cnames{jj});
        end
        % Store residuals
        for ii=1:Hnvar
            RESID(nlag+aux_fo+1:end-aux_lo,jj,ii) = VAR.(cnames{jj}).residuals(:,ii);
        end
        % If unstable don't compute IRFs
        if VAR.(cnames{jj}).maxEig>VARtol
            IR.(cnames{jj}) = NaN;
            VD.(cnames{jj}) = NaN;
            eigen(jj,1) = VAR.(cnames{jj}).maxEig;
            exclude(jj) = 1;
            exclude_text(jj) = {'eigen'};
        else
            VARopt.(cnames{jj}).nsteps = nsteps;
            VARopt.(cnames{jj}).pick = 1;
            [IR.(cnames{jj}), VAR.(cnames{jj})] = VARir(VAR.(cnames{jj}),VARopt.(cnames{jj}));
            [VD.(cnames{jj}), VAR.(cnames{jj})] = VARfevd(VAR.(cnames{jj}),VARopt.(cnames{jj}));
            eigen(jj,1) = VAR.(cnames{jj}).maxEig;
        end
    end
end
cross_section = sum(eigen<VARtol);

%% 3. RE-ORGANIZE IRFS (AND CENSOR) 
%==========================================================================
% Shock size: 1 std deviation
shock_size = nan(ncountry,1);
for ii=1:ncountry
    if ~isnan(IR.(cnames{ii}))
        aux = sqrt(diag(VAR.(cnames{ii}).sigma));
        shock_size(ii,1) = aux(shock);
    end
end
rescale = 100*nanmean(shock_size);
% Re-organize
for ii=1:ncountry
    if ~isnan(IR.(cnames{ii}))
        IR_all(:,:,ii) = rescale*IR.(cnames{ii})(:,:,shock)./IR.(cnames{ii})(1,shock,shock);
        IR_all_norescale(:,:,ii) = IR.(cnames{ii})(:,:,shock);
        VD_all(:,:,ii) = VD.(cnames{ii})(:,shock,:);
    end
end
% Censor IR
IRcens = prctile(IR_all,[97.5 2.5],3);
IR_all_cens = IR_all;
for jj=1:Hnvar
    for ii=1:nsteps
        IR_all_cens(ii,jj,IR_all(ii,jj,:)>IRcens(ii,jj,1)) = NaN;
        IR_all_cens(ii,jj,IR_all(ii,jj,:)<IRcens(ii,jj,2)) = NaN;
    end
end
% Censor VD
VDcens = prctile(VD_all,[97.5 2.5],3);
VD_all_cens = VD_all;
for jj=1:Hnvar
    for ii=1:nsteps
        VD_all_cens(ii,jj,VD_all(ii,jj,:)>VDcens(ii,jj,1)) = NaN;
        VD_all_cens(ii,jj,VD_all(ii,jj,:)<VDcens(ii,jj,2)) = NaN;
    end
end

%% 4. COMPUTE MG ESTIMATOR ON IR AND VD
%==========================================================================
% Compute IR mean group estimator
IR_uncens = nanmean(IR_all,3); % not used
N         = ncountry - sum(isnan(IR_all_cens),3);
IR_cens   = nansum(IR_all_cens,3)./N;
IR_std    = sqrt(nanvar(IR_all_cens,3)./(N-1));
IR_upp1   = IR_cens + 1*IR_std;
IR_low1   = IR_cens - 1*IR_std;
IR_upp2   = IR_cens + 2*IR_std;
IR_low2   = IR_cens - 2*IR_std;
% Compute VD mean group estimator
VD_uncens = nanmean(VD_all,3); % not used
N         = ncountry - sum(isnan(VD_all_cens),3);
VD_cens   = nansum(VD_all_cens,3)./N;
VD_std    = sqrt(nanvar(VD_all_cens,3)./(N-1));
VD_upp1   = VD_cens + 1*VD_std;
VD_low1   = VD_cens - 1*VD_std;
VD_upp2   = VD_cens + 2*VD_std;
VD_low2   = VD_cens - 2*VD_std;


%% 5. PLOT MG ESTIMATES (IR and VD)
%==========================================================================
% Plot IR panel estimator
FigSize(24,12)
row = 2; col = 3;
pnames = Hvnames_long;
for ii=shock:Hnvar
    subplot(row,col,ii-shock+1) 
    PlotSwathe(IR_cens(:,ii), [IR_upp2(:,ii) IR_low2(:,ii)], rgb('grey')); hold on;
    PlotSwathe(IR_cens(:,ii), [IR_upp1(:,ii) IR_low1(:,ii)], rgb('very dark blue')); hold on;
    plot(zeros(nsteps),'-k','LineWidth',0.5); hold on;
    title(pnames{ii})
    if ii==1||ii==col+1; ylabel('Percent'); end
    if ii>col; xlabel('Quarters'); end
    xlimits = [0 40]; set(gca,'xLim',xlimits,'xTick', [5 10 15 20 25 30 35 40]);
    axis tight
    set(gca,'Layer','top')
    axis tight
end
SaveFigure([outdir 'IRmg_' filename{1} '_' Hvnames{1}],highquality)
clf('reset')

% Plot VD panel estimator
FigSize(24,12)
pnames = Hvnames_long;
for ii=shock:Hnvar
    subplot(row,col,ii-shock+1) 
    PlotSwathe(VD_cens(:,ii), [VD_upp2(:,ii) VD_low2(:,ii)], rgb('grey')); hold on;
    PlotSwathe(VD_cens(:,ii), [VD_upp1(:,ii) VD_low1(:,ii)], rgb('very dark blue')); hold on;
    plot(zeros(nsteps),'-k','LineWidth',0.5); hold on;
    % plot(VD_uncens(:,ii),':','Color',rgb('dark red'),'LineWidth',2)
    title(pnames{ii})
    if ii==1||ii==col+1; ylabel('Percent'); end
    if ii>col; xlabel('Quarters'); end
    xlimits = [0 40]; set(gca,'xLim',xlimits,'xTick', [5 10 15 20 25 30 35 40]);
    set(gca,'Layer','top')
    axis tight
end
SaveFigure([outdir 'VDmg_' filename{1} '_' Hvnames{1}],highquality)
clf('reset')


%% 7. EXPORT TO EXCEL STATS (IR and VD)
%==========================================================================
% Export list of excluded countries
writecell([cnames exclude_text],[outdir 'OUT_' filename{1} '_' Hvnames{1} '.xls'],'Sheet','excluded')
% List of variabels to be exported
select = {'CONS','RHP','RXUSD'};
% Find position of selected in SOEvnames
for ii=1:length(select)
    select_num(ii) = find(strcmp(Hvnames,select{ii}));
end
% Find ER (to change its sign)
ERsel = find(strcmp(select,'RXUSD'));
% Create and poulate IRFstats matrix
IRFstats.imp = nan(ncountry,length(select));
IRFstats.av4 = nan(ncountry,length(select));
IRFstats.av8 = nan(ncountry,length(select));
IRFstats.max = nan(ncountry,length(select));
for ii=1:ncountry
    temp = IR_all_cens(1,select_num,ii);            IRFstats.imp(ii,:) = temp; 
    temp = nanmean(IR_all_cens(1:4,select_num,ii)); IRFstats.av4(ii,:) = temp; 
    temp = nanmean(IR_all_cens(1:8,select_num,ii)); IRFstats.av8(ii,:) = temp; 
    temp = nanmax(IR_all_cens(1:16,select_num,ii)); IRFstats.max(ii,:) = temp; IRFstats.max(ii,ERsel) = nanmin(IR_all_cens(1:16,select_num(ERsel),ii));
end
% Create and poulate VDstats matrix
VDstats.imp = nan(ncountry,length(select));
VDstats.av4 = nan(ncountry,length(select));
VDstats.av8 = nan(ncountry,length(select));
VDstats.max = nan(ncountry,length(select));
for ii=1:ncountry
    temp = VD_all_cens(1,select_num,ii); VDstats.imp(ii,:) = temp;
    temp = nanmean(VD_all_cens(1:4,select_num,ii)); VDstats.av4(ii,:) = temp;
    temp = nanmean(VD_all_cens(1:8,select_num,ii)); VDstats.av8(ii,:) = temp;
    temp = nanmax(VD_all_cens(1:16,select_num,ii)); VDstats.max(ii,:) = temp;
end

%% 8. EXPORT TO EXCEL LEV RESIDUALS AND PLOT CHART
%==========================================================================
writecell(TabPrint(RESID(:,:,shock),cnames',dates(VARfo:VARlo),5),[outdir 'OUT_' filename{1} '_' Hvnames{1} '.xls'],'Sheet','LEVres')
% Plot shock residuals;
FigSize(18,6)
plot(100.*RESID(:,:,shock),'LineWidth',0.5,'Color',cmap(2)); hold on;
plot(100.*nanmean(RESID(:,:,shock),2),'LineWidth',2,'Color',cmap(1)); hold on;
plot(rescale.*ones(length(RESID(:,:,shock))),':','LineWidth',0.5,'Color',cmap(1)); hold on;
plot(-rescale.*ones(length(RESID(:,:,shock))),':','LineWidth',0.5,'Color',cmap(1)); hold on;
plot(zeros(length(RESID(:,:,shock))),'-k','LineWidth',0.5); hold on;
axis tight; set(gca,'Layer','top')
ylabel('Percent'); 
DatesPlot(Date2Num({user_VARfo}),length(RESID(:,:,shock)),8)
SaveFigure([outdir 'RES_' filename{1} '_' Hvnames{1}],highquality)
clf('reset')
LEVshock = nanmean(RESID(:,:,shock),2);

%% 10. SAVE
%==========================================================================
save( [outdir 'VAR_' filename{1} '_' Hvnames{1}] );
close all

