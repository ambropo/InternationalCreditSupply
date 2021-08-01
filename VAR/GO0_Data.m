%% GO0_Data
%==========================================================================
% Creates data set according to the definitions in CodeVar.xls
%==========================================================================
clear all; clear session; close all; clc
warning off all


%% 1. Parameters specified by the user
%--------------------------------------------------------------------------
% lag for first differences
lag = 1; 
lambda = 1600;
frequency = 'Q';  % (Q or Y)
% User provided first observation (fo) and last observation (lo). For
% yearly data use Q4
user_fo = '1970Q1';
user_lo = '2012Q4';


%% 2. Import variables info
%--------------------------------------------------------------------------
outdir = 'OUT/'; mkdir(outdir);

% Import data to initialize remaining parameters
[xlsdata, xlstext] = xlsread('CodeVar.xls','Variable');
vnames_long = xlstext(2:end,1);
vnames      = xlstext(2:end,2);
[nvar, ~]   = size(vnames);
vtreat      = xlsdata(:,1);
clear xlsdata xlstext 
if max(vtreat)==0
    cyc_type = {''};
elseif max(vtreat)<=1
    cyc_type = {'Log'};
elseif max(vtreat)<=3
    cyc_type = {'Diff'};
elseif max(vtreat)<=5
    cyc_type = {'Trend'};
elseif max(vtreat)<=7
    cyc_type = {'Trend2'};
elseif max(vtreat)<=9
    cyc_type = {'HP'};
elseif max(vtreat)<=10
    cyc_type = {'Percent'};
end


%% 3. Import dates
%--------------------------------------------------------------------------
[~, xlstext] = xlsread('CodeVar.xls',['Date' frequency]);
dates_orig = xlstext;
clear xlstext 

% Find first and last observation
fo = find(strcmp(user_fo,dates_orig));
if isempty(fo)==1; error('Error: please chose a valid starting date'); end
lo = find(strcmp(user_lo,dates_orig));
if isempty(fo)==1; error('Error: please chose a valid starting date'); end

% Update initial observation
dates = dates_orig(fo:lo);
nobs = length(dates);


%% 4. Import regions info
%--------------------------------------------------------------------------
[country2region, label] = xlsread('CodeVar.xls','country2region');
country2region = Num2NaN(country2region);
nregion = max(country2region);
cnames_all = label(2:end,1);
cnames_all_long = label(2:end,2);
clear label

% Import region names
[~, rnames] = xlsread('CodeVar.xls','Region');
rnames = rnames(2:end,1);


%% 5. Loop for each region 
%--------------------------------------------------------------------------
for kk=1:nregion

    disp(['Running loop ' num2str(kk) ' ...'])
    
    % 1. Adjust the number of country according to the group
    ncountry = sum(country2region==kk); 
    cnames = cnames_all(country2region==kk,1);
    cnames_long = cnames_all_long(country2region==kk,1);
    
    % 2. Import data, get cyclical component, shorten the sample (as chosen above)
    for ii=1:nvar
        
        % Import the data, select region, get cyclical component
        [TEMP, ~] = xlsread(['Data' frequency '_57C.xls'],char(vnames(ii)));
        TEMP = Num2NaN(TEMP);
        TEMP = TEMP(fo:lo,country2region==kk);

        % Get cyclical 
        temp = GetCyclical(TEMP,vtreat(ii),lag,lambda);
       
        % Save cyclical data (with selected sample period)
        DATA.(vnames{ii}) = temp;
        clear temp
        if nregion>1
            save( [outdir 'Data' frequency '_' rnames{kk} num2str(kk) '_' cyc_type{1}]);
        else
            save( [outdir 'Data' frequency '_' rnames{kk} '_' cyc_type{1}] );
        end
        
    end
    clearex vtreat cnames_all cnames_all_long dates rnames vnames vnames vnames_long vtreat country2region ii lag ncountry nobs fo lo nregion nvar outdir lambda cyc_type frequency
end
disp('Done!')
 