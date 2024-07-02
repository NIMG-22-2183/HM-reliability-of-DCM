%  This code can be used to model LFP data using DCM for CSD
% Reliability of dynamic causal modelling of resting state magnetoencephalography 
%==========================================================================
clc
clear all
close all
addpath('C:\spm12')
spm('defaults','eeg')
%----------------------------Karl comments---------------------------------

% Folder set up
%==========================================================================
f3               = 'C:\DMN_EO_open' ; % Folder address
base             = {'BL', 'TW'};
subject_name     = {'P1??','P1??'}; % please add name ID! 
D                = {};
NumofSub         = 14; % number of subjects
% First level DCM
%==========================================================================
for i  =1:NumofSub
    j                    =1;
    t =1;
    folder               = f3;
    subject              = subject_name{1,i}   ;
    data_folder          = fullfile(folder,base{1,j},subject);
    meg_data             = spm_select('Fplist',data_folder,'^.*_rest_open_DMN.mat$');
    sub                  = append(base{1,j},'_', subject,'_',num2str(t));
    D1{i,1}              = DCM_CSD_HB(sub,meg_data,t);    
end

% Plot --------------------DCM results-------------------------------------

for i = 1: length(subject_name)
    
    DCM              = P1{i,1}; % Template DCM to get the name of channel
    DCM.xY.Ic        = 1;
    DCM.xY.name      = {'MPFC'};
    
    % DCM.xY.Ic      = 2;
    % DCM.xY.name    = {'PRC'};
   
    % DCM.xY.Ic      = 3;
    % DCM.xY.name    = {'LAG'};
    
    % DCM.xY.Ic      = 4;
    % DCM.xY.name    = {'RAG'};
    
    k         =1;
    Hz        =  DCM.M.Hz;
    cc        = DCM.xY.Ic;
    y_BL      = abs(P1{i,1}.Hc{k}(:,cc,cc)); %predicated BL
    y_TW      = abs(P2{i,1}.Hc{k}(:,cc,cc)); %predicated TW
    p_BL      = abs(P1{i,1}.Hc{k}(:,cc,cc)+ P1{i,1}.Rc{k}(:,cc,cc)); % observed response BL
    p_TW      = abs(P2{i,1}.Hc{k}(:,cc,cc)+ P2{i,1}.Rc{k}(:,cc,cc)); % observed response TW
    
    figure('color','white','units','centimeters','position',[2 2 5 5],'papersize',[15.6 4],'filename',[subject_name{1,i},'-', DCM.xY.name{1,1}])
    plot (Hz,y_TW,'k',Hz,p_TW,'--k', Hz,y_BL,'r', Hz,p_BL,'--r','LineWidth',2)
    xlabel('frequency (Hz)','fontsize',12)
    ylabel('PSD(\mu V^{2}/Hz)','fontsize',12)
    legend({'predicted TW', 'observed TW','predicted BL ','observed BL'}, 'Location','northeast')
    legend('boxoff')
    box off
    % saveas(gcf,['fit_',num2str(j)],'jpg')
    % close all   
end



function DCM        = DCM_CSD_HB(sub,meg_data,t)
%  DCM for CSD. This code can be used to model LFP data
% Reliability of dynamic causal modelling of resting state magnetoencephalography 
%==========================================================================
% Spesify (un-estiamted) DCM
%==========================================================================

DCM.options.analysis = 'CSD';       % DCM for cross-spectral density 
DCM.options.spatial  = 'LFP';       % virtual electrode data as data 
DCM.options.model    = 'CMM_NMDA';  % generative model is  spm_fx_cmc

tim_ax(1)            = 1;           % first sample of MEG data (time domain)
tim_ax(2)            = 501;        % last sample of MEG data (time domain)
DCM.options.Tdcm     = [tim_ax(1) tim_ax(2)];  % time axis in milli seconds

DCM.options.Fdcm     = [1   30];    % frequency range  
DCM.M.Hz             = (1 : 30)';
DCM.options.D        = 1;            % downsampling
DCM.Sname            = {'MPFC', 'PRC', 'LAG','RAG'}; % name(s) of the source

% Adress to data
%==========================================================================
 
DCM.xY.Dfile         = meg_data; 

% if model one LFP channel 
%==========================================================================

DCM.xY.modality      = 'LFP';
DCM.xY.Ic            =  1:4;
DCM.xY.name          = {'MPFC', 'PRC', 'LAG','RAG'}; 
DCM.options.trials   =  t; % even 1 and odd 2

% A, B ,C parameters spesficiation in the generative model  spm_fx_cmc
% for one region source  A =[0]. For one condtion B = [0] 
%==========================================================================
Ns                  = 4; 
MPFC                = 1 ; 
PRC                 = 2; 
LAG                 = 3; 
RAG                 = 4;  
%-------------------------------------------------------------------------
DCM.A{1}             = ones(Ns,Ns);
DCM.A{1}             = DCM.A{1};

DCM.A{2}             = ones(Ns,Ns);
DCM.A{2}             = DCM.A{2};

DCM.A{3}             = zeros(Ns,Ns);
DCM.B                    = {};
DCM.C                    = [ 0 0 0 0]; 
DCM.M.dipfit.Nm          = 1;
DCM.M.dipfit.model       = DCM.options.model;
DCM.M.dipfit.type        = DCM.options.spatial;
DCM.M.dipfit.Nc          = Ns;
DCM.M.dipfit.Ns          = Ns;
DCM.M.U                  = eye(Ns,Ns); 
DCM.name                 = sub ;     

% optional to define
%==========================================================================

[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);
DCM.M.pE   = pE;
DCM.M.pC   = pC;

% spm_dcm_csd_data needs to be called for one channel data within many lfps.
%==========================================================================

DCM.options.DATA     = 0;
DCM                  = spm_dcm_csd_data(DCM);

% Model inversion/estiamtion
%==========================================================================

DCM                  = spm_dcm_csd(DCM);

return
%**************************************************************************
% Results can be  plotted using the following :
%==========================================================================

%  spm_dcm_csd_results(DCM,'spectral data');%na
%  spm_dcm_csd_results(DCM,'Coupling (A)'); % not useful for 1 region
%  spm_dcm_csd_results(DCM,'Coupling (B)'); % not useful for 1 condtion
%  spm_dcm_csd_results(DCM,'Coupling (C)'); %na
%  spm_dcm_csd_results(DCM,'trial-specific effects');%na
%  spm_dcm_csd_results(DCM,'Input'); % Spectrum of 1/f b and c ch noise 
% spm_dcm_csd_results(DCM,'Transfer functions'); % Transfer function of CMC
%  spm_dcm_csd_results(DCM,'Cross-spectra (sources)'); 
%  spm_dcm_csd_results(DCM,'Cross-spectra (channels)'); 
%  spm_dcm_csd_results(DCM,'Coherence (sources)')% na
%  spm_dcm_csd_results(DCM,'Coherence (channels)') %na
%  spm_dcm_csd_results(DCM,'Covariance (sources)');
%  spm_dcm_csd_results(DCM,'Covariance (channels)');
end