%Example code for
%Financial Cycles: Characterisation and real-time measurement
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2021/11/15

clear all
close all
clc

%relevant output:
% T - time dimension
% date_us - cell with dates
% n_co - number of countries
%
% table_years (n_co x 9)    - for FCycle (broad, column 1-3), FCycle (narrow, column 4-6), BCycle (column 7-9): 
%                             [frequency window left, peak, frequency
%                             window right], similar to Table 3 in ESRB paper
%
% FCycle (T x n_co) - band pass filtered financial cycles (broad)
%
%This is the main real-time indicator:
% FCycle_rt (T x n_co) - real time ma filtered financial cycles (broad)
%
% FCycle_n (T x n_co) - band pass filtered financial cycles (narrow)  
% FCycle_n_rt (T x n_co) - real time ma filtered financial cycles (narrow)  
% BCycle (T x n_co) - band pass filtered business cycles   
% BCycle_rt (T x n_co) - real time ma filtered business cycles  
%
% FCycle_unsmooth (T x n_co) - unfiltered financial cycles (broad) (data pre-treated)  
% FCycle_unsmooth_rt (T x n_co) - real time unfiltered financial cycles (broad) (data not pre-treated)  
% FCycle_unsmooth_n (T x n_co) - unfiltered financial cycles (narrow)  (data pre-treated)
% FCycle_unsmooth_n_rt (T x n_co) - real time unfiltered financial cycles (narrow) (data not pre-treated) 
% BCycle_unsmooth (T x n_co) - unfiltered business cycles  (data pre-treated)
% BCycle_unsmooth_rt (T x n_co) - real time unfiltered business cycles  (data not pre-treated)
% 


addpath('functions');

% select countries
%1 - CA; 2 - DE, 3 - FR; 4 - IT; 5 - JP; 6 - UK; 7 - US
co = [7];
n_co = length(co);

% load data
%[dY_fin,fdY_fin,fY_fin,Y_fin,dY_bus,fdY_bus,fY_bus,Y_bus,co_string,T,d] = fn_data_g7(co,1);

load data\us.mat
%definition of variables is
% dY_fin - quarterly growth rates of credit, house, equity, and bond
% dY_bus - quarterly growth rates of gdp, consumption, investment, and hours
% fdY_fin - filtered dY_fin
% fdY_bus - filtered dY_bus
% for fdY: filter dY with fn_bpass_all(dY_fin_mu,2,200,0,0,0), where dY_fin_mu is
%           zero mean
% real financial variables derived by:   Y_fin = [[credit house equity]./[CPI CPI CPI] (1./(1+BONDYIELD/100))./[CPI]];



%% preliminaries
%%% Inputs
q_low = 5;                      % lower bound for spectral density, put 0 for lowest (10=benchmark)
q_high= 200;                    % upper bound for spectral density, put 0 for highest (0 = benchmark)
factor = 8;                     % precision of spectral density estimates
N = 1 ;                         % Number of peaks to select
pct = 0.67;                     % percentage rule for the area of power cohesion to select
corr_rol = 1 ;                  % Composite Cycle: if 1 - rolling corrleations; 0 - linear average
sgn_res = 1;                    % sgn restrictions CISS: emphasises postiviely related movements
idx_fc = [1 1 1 1];
idx_bc = [1 1 1 1];

%grid
[grid,step] = fn_grid(q_low,q_high,16);

table_freqband_fcycle = [];
table_freqband_fcycle_n = [];
table_freqband_bcycle = [];
grid_index_total = [];
grid_index_peak = [];

%FC broad
FCycle = zeros(T,length(co));
FCycle_rt = zeros(T,length(co));
FCycle_unsmooth = zeros(T,length(co));
FCycle_unsmooth_rt = zeros(T,length(co));
series_graph_fc = zeros(T,length(co)*5);
time_varying_weights =  zeros(T,length(co)*4);
%FC narrow
FCycle_n = zeros(T,length(co));
FCycle_n_rt = zeros(T,length(co));
series_graph_fc_n = zeros(T,length(co)*3);
FCycle_unsmooth_n = zeros(T,length(co));
FCycle_unsmooth_n_rt = zeros(T,length(co));
%BC
BCycle = zeros(T,length(co));
BCycle_rt = zeros(T,length(co));
series_graph_bc = zeros(T,length(co)*5);
BCycle_unsmooth = zeros(T,length(co));
BCycle_unsmooth_rt = zeros(T,length(co));


%%%%%%%%%%%%%%%% loop across countries: here, only US
for jip = 1:length(co)

dY_f = fdY_fin; %change for different countries in here... 
dY_b = fdY_bus;
dY_f_rt = dY_fin; 
dY_b_rt = dY_bus;
  
[pw_coh_fc,f_hats_fc] = pw_cohesion(dY_f,grid,factor,idx_fc);
[pw_coh_fc_n,f_hats_fc_n] = pw_cohesion(dY_f,grid,factor,[1 1 0 0]);
[pw_coh_bc,f_hats_bc] = pw_cohesion(dY_b,grid,factor,idx_bc);

peak_fc = zeros(1,N);
peak_fc_n = zeros(1,N);
peak_bc = zeros(1,N);
    
%%%%%%%%%%%%%%%%%%%% N HIGHEST EXTREMA in TOTAL
[temp1,~,~] = extrema(pw_coh_fc);
[temp15,~,~] = extrema(pw_coh_fc_n);
[temp2,~,~] = extrema(pw_coh_bc);
    

%find highest N
[~,sortIndex] = sort(pw_coh_fc(temp1(:,1)),'descend');                                                   
maxIndex = sortIndex(1:N); 
peak_fc(1,:) = (pi/2)./grid(temp1(maxIndex,1));

[~,sortIndex] = sort(pw_coh_fc_n(temp15(:,1)),'descend');                                                   
maxIndex = sortIndex(1:N); 
peak_fc_n(1,:) = (pi/2)./grid(temp15(maxIndex,1));


[~,sortIndex] = sort(pw_coh_bc(temp2(:,1)),'descend');                                                   
maxIndex = sortIndex(1:N); 
peak_bc(1,:) = (pi/2)./grid(temp2(maxIndex,1));


    
%derive area and distance interval
[min_distance_fc]=fn_integral_min_dist(pw_coh_fc,peak_fc(1,1),grid,pct);
[min_distance_fc_n]=fn_integral_min_dist(pw_coh_fc_n,peak_fc_n(1,1),grid,pct);
[min_distance_bc]=fn_integral_min_dist(pw_coh_bc,peak_bc(1,1),grid,pct);

%store freq windows and peaks
table_freqband_fcycle(jip,:) = [min_distance_fc(1,3) peak_fc(1,:) min_distance_fc(1,4)];
table_freqband_fcycle_n(jip,:) = [min_distance_fc_n(1,3) peak_fc_n(1,:) min_distance_fc_n(1,4)];
table_freqband_bcycle(jip,:) = [min_distance_bc(1,3) peak_bc(1,:) min_distance_bc(1,4)];

%%%
%Construct financial and business cycles
%%%

%%%SMOOTHED
%Financial Cycle: Broad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[series_graph_fc(:,5*jip-4:5*jip),FCycle(:,jip),~,~,FCycle_unsmooth(:,jip),time_varying_weights(:,(jip-1)*4+1:(jip-1)*4+4)] = fn_cycle_extract_ecdf_aggr(dY_f,idx_fc,min_distance_fc,corr_rol,sgn_res);
%Financial Cycle: Narrow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[series_graph_fc_n(:,3*jip-2:3*jip),FCycle_n(:,jip),~,~,FCycle_unsmooth_n(:,jip),~] = fn_cycle_extract_ecdf_aggr(dY_f,[1 1 0 0],min_distance_fc_n,0,sgn_res);
%Business Cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[series_graph_bc(:,5*jip-4:5*jip),BCycle(:,jip),~,~,BCycle_unsmooth(:,jip)] = fn_cycle_extract_ecdf_aggr(dY_b,idx_bc,min_distance_bc,corr_rol,sgn_res);


%%%REAL TIME (3 years training sample: can be adjusted by changing 12 to ??
%%% in fn_fcycle_dy_rt(~,~,~,~,12)
%Financial Cycle: Broad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[FCycle_unsmooth_rt(:,jip),FCycle_rt(:,jip),dY_ecdf,dY_ecdf_ma,weights] = fn_fcycle_dy_rt(dY_f_rt,[1 1 1 1],1,1,12); %3years training sample
%Financial Cycle: Narrow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[FCycle_unsmooth_n_rt(:,jip),FCycle_n_rt(:,jip),~,~,~] = fn_fcycle_dy_rt(dY_f_rt,[1 1 0 0],1,1,12); %3years training sample
%Business Cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BCycle_unsmooth_rt(:,jip),BCycle_rt(:,jip),~,~,~] = fn_fcycle_dy_rt(dY_b_rt,[1 1 1 1],1,1,12); %3years training sample

end

table_years = [table_freqband_fcycle table_freqband_fcycle_n table_freqband_bcycle];



%mark business cycle frequencies
freq_low = 2*pi/(2*4);
freq_high= 2*pi/(8*4);
bus_freq = (grid <= freq_low) & (grid >= freq_high);
freq_low_fc =  2*pi/(8*4);
freq_high_fc = 2*pi/(20*4);
fc_freq = (grid <= freq_low_fc) & (grid >= freq_high_fc);
grid_select = grid <= grid(end); %choose window size
      

 % plot each variable QUARTERLY
    Ynam = {'\Deltacr'; '\Deltap_h'; '\Deltap_e'; '\Deltap_b'; '\Deltaq'; '\Deltaco'; '\Deltai'; '\Deltah'}; 

 %FINANCIAL CROSS SPECTRA   
   size_plot = [8/2 8/2];
    ii=7;
     f_hat_fc_i = f_hats_fc;
     figure
     hold on
     area(grid(grid_select),bus_freq(grid_select)*max(f_hat_fc_i(:)))
     colormap([0.69 0.7686 0.87])
     area(grid(grid_select),fc_freq(grid_select)*max([f_hat_fc_i(:)]),'FaceColor',[230/255 230/255 250/255])
     a1=   plot(grid(grid_select),f_hat_fc_i(grid_select,1) ,'-', 'Color',[0 0 1],'LineWidth',1.5);
     a2=   plot(grid(grid_select),f_hat_fc_i(grid_select,2) ,'-','Color',[0 1 0],'LineWidth',1.5);
     a3=   plot(grid(grid_select),f_hat_fc_i(grid_select,3) ,'-','Color',[1 0 0],'LineWidth',1.5);
     a4=   plot(grid(grid_select),f_hat_fc_i(grid_select,4) ,'-','Color',[1 0 1],'LineWidth',1.5);
     a5=   plot(grid(grid_select),f_hat_fc_i(grid_select,5) ,'-','Color',[0 1 1],'LineWidth',1.5);
     a6=   plot(grid(grid_select),f_hat_fc_i(grid_select,6) ,'-','Color',[0 0 0],'LineWidth',1.5);
 
     hold off
     %L =length(grid);
     lab_x = [grid(18)  grid(56) grid(248)];
     set(gca, 'XTick', lab_x , 'XTickLabel', {'20'; '8'; '2'},'FontSize',10)
     l=   legend([a1,a2,a3,a4,a5,a6],[Ynam{1} '/' Ynam{2}],[Ynam{1} '/' Ynam{3}],[Ynam{1} '/' Ynam{4}],[Ynam{2} '/' Ynam{3}],[Ynam{2} '/' Ynam{4}],[Ynam{3} '/' Ynam{4}]); %'Location','northoutside','orientation','horizontal'
     set(l,'FontSize',10)
    legend('boxoff')
     ylabel('$|\hat{s}_{x_ix_j}|$','Interpreter','Latex','FontSize',15)
     xlabel('Years')
     axis tight
     set(gcf,'PaperPositionMode','manual')
     set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[0 0 size_plot]);
     figfile = fullfile('figures', ['f_hat_fc_' num2str(ii)]);
     saveas(gcf, figfile , 'epsc'); 
 

 
 %BUSINESS CYCLE CROSS SPECTRA
   size_plot = [8/2 8/2];
     f_hat_fc_i = f_hats_bc;
     figure
     hold on
     area(grid(grid_select),bus_freq(grid_select)*max(f_hat_fc_i(:)))
     colormap([0.69 0.7686 0.87])
     area(grid(grid_select),fc_freq(grid_select)*max([f_hat_fc_i(:)]),'FaceColor',[230/255 230/255 250/255])
     a1=   plot(grid(grid_select),f_hat_fc_i(grid_select,1) ,'-','Color',[0 0 1],'LineWidth',1.5);
     a2=   plot(grid(grid_select),f_hat_fc_i(grid_select,2) ,'-','Color',[0 1 0],'LineWidth',1.5);
     a3=   plot(grid(grid_select),f_hat_fc_i(grid_select,3) ,'-','Color',[1 0 0],'LineWidth',1.5);
     a4=   plot(grid(grid_select),f_hat_fc_i(grid_select,4) ,'-','Color',[1 0 1],'LineWidth',1.5);
     a5=   plot(grid(grid_select),f_hat_fc_i(grid_select,5) ,'-','Color',[0 1 1],'LineWidth',1.5);
     a6=   plot(grid(grid_select),f_hat_fc_i(grid_select,6) ,'-','Color',[0 0 0],'LineWidth',1.5);
 
     hold off
     %L =length(grid);
     lab_x = [grid(18)  grid(56) grid(248)];
     set(gca, 'XTick', lab_x , 'XTickLabel', {'20'; '8'; '2'},'FontSize',10)
     l=   legend([a1,a2,a3,a4,a5,a6],[Ynam{5} '/' Ynam{6}],[Ynam{5} '/' Ynam{7}],[Ynam{5} '/' Ynam{8}],[Ynam{6} '/' Ynam{7}],[Ynam{6} '/' Ynam{8}],[Ynam{7} '/' Ynam{8}]); %'Location','northoutside','orientation','horizontal'
     set(l,'FontSize',10)
    legend('boxoff')
     ylabel('$|\hat{s}_{x_ix_j}|$','Interpreter','Latex','FontSize',15)
     xlabel('Years')
     axis tight
     set(gcf,'PaperPositionMode','manual')
     set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[0 0 size_plot]);
     figfile = fullfile('figures', ['f_hat_bc_' num2str(ii)]);
     saveas(gcf, figfile , 'epsc'); 
 
 
 
 % print financial and real cycle power cohesion in one graph.
   size_plot = [8/2 8/2];
   jip=7;
        
%mark business cycle frequencies
freq_low = 2*pi/(2*4);
freq_high= 2*pi/(8*4);
bus_freq = (grid <= freq_low) & (grid >= freq_high);
freq_low_fc =  2*pi/(8*4);
freq_high_fc = 2*pi/(20*4);
fc_freq = (grid<= freq_low_fc) & (grid >= freq_high_fc);
grid_s = (grid <= grid(end)) & (grid>= grid(1)); %choose window size

      
        figure('Name','Power Cohesion')
        hold on
        area(grid(grid_s),bus_freq(grid_s)*max([pw_coh_fc(grid_s); pw_coh_bc(grid_s);pw_coh_fc_n(grid_s)]))
        colormap([0.69 0.7686 0.87])
        area(grid(grid_s),fc_freq(grid_s)*max([[pw_coh_fc(grid_s); pw_coh_bc(grid_s);pw_coh_fc_n(grid_s)]]),'FaceColor',[230/255 230/255 250/255])
        a1=plot(grid(grid_s),pw_coh_fc(grid_s,1),'-','Color', [0 0 0],'LineWidth',3); %yearly fc
        a2=plot(grid(grid_s),pw_coh_fc_n(grid_s,1),'-','Color', [0 76/255 153/255],'LineWidth',3); %yearly fc
        a3=plot(grid(grid_s),pw_coh_bc(grid_s,1),'-','Color',[0.7 0 0],'LineWidth',3); %yearly bc
        hold off
        lab_x = [grid(18)  grid(56) grid(248)];
        set(gca, 'XTick', lab_x , 'XTickLabel', {'20'; '8'; '2'},'FontSize',10)
        title([co_string{jip}],'FontSize', 11);
        l=   legend([a1,a2,a3],'Broad financial cycle','Narrow financial cycle','Business cycle');%'Location','northoutside','orientation','horizontal'
         set(l,'FontSize',10)
        legend('boxoff')
        ylabel('PCoh','FontSize',11)
        xlabel('Years')
        axis tight
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperPosition',[0 0 size_plot]);   
        figfile = fullfile('figures\', co_string{jip});
        saveas(gca, figfile , 'epsc'); 
      

   
