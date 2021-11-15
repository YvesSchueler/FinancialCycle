function     [series_graph,FCycle,dY_ecdf,C_out,FCyc_sav_unsmooth,weights] = fn_cycle_extract_ecdf_aggr(dY,idx,min_distance,corr_rol,sgn_res)
%Purpose: Create composite cycle indicators (using full sample information)
%
%Coherent financial cycles for G-7 countries: Why credit can be an asset
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2018/03/27
d = sum(idx);
y = dY(:,logical(idx));
temp = NaN(length(y),d);
fdY = NaN(length(y),d);
temp_indicator_smooth = NaN(length(y),d);
FCyc_sav_unsmooth = NaN(length(y),1);
FCycle = NaN(length(y),1);
weights = NaN(length(y),d);
C_out = NaN(length(y),d*(d+1)/2);
for ii=1:d
        idx_nan = isnan(y(:,ii));
        idx_nan = 1-idx_nan; % variable specific index
        temp2 = y(logical(idx_nan),ii); % drop NaNs
        [temp2_fcast,length_fcast] = fn_forecast_ar(temp2); 
        temp10 = fn_bpass_all(temp2_fcast-mean(temp2_fcast),min_distance(1,4)*4,min_distance(1,3)*4,0,0,0);
        temp10 = temp10(1:end-length_fcast,:)+mean(temp2);
        
        [xx,ff] = ecdf(temp2);
        [~,~,li] = unique(temp2);
        temp(logical(idx_nan),ii)  = xx(li+1,1);
        fdY(logical(idx_nan),ii) = temp10;
        
        x  = [1 0 ];
        x0 = [1 0];
        [x,fval1] = fminunc(@(x)myfun(x,temp10,xx(li+1,1)),x0,optimset('MaxFunEvals',50000, 'LargeScale', 'off' , 'MaxIter', 50000, 'TolX', 1.0e-5, 'TolFun', 1.0e-5,'Display','off'));
        temp_adjust = temp10*x(1)+x(2);

        temp_adjust(temp_adjust > 1) = 1; % cut off max and min
        temp_adjust(temp_adjust < 0) = 0;

        temp_indicator_smooth(logical(idx_nan),ii) = temp_adjust;
        
end
dY_ecdf = temp;
idx_cyc_sav = logical(min(~isnan(temp)')');
temp = temp(idx_cyc_sav,:); 
C = ones(d,d,size(temp,1));
C_vech = ones(size(temp,1),size(temp,2)*(size(temp,2)+1)/2); 
FCyc_sav = NaN(size(temp,1),1);
weights_sav = NaN(size(temp,1),d);

if corr_rol ==1
lam = .89; %smoothing parameter
sigma_ij = zeros(size(temp,1)+1,d*(d+1)/2);
sigma_ij(1,:) = vech(nancov(temp(1:8,:)))'; %initialize cov


 for ii = 1:length(temp)
   
   ttt = 0;
   sgn_sigma = zeros(1,d*(d+1)/2);
   for jj = 1:d
   for kk = jj:d
   ttt = ttt + 1;
   sigma_ij(ii+1,ttt) = lam*sigma_ij(ii,ttt) + (1-lam)*(temp(ii,jj)-0.5)*(temp(ii,kk)-0.5);
   if sigma_ij(ii+1,ttt) > 0
       sgn_sigma(1,ttt) = sigma_ij(ii+1,ttt);
   end
   end
   end
   
   if sgn_res == 1        
   temp_vech = ivech(sgn_sigma'); 
   else
   temp_vech = ivech(sigma_ij(ii+1,:)');
   end
   
   temp_cov = temp_vech + temp_vech' - diag(diag(temp_vech)); %COVARIANCE   
   C(:,:,ii) = sqrt(diag(diag(temp_vech)))\temp_cov/sqrt(diag(diag(temp_vech))); %CORRELATION
   C_vech(ii,:) = vech(squeeze(C(:,:,ii)))';
   FCyc_sav(ii,1) = sum(squeeze(C(:,:,ii)))*temp(ii,:)'/sum(sum(squeeze(C(:,:,ii)),1),2); %CORRELATION WEIGHTED INDICES
   weights_sav(ii,:) = sum(squeeze(C(:,:,ii)))/sum(sum(squeeze(C(:,:,ii)),1),2);
 end
    weights(idx_cyc_sav,:) = weights_sav; 
else
    
 FCyc_sav = mean(temp,2);

end
temp = FCyc_sav;
[xx,ff] = ecdf(temp);
[~,~,li] = unique(temp);
FCyc_sav_unsmooth(idx_cyc_sav,1) = xx(li+1,1);


[temp_fcast,length_fcast] = fn_forecast_ar(xx(li+1,1)); 
FCyc_sav_smooth = fn_bpass_all(temp_fcast-mean(temp_fcast),min_distance(1,4)*4,min_distance(1,3)*4,0,0,0);
xx_mu = mean(xx(li+1,1));
FCyc_sav_smooth = FCyc_sav_smooth(1:end-length_fcast,1)+xx_mu;


x  = [1 0 ];
x0 = [1 0];
[x,fval1] = fminunc(@(x)myfun(x,FCyc_sav_smooth,xx(li+1,1)),x0,optimset('MaxFunEvals',50000, 'LargeScale', 'off' , 'MaxIter', 50000, 'TolX', 1.0e-5, 'TolFun', 1.0e-5,'Display','off'));
temp_adjust = FCyc_sav_smooth*x(1)+x(2);

temp_adjust(temp_adjust > 1) = 1; % cut off max and min
temp_adjust(temp_adjust < 0) = 0;

FCycle(idx_cyc_sav,1) = temp_adjust;%zscore(FCyc_sav);
C_out(idx_cyc_sav,:) = C_vech;
series_graph = [FCycle  temp_indicator_smooth];

end
