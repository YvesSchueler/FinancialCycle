function     [FCyc_ecdf,FCyc_ecdf_ma,dY_ecdf,dY_ecdf_ma,weights] = fn_fcycle_dy_rt(dY,idx,corr_rol,sgn_res,T_fix_entire_sample)
%Purpose: Creat real time composite cycle indicators
%
%Coherent financial cycles for G-7 countries: Why credit can be an asset
% by Yves Schüler, Paul Hiebert, and Tuomas Peltonen
% written by Yves Schüler.
% last updated 2018/03/27
d = sum(idx);
y = dY(:,logical(idx));
weights = NaN(length(y),d);
FCyc_ecdf = NaN(length(y),1);
dY_ecdf = NaN(length(y),d);
dY_ecdf_ma = NaN(length(y),d);
idx_nan = min(~isnan(y)')'; %starts when all variables are available
T_nan_start = sum(idx_nan(1:T_fix_entire_sample)); %sample start minus nan (isnide loops)
T_nan_end = T_nan_start + sum(idx_nan(T_fix_entire_sample+1:end)); %sample end minus nan (inside loops)

for ii=1:d %loop ecdf indicators (expanding window)

    temp = y(logical(idx_nan),ii); % drop NaNs
    temp2 = [];   
    for jj = T_nan_start:T_nan_end
        
        [xx,ff] = ecdf(temp(1:jj));
        [~,~,li] = unique(temp(1:jj));
        
        if jj == T_nan_start
        temp2(1:jj,1) = xx(li+1,1);
        else
        temp2(jj,1) = xx(li(end),1);
        end
    end
    dY_ecdf(idx_nan,ii) = temp2;
end

%% time varying correlation aggregation
temp = dY_ecdf(idx_nan,:); 
C = ones(d,d,size(temp,1));
C_vech = ones(size(temp,1),size(temp,2)*(size(temp,2)+1)/2); 
FCyc_sav = NaN(size(temp,1),1);
weights_sav = NaN(size(temp,1),d);

if corr_rol ==1
lam = .89; %smoothing parameter 98
sigma_ij = zeros(size(temp,1)+1,d*(d+1)/2);
sigma_ij(1,:) = vech(nancov(temp(1:8,:)))';


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
    weights(idx_nan,:) = weights_sav; 
 else
    
 FCyc_sav = mean(temp,2); %w/o changing weights

end



%% Fcyc: pass it through ecdf
temp = FCyc_sav;
temp2=[];
for jj = T_nan_start:T_nan_end
        
        [xx,ff] = ecdf(temp(1:jj));
        [~,~,li] = unique(temp(1:jj));
        
        if jj == T_nan_start
        temp2(1:jj,1) = xx(li+1,1);
        else
        temp2(jj,1) = xx(li(end),1);
        end
end

%FCyc_ecdf(idx_nan,1) = FCyc_sav; %w/o ecdf second step
FCyc_ecdf(idx_nan,1) = temp2; %with ecdf second step



%% smoothing part
winz= 6;
weights_p = bartlett(winz*2+1);
weights_p_ = weights_p(winz+1:end-1);
temp = [] ;
w_p = length(weights_p_);
for ii = 1:w_p
temp = [temp FCyc_ecdf(w_p-ii+1:end-ii+1)];
end
FCyc_ecdf_ma = [nan(w_p-1,1) ; temp*weights_p_./sum(weights_p_)];


for ii = 1:d
    dY_ecdf_ma(:,ii) = [nan(5,1); [dY_ecdf(6:end,ii) dY_ecdf(5:end-1,ii) dY_ecdf(4:end-2,ii) dY_ecdf(3:end-3,ii) dY_ecdf(2:end-4,ii) dY_ecdf(1:end-5,ii)]*weights_p_./sum(weights_p_)];
end

%save output w/o training sample
dY_ecdf(1:T_fix_entire_sample,:)= nan(T_fix_entire_sample,d); %save w/o training sample
dY_ecdf_ma(1:T_fix_entire_sample,:)= nan(T_fix_entire_sample,d);
FCyc_ecdf(1:T_fix_entire_sample,:)= nan(T_fix_entire_sample,1);
FCyc_ecdf_ma(1:T_fix_entire_sample,:)= nan(T_fix_entire_sample,1);

end
