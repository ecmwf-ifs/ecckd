REF='/hugetmp/parr/ckdmip/evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_present.h5';
DIR='/hugetmp/parr/ecckd';
TOL=[0.16 0.08 0.04 0.02 0.01 0.005];
APP='global-nwp';
%APP='limited-area-nwp';
SCENARIO='present';
BANDSTR={'narrow','wide','fsck'};

nb = length(BANDSTR);
nt = length(TOL);
RRTMG='/hugetmp/parr/ckdmip/results/ecrad-rrtmg/lw_fluxes/ecrad-rrtmg_evaluation1_lw_climate_narrow-140_fluxes_present.nc';
figure(1)
stats = evaluate_ckd_lw_fluxes(REF,RRTMG,'RRTMG','Present');
rrtmg(1) = stats.toa_up_bias;
rrtmg(2) = stats.toa_up_rmse;
rrtmg(3) = stats.surf_dn_bias;
rrtmg(4) = stats.surf_dn_rmse;
rrtmg(5) = stats.heating_rate_high_rmse;
rrtmg(6) = stats.heating_rate_low_rmse;

SUFF='-t2';
RSUFF='-t2';
SUFF = '';
RSUFF = '';
if ~exist('ng','var')
  clear stat rstat
  for ib = 1:nb
    for it = 1:nt
      CKD=[DIR '/lw_ckd/lw_ckd_' APP '_' BANDSTR{ib} '_tol' num2str(TOL(it)) SUFF '.nc'];
      FLUX=[DIR '/lw_optical-depth/lw_raw-ckd_' APP '_' BANDSTR{ib} '-tol' num2str(TOL(it)) RSUFF '_optical-depth_' SCENARIO '.nc'];
      figure(1);
      stats = evaluate_ckd_lw_fluxes(REF,FLUX,'ecCKD','Present');
      drawnow

      ckd = loadnc(CKD)
      ng(it,ib) = size(ckd.gpoint_fraction,2);
      rstat(it,ib,1) = stats.toa_up_bias;
      rstat(it,ib,2) = stats.toa_up_rmse;
      rstat(it,ib,3) = stats.surf_dn_bias;
      rstat(it,ib,4) = stats.surf_dn_rmse;
      rstat(it,ib,5) = stats.heating_rate_high_rmse;
      rstat(it,ib,6) = stats.heating_rate_low_rmse;

      FLUX=[DIR '/lw_optical-depth/lw_ckd_' APP '_' BANDSTR{ib} '-tol' num2str(TOL(it)) SUFF '_optical-depth_' SCENARIO '.nc'];
      stats = evaluate_ckd_lw_fluxes(REF,FLUX,'ecCKD','Present');
      drawnow;
      stat(it,ib,1) = stats.toa_up_bias;
      stat(it,ib,2) = stats.toa_up_rmse;
      stat(it,ib,3) = stats.surf_dn_bias;
      stat(it,ib,4) = stats.surf_dn_rmse;
      stat(it,ib,5) = stats.heating_rate_high_rmse;
      stat(it,ib,6) = stats.heating_rate_low_rmse;
    end
  end
end

stat_title = {'Bias in TOA upwelling (W m^{-2})',...
	      'RMSE in TOA upwelling (W m^{-2})',...
	      'Bias in surface downwelling (W m^{-2})',...
	      'RMSE in surface downwelling (W m^{-2})',...
	      'RMSE in 4-1100 hPa heating rate (K d^{-1})',...
	      'RMSE in 0.02-4 hPa heating rate (K d^{-1})'};

nx = 2;
ny = 3;

cols = {'k','r','g'};

figure(5)
clf
for istat = 1:size(stat,3)
  subplot(ny,nx,istat)
  for ib = 1:nb
    plot(ng(:,ib), rstat(:,ib,istat), [cols{ib} '--'])
    hold on
    plot(ng(:,ib), stat(:,ib,istat), cols{ib})
    hold on
    plot(140, rrtmg(istat), 'mo');
  end
  ylabel(stat_title{istat});
  if istat == 1
    legend(BANDSTR);
  end
  xlim([0 150]);
end



