REF='/hugetmp/parr/ckdmip/evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_present.h5';
CKDr='/hugetmp/parr/ecckd/lw_optical-depth/lw_raw-ckd_global-nwp_fsck-tol0.04_optical-depth_present.nc';
CKD='/hugetmp/parr/ecckd/lw_optical-depth/lw_ckd_global-nwp_fsck-tol0.04_optical-depth_present.nc';
CKDrl='/hugetmp/parr/ecckd/lw_optical-depth/lw_raw-ckd_global-nwp_fsck-tol0.04-log_optical-depth_present.nc';
CKDl='/hugetmp/parr/ecckd/lw_optical-depth/lw_ckd_global-nwp_fsck-tol0.04-log_optical-depth_present.nc';
CKDrt='/hugetmp/parr/ecckd/lw_optical-depth/lw_raw-ckd_global-nwp_fsck_tol0.04-ln_optical-depth_present.nc';
CKDt='/hugetmp/parr/ecckd/lw_optical-depth/lw_ckd_global-nwp_fsck_tol0.04-ln_optical-depth_present.nc';

figure(1); evaluate_ckd_lw_fluxes(REF,CKD,'ecCKD','Present')
figure(2); evaluate_ckd_lw_fluxes(REF,CKDr,'Raw','Present')
figure(3); evaluate_ckd_lw_fluxes(REF,CKDl,'Logarithmic','Present')
figure(4); evaluate_ckd_lw_fluxes(REF,CKDrl,'Raw logarithmic','Present')
figure(5); evaluate_ckd_lw_fluxes(REF,CKDt,'T2','Present')
figure(6); evaluate_ckd_lw_fluxes(REF,CKDrt,'Raw T2','Present')
