REF='/hugetmp/parr/ckdmip/evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_present.h5';

bs = {'fsck','wide','narrow'};
tol = {'0.01','0.04'};
app = {'global-nwp','limited-area-nwp','climate'};
appshort = {'GNWP','LANWP','CLIM'};
figure(1)
for iapp = 3;%1:length(app)
  for ibs = 1:length(bs)
    for itol = 1:length(tol);
      BASE = ['lw_ckd_' app{iapp} '_' bs{ibs} '-tol' tol{itol}];
      if strcmp(app{iapp},'climate')
	NEWBASE = ['ecckd-ckd_evaluation1_lw_climate_' bs{ibs} '-tol' tol{itol}];
      else
	NEWBASE = BASE;
      end
      CKD=['/hugetmp/parr/ecckd/lw_optical-depth/' NEWBASE '_optical-depth_present.nc'];
      PNG=[BASE '_present.png'];
      tt = ['ecCLD/' appshort{iapp} '/' bs{ibs} '/' tol{itol}];
      clf
      evaluate_ckd_lw_fluxes(REF,CKD,tt,'Present')
      drawnow
      print_png(PNG,'120')
    end
  end
end
