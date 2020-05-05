if ~exist('d','var')
  w = loadnc('/hugetmp/parr/ecckd/lw_order/lw_order_fsck_composite.h5');
  wn = w.wavenumber;
  d=loadnc('/hugetmp/parr/ecckd/lw_gpoints/lw_gpoints_global-nwp_fsck_tol0.01.h5');  nx = 1; ny = 1;
  d=loadnc('/hugetmp/parr/ecckd/lw_gpoints/lw_gpoints_global-nwp_wide_tol0.005.h5');  nx = 3; ny = 2;
  d=loadnc('/hugetmp/parr/ecckd/lw_gpoints/lw_gpoints_global-nwp_narrow_tol0.01.h5');  nx = 5; ny = 3;
end

gases = {'composite','h2o','o3'};
gasname = {'Composite','H_2O','O_3'};
ngas = length(gases);

nband = length(d.wavenumber1_band);


ig = 1;
index_band = [];

for iband=1:nband
  ng_local = 1-ngas;
  for igas = 1:ngas
    ng_local = ng_local + d.([gases{igas} '_n_g_points'])(iband);
  end
  ng(iband) = ng_local;
  index_band(ig:(ig-1+ng(iband))) = iband;
  iig = find(d.g_point >= ig-1 & d.g_point < ig+ng(iband)-1);
  wn1(iband) = round(w.wavenumber(min(iig)));
  wn2(iband) = round(w.wavenumber(max(iig)));
  ig = ig + ng(iband);
end


clf

for iband=1:nband
  subplot(ny,nx,iband);
  cols = get(gca,'colororder');
  cols = cols.*0.5 + 0.5;
  ncol = size(cols,1);

  for igas = 1:ngas
    index = find(index_band == iband);
    gmin{igas} = d.([gases{igas} '_g_min'])(index);
    gmax{igas} = d.([gases{igas} '_g_max'])(index);

    gmin{igas} = gmin{igas}-gmin{igas}(1);
    gmax{igas} = gmax{igas}-gmax{igas}(1);
    nggas(igas) = d.([gases{igas} '_n_g_points'])(iband);
  end
  ngg = nggas;
  nggas(:) = 1;
  icol = 1;
  for ig = 1:ng(iband)
    if gmin{1}(ig) == 0
      h=patch([0 0 0 0],...
	      [gmin{2}(ig) gmin{2}(ig) gmax{2}(ig)+1 gmax{2}(ig)+1]./nggas(2),...
	      [gmin{3}(ig) gmax{3}(ig)+1 gmax{3}(ig)+1 gmin{3}(ig)]./nggas(3),'r');
      set(h,'facecolor',cols(icol,:),'edgecolor','k');
      hold on
    end
    if gmin{2}(ig) == 0
      h=patch([gmin{1}(ig) gmin{1}(ig) gmax{1}(ig)+1 gmax{1}(ig)+1]./nggas(1),...
	      [0 0 0 0],...
	      [gmin{3}(ig) gmax{3}(ig)+1 gmax{3}(ig)+1 gmin{3}(ig)]./nggas(3),'r');
      set(h,'facecolor',cols(icol,:),'edgecolor','k');
      hold on
    end
    if gmin{3}(ig) == 0
      h=patch([gmin{1}(ig) gmin{1}(ig) gmax{1}(ig)+1 gmax{1}(ig)+1]./nggas(1),...
	      [gmin{2}(ig) gmax{2}(ig)+1 gmax{2}(ig)+1 gmin{2}(ig)]./nggas(2),...
	      [0 0 0 0],'r');
      set(h,'facecolor',cols(icol,:),'edgecolor','k');
      hold on
    end
    icol = icol + 1;
    if icol > ncol
      icol = 1;
    end
  end
  set(gca,'proj','pers');
  view([-1 -1 -1])
  axis([0 ngg(1) 0 ngg(2) 0 ngg(3)]);
  daspect(ngg);
  camlight left
  xlabel(gasname{1});
  ylabel(gasname{2});
  zlabel(gasname{3});
  set(gca,'xticklabel','','yticklabel','','zticklabel','');
  title([num2str(wn1(iband)) '-' num2str(wn2(iband)) ' cm^{-1}: n_g=' num2str(ng(iband))]);
end
