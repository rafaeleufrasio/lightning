function NGC337_mstar_sfr_sfh_plots, results
;+
; Name
; ----
;   NGC337_mstar_sfr_sfh_plots
;
; Purpose
; -------
;   Make a figure of the stellar mass and SFR posteriors,
;   along with the SFH.
;
; Input
; -----
;   results : structure
;       The IDL structure containing the fit to NGC 337.
;
; Output
; ------
;   Returns an IDL graphic object containing the plots
;
;-


    ; These are suffixed '_chain' to denote that they are
    ; vectors of 1000 samples
    mstar_chain = results.MSTAR ; Total surviving stellar mass
    mstar_percentiles = percentile(mstar_chain, [0.16, 0.50, 0.84])
    SFR_chain = 0.1 * results.PSI[0,*] + 0.9 * results.PSI[1,*] ; SFR averaged over the last 100 Myr
    SFR_percentiles = percentile(SFR_chain, [0.16, 0.50, 0.84])
    SFH_percentiles = results.PSI_PERCENTILES
    age_bin_edges = results.STEPS_BOUNDS

    ; Mass histogram -- bin size is determined by Scott's normal reference rule
    logmstar_chain = alog10(mstar_chain)
    mstar_binsize = 3.49d * stddev(logmstar_chain) / (double(n_elements(logmstar_chain)))^(1/3.d)
    mstar_nbins = ceil((max(logmstar_chain) - min(logmstar_chain)) / mstar_binsize)

    mstar_pdf = histogram(logmstar_chain,$
                          min=min(logmstar_chain),$
                          max=max(logmstar_chain),$
                          nbins=mstar_nbins,$
                          locations=mstar_binloc)

    mstar_pdf = mstar_pdf / mstar_binsize / n_elements(logmstar_chain)

    p1 = plot(mstar_binloc,$
             mstar_pdf,$
             linestyle='-',$
             thick=2,$
             /stair)

    ; Plot the percentiles on top
    yrange = p1.yrange
    p1 = plot(alog10([mstar_percentiles[0], mstar_percentiles[0]]),$
              yrange,$
              linestyle='--',$
              color='orange',$
              thick=2,$
              /overplot)

    p1 = plot(alog10([mstar_percentiles[1], mstar_percentiles[1]]),$
              yrange,$
              linestyle='-',$
              color='orange',$
              thick=2,$
              /overplot)

    p1 = plot(alog10([mstar_percentiles[2], mstar_percentiles[2]]),$
              yrange,$
              linestyle='--',$
              color='orange',$
              thick=2,$
              /overplot)

    hist_width = 0.35 ; in percentage of the figure width

    p1.position = [0.1, 0.9 - hist_width, 0.1 + hist_width, 0.9]
    p1.xtitle = 'log M$_{\star}\ [M_{\odot}]$'
    p1.ytitle = 'PDF'

    ; SFR histogram
    SFR_binsize = 3.49d * stddev(SFR_chain) / (double(n_elements(SFR_chain)))^(1/3.d)
    SFR_nbins = ceil((max(SFR_chain) - min(SFR_chain)) / SFR_binsize)

    SFR_pdf = histogram(SFR_chain,$
                        min=min(SFR_chain),$
                        max=max(SFR_chain),$
                        nbins=SFR_nbins,$
                        locations=SFR_binloc)

    SFR_pdf = SFR_pdf / SFR_binsize / n_elements(SFR_chain)

    p2 = plot(SFR_binloc,$
             SFR_pdf,$
             linestyle='-',$
             thick=2,$
             /stair,$
             /current)

    ; Plot the percentiles on top
    yrange = p2.yrange
    p2 = plot([SFR_percentiles[0], SFR_percentiles[0]],$
              yrange,$
              linestyle='--',$
              color='orange',$
              thick=2,$
              /overplot)

    p2 = plot([SFR_percentiles[1], SFR_percentiles[1]],$
              yrange,$
              linestyle='-',$
              color='orange',$
              thick=2,$
              /overplot)

    p2 = plot([SFR_percentiles[2], SFR_percentiles[2]],$
              yrange,$
              linestyle='--',$
              color='orange',$
              thick=2,$
              /overplot)

    p2.xrange = [0, p2.xrange[1]]

    p2.position = [0.1, 0.1, 0.1 + hist_width, 0.1 + hist_width]
    p2.xtitle = 'SFR [$M_{\odot}\ yr^{-1}$]'
    p2.ytitle = 'PDF'

    ; We have to do some weird stuff to get the SFH plotted nicely
    ; with a connected stairstep plot
    bounds = ([age_bin_edges, age_bin_edges])[sort([age_bin_edges,age_bin_edges])]
    bounds = bounds[1:-2]
    bounds[0] = (bounds[1])/10.d0

    Nsteps = n_elements(age_bin_edges) - 1
    SFH_plot = dblarr(Nsteps * 2, 3)
    for i=0, Nsteps - 1 do begin
        SFH_plot[(2*i):(2*i)+1, *] = rebin(SFH_percentiles[i, *], 2, 3)
    endfor

    ; Plot the uncertainty range
    p3 = Fillplot(bounds,$
                  [[SFH_plot[*, 0]], [SFH_plot[*, 2]]],$
                  color='gray',$
                  /current,$
                  /xlog)

    ; Plot the SFH on top
    p3 = plot(bounds,$
              SFH_plot[*, 1],$
              /overplot,$
              color='black',$
              thick=2,$
              /xlog)

    p3.xrange = [bounds[0], bounds[-1]]
    p3.position = [0.55, 0.1, 0.9, 0.9]

    p3.xtitle = 'Lookback Time $t$ [yr]'
    p3.ytitle = 'SFR(t) [$M_{\odot}\ yr^{-1}$]'

    return, p3

end
