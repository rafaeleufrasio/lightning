; Load the postprocessed results
res = mrdfits('lightning_output/GOODSN_normal_galaxies.fits.gz',1)

print, '//Convergence for the first galaxy fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(res[0].ACCEPTANCE_FRAC), 2)
print, 'Convergence flag: ', strtrim(res[0].CONVERGENCE_FLAG, 2)
print, 'Short chain flag: ', strtrim(res[0].SHORT_CHAIN_FLAG, 2)
print, 'Number of "stranded" walkers: ', strtrim(total(res[0].STRANDED_FLAG), 2)
print, 'PPC p-value: ', strtrim(res[0].PVALUE, 2)
print, ''
print, '//Convergence for the second galaxy fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(res[1].ACCEPTANCE_FRAC), 2)
print, 'Convergence flag: ', strtrim(res[1].CONVERGENCE_FLAG, 2)
print, 'Short chain flag: ', strtrim(res[1].SHORT_CHAIN_FLAG, 2)
print, 'Number of "stranded" walkers: ', strtrim(total(res[1].STRANDED_FLAG), 2)
print, 'PPC p-value: ', strtrim(res[1].PVALUE, 2)

p = sed_plots(res)

p.save, 'GOODSN_normal_galaxies_SEDs.png', dpi=400
p.close
