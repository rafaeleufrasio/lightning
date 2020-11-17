function dl07_fPDR,dl07_alpha,dl07_Umin,dl07_Umax,dl07_gamma
;computes fPDR according to Eq. 4 of Aniano+2020
; written by Rafael Eufrasio in August 19 2020

hund_umax = 1d2       / dl07_Umax
umin_umax = dl07_Umin / dl07_Umax
alpha_minus_one = dl07_alpha - 1d0
two_minus_alpha = 2d0 - dl07_alpha
gam = dl07_gamma

N = dl07_gamma.length
fPDR = dblarr(N)
num = dblarr(N)
den = dblarr(N)

heq2 = where(dl07_alpha eq 2d0,Neq2,/null)
hne2 = where(dl07_alpha ne 2d0,Nne2,/null)

;for alpha = 2
IF Neq2 ne 0 THEN BEGIN
  num[heq2] = -gam[heq2]*alog(hund_umax[heq2])
  den[heq2] = (1d0 - gam[heq2])*(1d0 - Umin_Umax[heq2]) - gam[heq2]*alog(Umin_Umax[heq2])
ENDIF

;for alpha != 2
IF Nne2 ne 0 THEN BEGIN
  num[hne2] = gam[hne2] * (1d0 - (hund_umax[hne2])^two_minus_alpha[hne2])
  den[hne2] = (1d0-gam[hne2]) * two_minus_alpha[hne2]/alpha_minus_one[hne2]   $
                              * (Umin_Umax[hne2]^two_minus_alpha[hne2])       $
                              * (1d0 - umin_umax[hne2]^alpha_minus_one[hne2]) $
               + gam[hne2]    * (1d0 - umin_umax[hne2]^two_minus_alpha[hne2])
ENDIF

return,num/den

end