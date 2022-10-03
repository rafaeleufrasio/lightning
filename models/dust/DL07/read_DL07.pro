Umin_grid = []
Umin_grid = [Umin_grid, '0.10'] ; Umax=[ '0.10', 1e2 ,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '0.15'] ; Umax=[ '0.15', 1e2 ,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '0.20'] ; Umax=[ '0.20', 1e2 ,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '0.30'] ; Umax=[ '0.30',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '0.40'] ; Umax=[ '0.40',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '0.50'] ; Umax=[ '0.50',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '0.70'] ; Umax=[ '0.70',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '0.80'] ; Umax=[ '0.80',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '1.00'] ; Umax=[ '1.00',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '1.20'] ; Umax=[ '1.20',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '1.50'] ; Umax=[ '1.50',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '2.00'] ; Umax=[        ?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '2.50'] ; Umax=[ '2.50',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '3.00'] ; Umax=[ '3.00',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '4.00'] ; Umax=[ '4.00',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '5.00'] ; Umax=[ '5.00',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '7.00'] ; Umax=[ '7.00',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid, '8.00'] ; Umax=[ '8.00',?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid,'12.0' ] ; Umax=['12.0' ,?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid,'15.0' ] ; Umax=['15.0' ,?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid,'20.0' ] ; Umax=['20.0' ,?1e2?,1e3,1e4,1e5,1e6]
Umin_grid = [Umin_grid,'25.0' ] ; Umax=['25.0' ,?1e2?,1e3,1e4,1e5,1e6]
;Umin_grid = [Umin_grid,'1e2'  ] ; Umax=[1e2]
;Umin_grid = [Umin_grid,'3e2'  ] ; Umax=[3e2]
;Umin_grid = [Umin_grid,'1e3'  ] ; Umax=[1e3]
;Umin_grid = [Umin_grid,'3e3'  ] ; Umax=[3e3]
;Umin_grid = [Umin_grid,'1e4'  ] ; Umax=[1e4]
;Umin_grid = [Umin_grid,'3e4'  ] ; Umax=[3e4]
;Umin_grid = [Umin_grid,'1e5'  ] ; Umax=[1e5]
;Umin_grid = [Umin_grid,'3e5'  ] ; Umax=[3e5]


Umax_grid = []
;Umax_grid = [Umax_grid, '1e2']
Umax_grid = [Umax_grid, '1e3']
Umax_grid = [Umax_grid, '1e4']
Umax_grid = [Umax_grid, '1e5']
Umax_grid = [Umax_grid, '1e6']

modl_grid = []
modl_grid = [modl_grid,'MW3.1_00'] ;qPAH = 0.47 %
modl_grid = [modl_grid,'MW3.1_10'] ;qPAH = 1.12 %
modl_grid = [modl_grid,'MW3.1_20'] ;qPAH = 1.77 %
modl_grid = [modl_grid,'MW3.1_30'] ;qPAH = 2.50 %
modl_grid = [modl_grid,'MW3.1_40'] ;qPAH = 3.19 %
modl_grid = [modl_grid,'MW3.1_50'] ;qPAH = 3.90 %
modl_grid = [modl_grid,'MW3.1_60'] ;qPAH = 4.58 %
modl_grid = [modl_grid,'LMC2_00']  ;qPAH = 0.75 %
modl_grid = [modl_grid,'LMC2_05']  ;qPAH = 1.49 %
modl_grid = [modl_grid,'LMC2_10']  ;qPAH = 2.37 %
modl_grid = [modl_grid,'smc']      ;qPAH = 0.10 %


;  jm     name      q_PAH
;  --   --------    ------
;   1   MW3.1_00    0.47 %
;   2   MW3.1_10    1.12 %
;   3   MW3.1_20    1.77 %
;   4   MW3.1_30    2.50 %
;   5   MW3.1_40    3.19 %
;   6   MW3.1_50    3.90 %
;   7   MW3.1_60    4.58 %
;   8   LMC2_00     0.75 %
;   9   LMC2_05     1.49 %
;  10   LMC2_10     2.37 %
;  11   smc         0.10 %


; IF ALPHA = 2 :
; format: '~/Dropbox/DL07/DL07spec/'+Umin+'/U'+Umin+'_'+Umax+'_'+modl+'.txt'
readcol,'~/Dropbox/DL07/DL07spec/U0.10/U0.10_0.10_MW3.1_00.txt',lambda,nuLnu,jnu
l0 = where(lambda eq 1e4,/null)
lambda_min=lambda[l0:*]                                                             
nuLnu_min=nuLnu[l0:*]  
jnu_min=jnu[l0:*]      
readcol,'~/Dropbox/DL07/DL07spec/U0.10/U0.10_1e4_MW3.1_00.txt',lambda,nuLnu,jnu
l0 = where(lambda eq 1e4,/null)
lambda_pow=lambda[l0:*]                                                             
nuLnu_pow=nuLnu[l0:*]  
jnu_pow=jnu[l0:*]

device,decompose=0
loadct,39         

plot_oo,lambda_pow,nuLnu_pow
  for i=0,10 do oplot,lambda_pow,0.1*i*nuLnu_min + (1-0.1*i)*nuLnu_pow,color=54+200*i/10
  for i=0,10 do oplot,lambda_pow,0.1*i*nuLnu_min,color=54+200*i/10,linestyle=1
  for i=0,10 do oplot,lambda_pow,(1-0.1*i)*nuLnu_pow,color=54+200*i/10,linestyle=2




foreach umin, Umin_grid do begin $
  foreach umax, Umax_grid do begin $
    foreach modl, modl_grid do begin $
      readcol,'~/Dropbox/DL07/DL07spec/U'+Umin+'/U'+Umin+'_'+Umin+'_'+modl+'.txt',lambda_min,nuLnu_min,jnu_min,/silent &$
      l0 = where(lambda_min eq 1e4,/null) &$
      readcol,'~/Dropbox/DL07/DL07spec/U'+Umin+'/U'+Umin+'_'+Umax+'_'+modl+'.txt',lambda_pow,nuLnu_pow,jnu_pow,/silent &$
      l1 = where(lambda_pow eq 1e4,/null) &$
      print,umin,umax,modl,l0,l1,form='(3A10,I10,I10)' &$
      for i=0,10 do oplot,lambda_pow[l1:*],0.1*i*nuLnu_min[l0:*] + (1-0.1*i)*nuLnu_pow[l1:*],color=54+200*i/10



;Umin = 0.10
;U = [0.1, 1.e2, 1.e3, 1.e4, 1.e5]
Ustr = ['0.10','0.15','0.20','0.30','0.40','0.50','0.70','0.80','1.00','1.20','1.50', $
        '2.50','3.00','4.00','5.00','7.00','8.00','12.0','15.0','20.0','25.0',        $
        '1e2','3e2','1e3','3e3','1e4','3e4','1e5','3e5']

Udbl = double(Ustr)

n_U = Ustr.length
lambda = dblarr(n_U,1001)
; => lambdas are all exactly the same for all delta functions
nuLnu = dblarr(n_U,1001)

for k=0,(n_U-1) do begin                          $
  U = Ustr[k]                                    &$
  print,k,'  ',U                                 &$
  readcol,'~/Dropbox/DL07/DL07spec/U'+U+'/U'+U+'_'+U+'_MW3.1_00.txt',lambda_min,nuLnu_min,jnu_min,/silent &$
  l0 = ((where(lambda_min eq 1e4,/null))[-1])[0] &$
  lambda = lambda_min[l0:*]                      &$
  nuLnu[k,*] = nuLnu_min[l0:*]


readcol,'~/Dropbox/DL07/DL07spec/U0.10/U0.10_1e5_MW3.1_00.txt',lambda_pow,nuLnu_pow,jnu_pow,/silent
;readcol,'~/Dropbox/DL07/DL07spec/U0.10/U0.10_1e4_MW3.1_00.txt',lambda_pow,nuLnu_pow,jnu_pow,/silent
;readcol,'~/Dropbox/DL07/DL07spec/U0.10/U0.10_1e3_MW3.1_00.txt',lambda_pow,nuLnu_pow,jnu_pow,/silent
l0 = ((where(lambda_pow eq 1e4,/null))[-1])[0]
lambda_pow = lambda_pow[l0:*]
nuLnu_pow = nuLnu_pow[l0:*]


U = Udbl
Umin = 0.1
Umax = 1e4

alpha = 2
fU = U^(-alpha)
fU /= trap_int(U,fU,xrange=[Umin,Umax],/pow)

nuLnu_pow_stitched = dblarr(1001)
for k=0,lambda.length-1 do nuLnu_pow_stitched[k] = trap_int(U,nuLnu[*,k]*fU,xrange=[Umin,Umax],/pow)


plot_oo,lambda,nuLnu[-1,*]
  for k=0,(n_U-1) do oplot,lambda,nuLnu[k,*]
  oplot,lambda_pow,nuLnu_pow,color=54,thick=3
  oplot,lambda,nuLnu_pow_stitched,color=254,linestyle=2,thick=3


nu = 1.d4* !cv.clight/lambda


print,minmax(nuLnu_pow/nuLnu_pow_stitched),trap_int(nu,nuLnu_pow/nu),trap_int(nu,nuLnu_pow_stitched/nu)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; alpha=2, Umin=0.1, Umax=1e5, MW3.1_00
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;
;/pow
;;;;;;;;
;      0.82369260       1.3171687
;   6.1682818e-24   6.1298851e-24
;;;;;;;;
;/exp
;;;;;;;;
;      0.78943574       1.2675572
;   6.1682818e-24   6.4687332e-24
;;;;;;;;
;LINEAR
;;;;;;;;
;      0.78859554       1.2609870
;   6.1682818e-24   6.7134078e-24
;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; alpha=2, Umin=0.1, Umax=1e4, MW3.1_00
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;
;/pow
;;;;;;;;
;      0.75564881       1.3334085
;   5.1428762e-24   5.1142739e-24
;;;;;;;;
;/exp
;;;;;;;;
;      0.72041441       1.2756716
;   5.1428762e-24   5.3594158e-24
;;;;;;;;
;LINEAR
;;;;;;;;
;      0.71151744       1.2530990
;   5.1428762e-24   5.5234496e-24
;;;;;;;;




