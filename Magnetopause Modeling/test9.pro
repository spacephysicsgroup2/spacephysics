;讀資料

openr,1,'test.txt'

n_l=file_lines('test.txt')

day=fltarr(n_l) & month=fltarr(n_l) & year=fltarr(n_l) & hour=fltarr(n_l) & minute=fltarr(n_l) & second=fltarr(n_l) & Bx_GSE=fltarr(n_l) & By_GSE=fltarr(n_l) & Bz_GSE=fltarr(n_l) & flow_speed=fltarr(n_l)

By_GSM=fltarr(n_l) & Bz_GSM=fltarr(n_l) & vx_GSE=fltarr(n_l) & vy_GSE=fltarr(n_l) & vz_GSE=fltarr(n_l) & np=fltarr(n_l) & T=fltarr(n_l) & Dp=fltarr(n_l) 

m=0L

while ~eof(1) do begin

  readf,1,dy,mo,yr,hr,mi,sc,Bx_e,By_e,Bz_e,By_m,Bz_m,flow_s,vx,vy,vz,n,Ti,Pi,format='(G2,1X,G2,1X,G4,1X,G2,1X,G2,1X,G6,G14,G14,G14,G14,G14,8X,G8,8X,G9,8X,G9,8X,G9,8X,A7,3X,G11,6X,G8)'

  n=float(n) 

  day[m]=dy & month[m]=mo & year[m]=yr & hour[m]=hr & minute[m]=mi & second[m]=sc & Bx_GSE[m]=Bx_e & By_GSE[m]=By_e & Bz_GSE[m]=Bz_e & By_GSM[m]=By_m & Bz_GSM[m]=Bz_m & flow_speed[m]=flow_s

  vx_GSE[m]=vx & vy_GSE[m]=vy & vz_GSE[m]=vz & np[m]=n & T[m]=Ti & Dp[m]=Pi 

  m=m+1L

endwhile

close,1

B=sqrt(Bx_GSE^2+By_GSE^2+Bz_GSE^2)

Bx_GSM=sqrt(B^2-By_GSM^2-Bz_GSM^2)

;設置繪圖環境

device,decomposed=0

!p.background=255 & !p.color=0

;分析磁層頂

mu_0=4*!pi*1e-7

cs_coeff=(5./3.)*1.38e-23/(0.84*1.66e-27)

theta=[0.:360.:10.]*!dtor

x_c=cos(theta) & z_c=sin(theta)

usersym,x_c,z_c,/fill

window,xsize=800.,ysize=800.

Bz_GSM_csv=fltarr(1)

Dp_csv=fltarr(1)

x_pause_csv=fltarr(1)

z_pause_csv=fltarr(1)

time_csv_1=fltarr(1)

time_csv_2=fltarr(1)

x_Shue_csv=fltarr(1)

y_Shue_csv=fltarr(1)

for i=0,n_l-1 do begin
  
  geopack_recalc,year[i],month[i],day[i],hour[i],minute[i],second[i],/date,tilt=d_tilt_angle
  
  plot,[0.,1.],/nodata,xrange=[-25.,30.],yrange=[-25.,30.],xstyle=1,ystyle=1,xtitle='x(Re)',ytitle='z(Re)',title='Magnetopause :'+string(year[i],format='(I4.4)')+'/'+string(month[i],format='(I2.2)')+'/'+string(day[i],format='(I2.2)')+' '+string(hour[i],format='(I2.2)')+':'+string(minute[i],format='(I2.2)')+':'+string(second[i],format='(I2.2)')
  
  xyouts,-5.,28.,' Bz ='+string(Bz_GSM[i],format='(f6.3)')+' nT'+' ,Dp ='+string(Dp[i],format='(f6.3)')+' nPa'+' ,dipole tilt ='+string(d_tilt_angle,'(f7.3)')
  
  oplot,[14.,14.],[26.,26.],psym=8,symsize=0.4,color=fsc_color('green')
  
  xyouts,15.,26.,'Equllibrium Point'
  
  oplot,[9.,14.],[24.,24.],thick=2,color=fsc_color('green')
  
  xyouts,15.,24.,'Polynomial fitting'
  
  oplot,[9.,14.],[22.,22.],thick=2,color=fsc_color('blue')

  xyouts,15.,22.,'Shue et al. 1998 Model'
    
  z_max_array=fltarr(60)

  x_max_array=fltarr(60)  
  
  z_min_array=fltarr(60)

  x_min_array=fltarr(60)

  x_right_array=fltarr(60)
  
  z_right_array=fltarr(60)
    
  for j=1,60 do begin
    
    geopack_sphcar,1.,j*3.,180.,xi,yi,zi,/to_rect,/degree
    
    parm=[Dp[i],0.,By_GSM[i],Bz_GSM[i],0.,0.,0.,0.,0.,0.]
    
    geopack_trace,xi,yi,zi,1,parm,xf,yf,zf,fline=coordinate,tilt=d_tilt_angle,/TS04
    
    x_m=coordinate[*,0] & z_m=coordinate[*,2]
    
    if n_elements(z_m[where(z_m gt 0.)]) ne 1 then begin
      
       oplot,x_m,z_m
      
    endif
    
    z_max_array[j-1]=max(z_m,max_index)
    
    x_max_array[j-1]=x_m[max_index]
    
    x_right_array[j-1]=max(x_m,right_index)

    z_right_array[j-1]=z_m[right_index]
     
  endfor
  
  oplot,x_c,z_c
  
  z_max_array2=z_max_array[where(z_max_array gt 0. and x_max_array lt 0.)]
  
  x_right_array2=x_right_array[where(x_right_array gt 0. and min(abs(z_right_array)) lt 1.)]
  
  j_max=where(z_max_array eq max(z_max_array2))+1
  
  j_right=where(x_right_array eq max(x_right_array2))+1
  
  ;Upper Boundary theta
  
  geopack_sphcar,1.,j_max*3.,180.,xi_max,yi_max,zi_max,/to_rect,/degree

  parm=[Dp[i],0.,By_GSM[i],Bz_GSM[i],0.,0.,0.,0.,0.,0.]

  geopack_trace,xi_max,yi_max,zi_max,1,parm,xf_max,yf_max,zf_max,fline=coordinate_max,tilt=d_tilt_angle,/TS04
  
  x_max=coordinate_max[*,0] & y_max=coordinate_max[*,1] & z_max=coordinate_max[*,2]

  z_max_index=where(z_max gt 1.)

  x_max=x_max[z_max_index] & y_max=y_max[z_max_index] & z_max=z_max[z_max_index]
  
  geopack_dip,x_max,y_max,z_max,Bx_dip,By_dip,Bz_dip,tilt=d_tilt_angle
  
  geopack_TS04,parm,x_max,y_max,z_max,Bx_max_T96,By_max_T96,Bz_max_T96,tilt=d_tilt_angle
  
  Bx_max_total=Bx_dip+Bx_max_T96 & By_max_total=By_dip+By_max_T96 & Bz_max_total=Bz_dip+Bz_max_T96
  
  theta_max=fltarr(n_elements(Bx_max_total))
  
  for k=0,n_elements(Bx_max_total)-1 do begin
    
    B_cross=crossp([Bx_GSM[i],By_GSM[i],Bz_GSM[i]],[Bx_max_total[k],By_max_total[k],Bz_max_total[k]])
    
    unit_normal=B_cross/norm(B_cross)
    
    theta_max[k]=180.-acos(unit_normal[0])*!radeg
    
  endfor
 
  ;Right Boundary theta
  
  geopack_sphcar,1.,j_right*3.,180.,xi_min,yi_min,zi_min,/to_rect,/degree

  parm=[Dp[i],0.,By_GSM[i],Bz_GSM[i],0.,0.,0.,0.,0.,0.]

  geopack_trace,xi_min,yi_min,zi_min,1,parm,xf_min,yf_min,zf_min,fline=coordinate_right,tilt=d_tilt_angle,/TS04

  x_right=coordinate_right[*,0] & y_right=coordinate_right[*,1] & z_right=coordinate_right[*,2]

  x_right_index=where(x_right gt 1.)

  x_right=x_right[x_right_index] & y_right=y_right[x_right_index] & z_right=z_right[x_right_index]

  geopack_dip,x_right,y_right,z_right,Bx_dip,By_dip,Bz_dip,tilt=d_tilt_angle

  geopack_TS04,parm,x_right,y_right,z_right,Bx_right_T96,By_right_T96,Bz_right_T96,tilt=d_tilt_angle

  Bx_right_total=Bx_dip+Bx_right_T96 & By_right_total=By_dip+By_right_T96 & Bz_right_total=Bz_dip+Bz_right_T96

  theta_right=fltarr(n_elements(Bx_right_total))

  for k=0,n_elements(Bx_right_total)-1 do begin

    B_cross=crossp([Bx_GSM[i],By_GSM[i],Bz_GSM[i]],[Bx_right_total[k],By_right_total[k],Bz_right_total[k]])

    unit_normal=B_cross/norm(B_cross)

    theta_right[k]=180.-acos(unit_normal[0])*!radeg
    
  endfor
  
  ;Newtonian Approximation
  
  cs=sqrt(cs_coeff*T[i])*1e-3
  
  Mach=flow_speed[i]/cs
  
  ga=5./3.
  
  K_coeff=((ga+1.)/2)^(ga^2-1)*1./(ga*(ga-(ga-1.)/(2*Mach^2))^(1./(ga-1)))-Np[i]*(1e6)*(1.38e-23)*(1e-9)*T[i]/Dp[i]
  
  ;找到壓力平衡處:Upper Boundary
 
  x_limit=x_max[where(x_max gt -25.)]
  
  z_limit=z_max[where(x_max gt -25.)]
  
  x_equi=x_limit 
  
  z_equi=fltarr(n_elements(x_equi))
  
  p_difference_min=fltarr(n_elements(x_equi))
  
  P_sh=K_coeff*Dp[i]*(cos(theta_max*!dtor))^2+Np[i]*(1e6)*(1.38e-23)*T[i]*(1e-9)
  
  P_Sh_limit=P_sh[where(x_max gt -25.)]
 
  z_range=[0.:ceil(max(z_max))+20.:0.1]
  
  for k=0,n_elements(x_limit)-1 do begin
    
    x_range=replicate(x_limit[k],n_elements(z_range)) & y_range=fltarr(n_elements(z_range))
    
    geopack_dip,x_range,y_range,z_range,Bx_dip,By_dip,Bz_dip,tilt=d_tilt_angle

    geopack_TS04,parm,x_range,y_range,z_range,Bx_range_T96,By_range_T96,Bz_range_T96,tilt=d_tilt_angle

    Bx_range_total=Bx_dip+Bx_range_T96 & By_range_total=By_dip+By_range_T96 & Bz_range_total=Bz_dip+Bz_range_T96

    B_range_total=sqrt(Bx_range_total^2+By_range_total^2+Bz_range_total^2)
    
    P_B_range=((B_range_total*1e-9)^2)/(2*mu_0)*1e9
    
    P_difference=P_B_range-P_sh_limit[k]
    
    if min(abs(p_difference)) lt 0.1 then begin
      
      p_difference_min[k]=min(abs(p_difference),min_range_index)
      
      z_equi[k]=z_range[min_range_index]
      
    endif
    
  endfor
  
  equi_index=where(z_equi ne 0. and z_equi-mean(z_equi) lt 10. and z_equi gt 5.)
  
  z_equi=z_equi[equi_index] & x_equi=x_equi[equi_index]
  
  x_pause_upper=x_equi & z_pause_upper=z_equi
  
  ;找到壓力平衡處:Right Upper Boundary

  x_limit=x_right[where(x_right gt 1. and z_right gt 0.)]
  
  z_right_upper=z_right[where(x_right gt 1. and z_right gt 0.)]

  x_equi=x_limit

  z_equi=fltarr(n_elements(x_equi))

  p_difference_min=fltarr(n_elements(x_equi))

  P_sh=K_coeff*Dp[i]*(cos(theta_right*!dtor))^2+Np[i]*(1e6)*(1.38e-23)*T[i]*(1e-9)

  P_Sh_limit=P_sh[where(x_right gt 1. and z_right gt 0.)]

  z_range=[0.:ceil(max(z_right_upper))+30.:0.1]

  for k=0,n_elements(x_limit)-1 do begin

    x_range=replicate(x_limit[k],n_elements(z_range)) & y_range=fltarr(n_elements(z_range))

    geopack_dip,x_range,y_range,z_range,Bx_dip,By_dip,Bz_dip,tilt=d_tilt_angle

    geopack_TS04,parm,x_range,y_range,z_range,Bx_range_T96,By_range_T96,Bz_range_T96,tilt=d_tilt_angle

    Bx_range_total=Bx_dip+Bx_range_T96 & By_range_total=By_dip+By_range_T96 & Bz_range_total=Bz_dip+Bz_range_T96

    B_range_total=sqrt(Bx_range_total^2+By_range_total^2+Bz_range_total^2)

    P_B_range=((B_range_total*1e-9)^2)/(2*mu_0)*1e9

    P_difference=P_B_range-P_sh_limit[k]

    if min(abs(p_difference)) lt 0.1 then begin

      p_difference_min[k]=min(abs(p_difference),min_range_index)

      z_equi[k]=z_range[min_range_index]

    endif

  endfor

  equi_index=where(z_equi ne 0. and z_equi-mean(z_equi) lt 10.)
  
  z_equi=z_equi[equi_index] & x_equi=x_equi[equi_index]
  
  x_pause_right_upper=x_equi & z_pause_right_upper=z_equi
  
  x_pause_right_lower=x_equi & z_pause_right_lower=-z_equi
  
  ;磁層頂資料點
  
  x_pause_lower=x_pause_upper & z_pause_lower=-z_pause_upper
  
  x_pause=[x_pause_lower,x_pause_right_lower,x_pause_right_upper,x_pause_upper] &  z_pause=[z_pause_lower,z_pause_right_lower,z_pause_right_upper,z_pause_upper]
  
  oplot,x_pause,z_pause,symsize=0.4,psym=8,color=fsc_color('green')
  
  ;磁層頂擬合
  
  coeff=poly_fit(z_pause,x_pause,2)
  
  z_magnetopause=[-30.:30.:0.1]
  
  x_magnetopause=poly(z_magnetopause,coeff)
  
  oplot,x_magnetopause,z_magnetopause,thick=2,color=fsc_color('green')
  
  ;Shue et al. 1998 Model
  
  r0=(10.22+1.29*tanh(0.184*(Bz_GSM[i]+8.14)))*Dp[i]^(-1/6.6)
  
  alpha=(0.58-0.007*Bz_GSM[i])*(1+0.024*alog(Dp[i]))
  
  theta_Shue=[-30.:30.:0.1]
  
  r_Shue=r0*(2/(1+cos(theta_Shue)))^alpha
  
  oplot,r_Shue*cos(theta_Shue),r_Shue*sin(theta_Shue),thick=2,color=fsc_color('blue')
  
  ;Screenshot
  
  image1=tvrd(true=1)
  
  filename='theoritical_'+string(year[i],format='(I4.4)')+string(month[i],format='(I2.2)')+string(day[i],format='(I2.2)')+string(hour[i],format='(I2.2)')+string(minute[i],format='(I2.2)')+string(second[i],format='(I2.2)')+'_Magnetopause.png'
  
  z_magnetopause_min=min(abs(z_magnetopause),x_concavity_index)
  
  if x_magnetopause[x_concavity_index] gt 0. then write_png,filename,image1
  
  if x_magnetopause[x_concavity_index] gt 0. then write_png,'theoritical'+string(i)+'.png',image1
 
  time_mark_1=replicate(i+1,n_elements(x_pause))
  
  Bz_GSM_csv_1=replicate(Bz_GSM[i],n_elements(x_pause))
  
  Dp_csv_1=replicate(Dp[i],n_elements(x_pause))
  
  time_csv_1=[time_csv_1,time_mark_1]
  
  x_pause_csv=[x_pause_csv,x_pause]
  
  z_pause_csv=[z_pause_csv,z_pause]
  
  Bz_GSM_csv=[Bz_GSM_csv,Bz_GSM_csv_1]
  
  Dp_csv=[Dp_csv,Dp_csv_1]
  
  time_mark_2=replicate(i+1,n_elements(r_Shue))

  x_shue=r_Shue*cos(theta_Shue) & y_Shue=r_Shue*sin(theta_Shue)
  
  time_csv_2=[time_csv_2,time_mark_2]
  
  x_Shue_csv=[x_Shue_csv,x_Shue]
  
  y_Shue_csv=[y_Shue_csv,y_shue]
    
endfor

write_csv,'example1_theory.csv',time_csv_1,x_pause_csv,z_pause_csv,Bz_GSM_csv,Dp_csv

write_csv,'example2_theory.csv',time_csv_2,x_shue_csv,y_shue_csv

end