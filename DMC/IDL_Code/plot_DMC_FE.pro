pro plot_krig_FE, expt, pngorps

; Setup plotting
if (pngorps eq 'png') then begin
   set_plot, 'win'
   device, decomposed=0
endif else begin
   set_plot, 'ps'
endelse
plot_setup, pngorps, scale, charscale, linescale
device, set_font='helvetica', /tt_font

; Read PRISM3 locations to mark on plots
get_lun, unit
openr, unit, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\observation_data\'+expt+'_prism.txt'
data=read_ascii('C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\observation_data\'+expt+'_prism.txt',record_start=1)
close, unit
free_lun, unit

lon_prism =reform(data.field1(0,*))
lat_prism =reform(data.field1(1,*))
temp_prism =reform(data.field1(2,*))
SD_prism =reform(data.field1(3,*))

; Read original UM output

dir = 'O:\Documents\11_Projects\03_PhD_Pliocene_CO2\HadCM3_jobs\webpages\'
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','longitude',lon
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','latitude',lat
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','temp_mm_uo',temp
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfsdann.nc','temp_mm_uo',temp_sd
read_ncvar, dir+'tczyc'+path_sep()+'climate'+path_sep()+'tczyco.pfclann.nc','temp_mm_uo',temp_cont

; shift to -180/180 lons
ids = where(lon ge 180)
lon(ids) = lon(ids)-360.
lon = shift(lon,144)
temp = shift(temp,144,0)
temp_cont = shift(temp_cont,144,0)
temp_sd = shift(temp_sd,144,0)
temp_anom = temp-temp_cont


; Read FE gridded output of krigged field
; FE 2 grid conversion runs from min to max lon and lat found in output and splits evenly into 'detail' points
; Find min and max and n_elements for lon/lat arrays, then can convert to indices to create mean/std arrays

; Read mean and std from text file
get_lun, unit
openr, unit, "C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\"+expt+"_prism_grid.txt"
data=read_ascii("C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\"+expt+"_prism_grid.txt", missing_value=-999, data_start=1)
close, unit
free_lun, unit

lon_fit  = reform(data.field1(1,*)*180./!pi)
lat_fit  = reform(data.field1(2,*)*180./!pi)
temp_fit = reform(data.field1(3,*))
SD_fit   = reform(data.field1(6,*))

; Find unique values of lon and lat
lons = lon_fit(UNIQ(lon_fit, SORT(lon_fit)))
lats = lat_fit(UNIQ(lat_fit, SORT(lat_fit)))
; Convert lon_fit and lat_fit to indices for assigning temp and SD to arrays
lon_ind= round((lon_fit-min(lon_fit))/(lons[1]-lons[0]))
lat_ind= round((lat_fit-min(lat_fit))/(lats[1]-lats[0]))

; Generate temp array and populate
temp_grid=make_array(n_elements(lons),n_elements(lats),/float,value=-9999)
temp_grid(lon_ind,lat_ind)=temp_fit
SD_grid=make_array(n_elements(lons),n_elements(lats),/float,value=-9999)
SD_grid(lon_ind,lat_ind)=SD_fit

; Plot fit output


; Set plotting levels
get_levels, 'SST_del', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, temp_grid, temp_grid_out, levels
bin_data, temp_anom, temp_anom_out, levels
bin_data, temp_prism, temp_prism_out, levels

ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'

; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\'+expt+'_SST_anom_fit'
title='SST anomaly - Kriged Field!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, temp_grid_out[1:n_elements(lons)-1,*], lons[1:n_elements(lons)-1], lats, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map

; Plots circles to show data locations, with original values
for i=0,n_elements(lon_prism)-1 do begin
   xx=lon_prism(i) & yy=lat_prism(i)
         polyfill, circle(xx,yy,2,1),/data,color=temp_prism_out(i)
         plots,    circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor

;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse
;     END PLOT

; Set plotting levels
get_levels, 'SST_SD', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, SD_grid, SD_grid_out, levels
bin_data, SD_prism, SD_prism_out, levels
map_colour_table, 99, ncol, 0, 0, "white"


; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\'+expt+'_SST_SD_fit'
title='SST SD - Kriged Field!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, SD_grid_out[1:n_elements(lons)-1,*], lons[1:n_elements(lons)-1], lats, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
for i=0,n_elements(lon_prism)-1 do begin
   xx=lon_prism(i) & yy=lat_prism(i)
         polyfill, circle(xx,yy,2,1),/data,color=SD_prism_out(i)
         plots,    circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor

;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse
;     END PLOT

; Regrid kriged field onto UM grid for comparison
; Find array positions in kriged field of UM grid locations
lon_count=findgen(n_elements(lons))
lat_count=findgen(n_elements(lats))
lon_interp=interpol(lon_count,lons,lon)
lat_interp=interpol(lat_count,lats,lat)

temp_grid_um = interpolate(temp_grid, lon_interp, lat_interp, /GRID)
SD_grid_um = interpolate(SD_grid, lon_interp, lat_interp, /GRID)

; Set plotting levels
get_levels, 'SST_del', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, temp_grid_um, temp_grid_um_out, levels

ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'

; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\'+expt+'_SST_anom_fit_um'
title='SST anomaly - Kriged Field, UM Resolution!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, temp_grid_um_out, lon, lat, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
for i=0,n_elements(lon_prism)-1 do begin
   xx=lon_prism(i) & yy=lat_prism(i)
         plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor

;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse
;     END PLOT

; Set plotting levels
get_levels, 'SST_SD', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, SD_grid_um, SD_grid_um_out, levels
map_colour_table, 99, ncol, 0, 0, "white"


; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\'+expt+'_SST_SD_fit_um'
title='SST SD - Kriged Field, UM Resolution!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, SD_grid_um_out, lon, lat, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
for i=0,n_elements(lon_prism)-1 do begin
   xx=lon_prism(i) & yy=lat_prism(i)
         plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor

;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse
;     END PLOT

; Now can find difference between model and kriged version
diff=temp_grid_um-temp_anom
diff(where(temp gt 100)) = -999

; Set plotting levels
get_levels, 'SST_del', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, diff, diff_out, levels

ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'

; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\'+expt+'_SST_anom_fit_diff'
title='SST anomaly: UM-Kriged Field!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, diff_out, lon, lat, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
for i=0,n_elements(lon_prism)-1 do begin
   xx=lon_prism(i) & yy=lat_prism(i)
         plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor

;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse
;     END PLOT


return
end