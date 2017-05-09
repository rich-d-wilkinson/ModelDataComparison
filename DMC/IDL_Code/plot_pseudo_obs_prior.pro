pro plot_pseudo_obs_prior, expt, pngorps, leg_plot

; Plot results from expts with various priors using full grid as a base dataset and pseudo-observations
; Setup plotting
plot_setup, pngorps, scale, charscale, linescale
leg=''
if (leg_plot) then leg = '_leg'

expt_group='tdgtg'
dir = 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\
model_dir = dir+'Model_data\'+expt_group+'\'
obs_dir = dir+'Observation_Data\'+expt_group+'\'
out_dir = dir+'Output\O-tdgtg_M-tdgtg_ME-Egg-0.05_SM-0.3\'

case expt of
'globav'  :  colour='SST_del2'
'latav'   :  colour='SST_del'
'cosine'  :  colour='SST_del'
'cosine2' :  colour='SST_del'
'anom'    :  colour='SST_del'
'abs'     :  colour='SST'
else      :  colour='SST_del'
endcase

; Read original datasets
; ----------------------------------
if (expt_group eq 'HadISST') then begin
   lon_model=findgen(360)-179.5
   lat_model=-1.0*(findgen(180)-89.5)
endif else begin
   lon_model=1.25*findgen(288)-180.
   lat_model=1.25*findgen(144)-89.375
endelse      

model_filename = model_dir+expt+'.txt'
get_lun, unit
openr, unit, model_filename
data=read_ascii(model_filename)
t_in=data.field1
close, unit
free_lun, unit

get_lun, unit
mask_filename = model_dir+'mask.txt'
openr, unit, mask_filename
data=read_ascii(mask_filename)
mask=data.field1
close, unit
free_lun, unit

temp=make_array(n_elements(lon_model),n_elements(lat_model),/float,value=-999)
; Need to write elements of t_in to appropriate place in temp array...
ilat=floor(mask/n_elements(lon_model))
ilon=mask-n_elements(lon_model)*ilat-1
temp(ilon,ilat)=t_in

; Read input points for data location plotting and comparison with output field
; ----------------------------------
get_lun, unit
openr, unit, obs_dir+'prism_'+expt+'.txt'
data=read_ascii(obs_dir+'prism_'+expt+'.txt', data_start=1)
close, unit
free_lun, unit
lon_obs =reform(data.field1(0,*))
lat_obs =reform(data.field1(1,*))
temp_obs =reform(data.field1(2,*))
sd_obs =reform(data.field1(3,*))

read_DMC_output, out_dir, expt, lon_grid, lat_grid, temp_grid, SD_grid, mask_grid

; Plot fit output


; Set plotting levels
get_levels, colour, levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, temp_grid, temp_grid_out, levels
if (colour eq 'SST') then begin
   map_colour_table, 99, ncol, 0, 0, "white"
endif else begin   
   ncol=n_elements(levels)
   ncol_side=ncol/2
   one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
endelse
; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename=out_dir+expt+'_SST_fit'+leg
title='Fit Field SST!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, temp_grid_out[1:n_elements(lon_grid)-1,*], lon_grid[1:n_elements(lon_grid)-1], lat_grid, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  legend=leg_plot, /map
for i=0,n_elements(lon_obs)-1 do begin
   xx=lon_obs(i) & yy=lat_obs(i)
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
bin_data, SD_grid, SD_grid_out, levels
map_colour_table, 99, ncol, 0, 0, "white"

;Create filename stem for plot
filename=out_dir+expt+'_SST_SD_fit'+leg
title='Fit Field SD!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, SD_grid_out[1:n_elements(lon_grid)-1,*], lon_grid[1:n_elements(lon_grid)-1], lat_grid, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  legend=leg_plot, /map
for i=0,n_elements(lon_obs)-1 do begin
   xx=lon_obs(i) & yy=lat_obs(i)
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

get_levels, 'lsm_frac', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, mask_grid, mask_grid_out, levels
map_colour_table, 99, ncol, 0, 0, "white"
filename=out_dir+expt+'_mask_fit'+leg
title='Mask!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
plot_open, pngorps, filename
map_discontinuous, mask_grid_out[1:n_elements(lon_grid)-1,*], lon_grid[1:n_elements(lon_grid)-1], lat_grid, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  legend=leg_plot, /map
for i=0,n_elements(lon_obs)-1 do begin
   xx=lon_obs(i) & yy=lat_obs(i)
         plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor
xyouts, 0.5, 0.85-leg_plot*0.05, 'Mask has been interpolated onto an FE mesh and back again!', alignment=0.5, charsize=1.0*charscale, color=fsc_color("Black"),/normal
plot_close, pngorps, filename

; Regrid fitted field onto model grid for comparison with original field
; ------------------------------------------------------------------------
; Find array positions in fitted field of model grid locations
lon_count=findgen(n_elements(lon_grid))
lat_count=findgen(n_elements(lat_grid))
lon_interp=interpol(lon_count,lon_grid,lon_model)
lat_interp=interpol(lat_count,lat_grid,lat_model)

temp_grid_model = interpolate(temp_grid, lon_interp, lat_interp, /GRID)
SD_grid_model = interpolate(SD_grid, lon_interp, lat_interp, /GRID)

; Set plotting levels
get_levels, colour, levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, temp_grid_model, temp_grid_model_out, levels

if (colour eq 'SST') then begin
   map_colour_table, 99, ncol, 0, 0, "white"
endif else begin   
   ncol=n_elements(levels)
   ncol_side=ncol/2
   one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
endelse

; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename=out_dir+expt+'_SST_fit_model'+leg
title='Fit Field SST: Input Resolution!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, temp_grid_model_out, lon_model, lat_model, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  legend=leg_plot, /map
for i=0,n_elements(lon_obs)-1 do begin
   xx=lon_obs(i) & yy=lat_obs(i)
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
bin_data, SD_grid_model, SD_grid_model_out, levels
map_colour_table, 99, ncol, 0, 0, "white"


; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename=out_dir+expt+'_SST_SD_fit_model'+leg
title='Fit Field SD: Input Resolution!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, SD_grid_model_out, lon_model, lat_model, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  legend=leg_plot, /map
for i=0,n_elements(lon_obs)-1 do begin
   xx=lon_obs(i) & yy=lat_obs(i)
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

; Now can find difference between model and kriged fields
diff=temp_grid_model-temp
diff(where(temp lt -100)) = -999.
diff_norm=(temp_grid_model-temp)/sd_grid_model
diff_norm(where(temp lt -100)) = -999.

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
filename=out_dir+expt+'_SST_fit_del'+leg
title='SST anomaly: Fit Field-Input!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   
map_discontinuous, diff_out, lon_model, lat_model, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  legend=leg_plot, /map
for i=0,n_elements(lon_obs)-1 do begin
   xx=lon_obs(i) & yy=lat_obs(i)
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

get_levels, 'SST_SD_del', levels, ncol, format, ticknames
one_two_tone_MO, 2, 'blue', ncol/2, 'red', ncol/2, 'white'
bin_data, diff_norm, diff_norm_out, levels
filename=out_dir+expt+'_SST_fit_del_norm'+leg
title='SST anomaly: (Fit Field-Input)/Fit SD!c'+expt
xtitle='lon'
ytitle='lat'
extra_leg=0
plot_open, pngorps, filename
map_discontinuous, diff_norm_out, lon_model, lat_model, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  legend=leg_plot, /map
for i=0,n_elements(lon_obs)-1 do begin
   xx=lon_obs(i) & yy=lat_obs(i)
         plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor
plot_close, pngorps, filename

; Find how closely fitted field matches input data points, absolute and /SD
; Sample output field at PRISM points and compare with input

x = VALUE_LOCATE(lon_grid-((lon_grid[1]-lon_grid[0])/2),lon_obs) 
y = VALUE_LOCATE(lat_grid-((lat_grid[1]-lat_grid[0])/2),lat_obs)
temp_grid_obs=temp_grid[x,y]
fit_del=temp_grid_obs-temp_obs
fit_del_norm=(temp_grid_obs-temp_obs)/sd_obs
SD_grid_obs=SD_grid[x,y]
fit_del_norm2=(temp_grid_obs-temp_obs)/sd_grid_obs
SD_grid_obs_norm=SD_grid_obs/sd_obs

get_levels, 'SST_SD_del',  levels, ncol, fmt, ticknames
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
bin_data, fit_del, fit_del_out, levels
bin_data, fit_del_norm, fit_del_norm_out, levels
bin_data, fit_del_norm2, fit_del_norm2_out, levels

filename = out_dir+expt+'_SST_fit_del_points'+leg
plot_open, pngorps, filename
map_points,    fit_del_out, lon_obs, lat_obs, ncol, 'Fit Field-Input Points!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
plot_close, pngorps, filename
print, 'Fit field - input: min/max ',min(fit_del), max(fit_del)

filename = out_dir+expt+'_SST_fit_del_norm_points'+leg
plot_open, pngorps, filename
map_points,    fit_del_norm_out, lon_obs, lat_obs, ncol, '(Fit Field-Input Points)/Input SD!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
plot_close, pngorps, filename
print, 'Fit field - input/input SD: min/max ',min(fit_del_norm), max(fit_del_norm)

filename = out_dir+expt+'_SST_fit_del_norm2_points'+leg
plot_open, pngorps, filename
map_points,    fit_del_norm2_out, lon_obs, lat_obs, ncol, '(Fit Field-Input Points)/Fit SD!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
plot_close, pngorps, filename
print, 'Fit field - input/fit SD: min/max ',min(fit_del_norm2), max(fit_del_norm2)

get_levels, 'SST_SD', levels, ncol, format, ticknames
map_colour_table, 99, ncol, 0, 0, "white"
bin_data, SD_grid_obs, SD_grid_obs_out, levels
bin_data, SD_grid_obs_norm, SD_grid_obs_norm_out, levels

filename = out_dir+expt+'_SD_fit_points'+leg
plot_open, pngorps, filename
map_points,    SD_grid_obs_out, lon_obs, lat_obs, ncol, 'Fit SD!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
plot_close, pngorps, filename
print, 'Fit SD: min/max ',min(SD_grid_obs), max(SD_grid_obs)

filename = out_dir+expt+'_SD_fit_norm_points'+leg
plot_open, pngorps, filename
map_points,    SD_grid_obs_norm_out, lon_obs, lat_obs, ncol, 'Fit SD/Input SD!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
plot_close, pngorps, filename
print, 'Fit SD/ input SD: min/max ',min(SD_grid_obs_norm), max(SD_grid_obs_norm)


return
end