pro plot_krig_FE_data, output, pngorps

dir = 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\'
output_dir = dir+'Output\'+output+path_sep()
obs = file_basename(file_search(output_dir,'*_grid_output.txt'),'_grid_output.txt')
print, 'Files found: ', obs
dummy = strsplit(output,'O-', /regex, /extract)
dummy = strsplit(dummy[0],'_M-', /regex, /extract)
obs_set=dummy[0]
; this could be egg of sphere, need to test for this in some way
dummy = strsplit(dummy[1],'_Me-', /regex, /extract)
model_set=dummy[0]
obs_dir = dir+'Observation_data'+path_sep()+obs_set+path_sep()
model_dir = dir+'Model_data'+path_sep()+model_set+path_sep()

; Setup plotting
if (pngorps eq 'png') then begin
   set_plot, 'win'
   device, decomposed=0
endif else begin
   set_plot, 'ps'
endelse
plot_setup, pngorps, scale, charscale, linescale
device, set_font='helvetica', /tt_font

for i=0,n_elements(obs)-1 do begin
   ; Read observation data to mark on plots
   filename=obs_dir+obs[i]+'.txt'
   get_lun, unit
   openr, unit, filename
   data=read_ascii(filename,record_start=1)
   close, unit
   free_lun, unit
   
   lon_obs =reform(data.field1(0,*))
   lat_obs =reform(data.field1(1,*))
   temp_obs =reform(data.field1(2,*))
   SD_obs =reform(data.field1(3,*))
   
   
   ; Read FE gridded output of krigged field
   ; FE 2 grid conversion runs from min to max lon and lat found in output and splits evenly into 'detail' points
   ; Find min and max and n_elements for lon/lat arrays, then can convert to indices to create mean/std/mask arrays
   
   ; Read mean, std and mask from text file
   filename = output_dir+obs(i)+'_grid_output.txt'
   get_lun, unit
   openr, unit, filename
   data=read_ascii(filename, missing_value=-999, data_start=1)
   close, unit
   free_lun, unit
   
   lon_fit  = reform(data.field01(1,*)*180./!pi)
   lat_fit  = reform(data.field01(2,*)*180./!pi)
   temp_fit = reform(data.field01(3,*))
   SD_fit   = reform(data.field01(6,*))
   mask_fit = reform(data.field01(9,*))
   
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
   mask_grid=make_array(n_elements(lons),n_elements(lats),/float,value=-9999)
   mask_grid(lon_ind,lat_ind)=mask_fit
   
   ; Plot fit output
   
   
   ; Set plotting levels
   if (strmatch(obs_set,'*_anom')) then begin
      get_levels, 'SST_del', levels, ncol, format, ticknames
      one_two_tone_MO, 2, 'blue', ncol/2, 'red', ncol/2, 'white'
   endif else begin   
      get_levels, 'SST', levels, ncol, format, ticknames
      map_colour_table, 99, ncol, 0, 0, "white"
   endelse   
   ; Bin data for each dataset if gridded data required
   bin_data, temp_grid, temp_grid_out, levels
   bin_data, temp_obs, temp_obs_out, levels
   
   ncol=n_elements(levels)
   ncol_side=ncol/2
   
   
   ; FOR EACH PLOT REQUIRED
   ; CALL SPEICALIST PLOTTING SUBROUTINE
   ; Set plot attributes
   ;Create filename stem for plot
   filename=output_dir+obs(i)+'_SST_obs_fit'
   title='SST Observations - Kriged Field!c'+obs(i)
   xtitle='lon'
   ytitle='lat'
   extra_leg=0
   if (pngorps eq 'ps') then begin
      device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
   endif   
   
   map_discontinuous, temp_grid_out[1:n_elements(lons)-1,*], lons[1:n_elements(lons)-1], lats, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
   ;contour, mask_grid, lons, lats, levels=0.5,  color=fsc_color("Black"), thick=1.5, /isotropic, /overplot
   ; Plots circles to show data locations, with original values
   for m=0,n_elements(lon_obs)-1 do begin
      xx=lon_obs(m) & yy=lat_obs(m)
            polyfill, circle(xx,yy,2,1),/data,color=temp_obs_out(m)
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
   bin_data, SD_obs, SD_obs_out, levels
   map_colour_table, 99, ncol, 0, 0, "white"
   
   
   ; FOR EACH PLOT REQUIRED
   ; CALL SPEICALIST PLOTTING SUBROUTINE
   ; Set plot attributes
   ;Create filename stem for plot
   filename=output_dir+obs(i)+'_Uncertainty_fit'
   title='Uncertainty - Kriged Field!c'+obs(i)
   xtitle='lon'
   ytitle='lat'
   extra_leg=0
   if (pngorps eq 'ps') then begin
      device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
   endif   
   
   map_discontinuous, SD_grid_out[1:n_elements(lons)-1,*], lons[1:n_elements(lons)-1], lats, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
   ;contour, mask_grid, lons, lats, levels=0.5,  color=fsc_color("Black"), thick=1.5, /isotropic, /overplot
   for m=0,n_elements(lon_obs)-1 do begin
      xx=lon_obs(m) & yy=lat_obs(m)
            ;polyfill, circle(xx,yy,2,1),/data,color=SD_prism_out(m)
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
   get_levels, 'lsm_frac', levels, ncol, format, ticknames
   ; Bin data for each dataset if gridded data required
   bin_data, mask_grid, mask_grid_out, levels
   map_colour_table, 99, ncol, 0, 0, "white"
   
   
   ; FOR EACH PLOT REQUIRED
   ; CALL SPEICALIST PLOTTING SUBROUTINE
   ; Set plot attributes
   ;Create filename stem for plot
   filename=output_dir+obs(i)+'_Mask_fit'
   title='Mask - Kriged Field!c'+obs(i)
   xtitle='lon'
   ytitle='lat'
   extra_leg=0
   if (pngorps eq 'ps') then begin
      device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
   endif   
   
   map_discontinuous, mask_grid_out[1:n_elements(lons)-1,*], lons[1:n_elements(lons)-1], lats, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
   for m=0,n_elements(lon_obs)-1 do begin
      xx=lon_obs(m) & yy=lat_obs(m)
            ;polyfill, circle(xx,yy,1,1),/data,color=SD_prism_out(m)
            plots,    circle(xx,yy,2,1),/data,color=fsc_color("Black")
   endfor
   xyouts, 0.5, 0.8, 'Mask has been interpolated onto an FE grid and back again!', alignment=0.5, charsize=1.5*charscale, color=fsc_color("Black"),/normal
   ;     Save plot
   if (pngorps eq 'png') then begin
      write_png, filename+'.png', tvrd(true=1)
      wdelete
   endif else begin
      device, /close
   endelse
   ;     END PLOT
endfor

return
end