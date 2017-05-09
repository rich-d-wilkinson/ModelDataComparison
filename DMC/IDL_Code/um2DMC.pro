pro um2krig, expt, analysis, anom=anom, obs=obs, model=model

; Reads SST data from UM expt and writes files for kriging process
; expt - experiment id
; analysis - title of analysis = sub-directory name for files
; /anom - calculate anomalies with expt tczyc
; /obs - write out as observations at PRISM locations
; /model - write out as model data

;Averaging period for standard error of mean when using model as "observations"
n_av=50

; Read SST data for expt
dir_um = 'O:\Documents\11_Projects\03_PhD_Pliocene_CO2\HadCM3_jobs\webpages\'
read_ncvar, dir_um+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','longitude',lon
read_ncvar, dir_um+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','latitude',lat
read_ncvar, dir_um+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','temp_mm_uo',temp
temp2=temp

; Anomaly data
if(keyword_set(anom)) then begin
   read_ncvar, dir_um+'tczyc'+path_sep()+'climate'+path_sep()+'tczyco.pfclann.nc','temp_mm_uo',temp_cont
   temp=temp-temp_cont
   analysis=analysis+'_anom'
endif

; Shift hemispheres
ids = where(lon ge 180)
lon(ids) = lon(ids)-360.
lon = shift(lon,144)
temp = shift(temp,144,0)
temp2 = shift(temp2,144,0)

dir = 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\'

; Write gridded model data
if (keyword_set(model)) then begin
   ;Write ocean only points to file
   ;It turns out it doesn't want the lat/lon values as the position in the array is all that matters.
   dir_mod = dir+'Model_data\'+analysis+path_sep()
   if(~ file_test(dir_mod)) then file_mkdir,dir_mod
   get_lun, unit
   openw, unit, dir_mod+expt+'.txt'
   get_lun, unit_mask
   openw, unit_mask, dir_mod+'mask.txt'
   
   for j = 0, n_elements(lat)-1 do begin
      for i = 0, n_elements(lon)-1 do begin
         n=n_elements(lon)*j+i+1
         if(temp2(i,j)lt 1000.) then begin
            printf, unit, temp(i,j)
            printf, unit_mask, n
         endif
      end
   end
   
   close, unit
   free_lun, unit
   close, unit_mask
   free_lun, unit_mask
   
   print, 'UM2Krig: ', n_elements(lon), min(lon), max(lon), n_elements(lat), min(lat), max(lat)
   
   set_plot, 'win'
   device, decomposed=0
   if (keyword_set(anom)) then begin
      get_levels, 'SST_del',  levels, ncol, fmt, ticknames
      one_two_tone_MO, 2, 'blue', ncol/2, 'red', ncol/2, 'white'
   endif else begin
      get_levels, 'SST',  levels, ncol, fmt, ticknames
      map_colour_table, 99, ncol, 0, 0, "white"
   endelse      
   bin_data, temp, temp_out, levels
   map_discontinuous, temp_out, lon, lat, ncol, 'Model SST!c'+analysis+' '+expt, 'lon', 'lat', ticknames, 0, 0, /legend, /map
   write_png, dir_mod+expt+'.png', tvrd(true=1)
   wdelete
   
endif

if (keyword_set(obs)) then begin
   ; Model data as observations
   ; Interpolate onto PRISM data locations
   dir_obs = dir+'Observation_data\'+analysis+path_sep()
   if(~ file_test(dir_obs)) then file_mkdir,dir_obs
   
   ; Uncertainty data for /obs case
   ; Webpages calculate sample standard deviation over specified averaging period
   ; Standard error of mean = s/sqrt(n)
   read_ncvar, dir_um+expt+path_sep()+'climate'+path_sep()+expt+'o.pfsdann.nc','temp_mm_uo',temp_sd
   temp_sd = shift(temp_sd,144,0)
   temp_se=temp_sd/sqrt(n_av)
   ; Probably ought to do something about combining uncertainties in anomaly cases!
   
   ; ; Read PRISM3 locations
   get_lun, unit
   openr, unit, dir+'Observation_data\PRISM3+_lon_lat.txt'
   data=read_ascii(dir+'Observation_data\PRISM3+_lon_lat.txt')
   close, unit
   free_lun, unit
   
   lon_prism =reform(data.field1(0,*))
   lat_prism =reform(data.field1(1,*))
   lon_prism_shift =reform(data.field1(2,*))
   lat_prism_shift =reform(data.field1(3,*))
   
   if keyword_set(anom) then begin
      colours = 'SST_del'
   endif else begin
      colours = 'SST'
      endelse      

   plot_setup, 'png', scale, charscale, linescale
   
   dig_points, lon_prism_shift, lat_prism_shift, lon, lat, temp, 2e+020, points, colours
   xyouts,0.5,0.93,'Model SST - PRISM3 Locations!c'+analysis+' '+expt,alignment=0.5,/normal, color=fsc_color("Black"), font=1, charsize=3*charscale
   write_png, dir_obs+expt+'.png', tvrd(true=1)
   wdelete
   
   dig_points, lon_prism_shift, lat_prism_shift, lon, lat, temp_se, 2e+020, points_se, 'SST_SD2'
   xyouts,0.5,0.93,'Model SST S.Err - PRISM3 Locations!c'+analysis+' '+expt,alignment=0.5,/normal, color=fsc_color("Black"), font=1, charsize=3*charscale
   write_png, dir_obs+expt+'_SE.png', tvrd(true=1)
   wdelete
   
   ;Write out to file
   get_lun, unit
   openw, unit, dir_obs+expt+'.txt'
      printf, unit, "x      y      z      std"
   for i=0,n_elements(points)-1 do begin
      ;              lon,         lat,          anomaly,   uncertainty
      printf, unit, lon_prism(i), lat_prism(i), points(i), points_se(i)
   endfor
   close, unit
   free_lun, unit

endif

return
end