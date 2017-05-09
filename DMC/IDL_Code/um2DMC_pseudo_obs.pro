pro um2DMC_pseudo_obs, expt, pngorps, leg_plot

; Sample gridded climate data (modern observations or model results) at data locations
; Create inputs for DMC process, both as models and as "pseudo-observations"
; Test several types of prior

plot_setup, pngorps, scale, charscale, linescale
leg=''
if (leg_plot) then leg='_leg'

; Read SST data for expt
if (expt eq 'HadISST') then begin
   read_hadisst, temp, temp_sd, lon, lat
   ; 2 pts in Canadian Arctic end up NaN SD and are effectively ignored
   ; set to 3 to check influence of these points
   ;temp_sd(where(temp_sd gt 100)) = 3.
   nlon=n_elements(lon)
   nlat=n_elements(lat)
endif else begin  
   dir = 'O:\Documents\11_Projects\03_PhD_Pliocene_CO2\HadCM3_jobs\webpages\'
   read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','longitude',lon
   read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','latitude',lat
   read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','temp_mm_uo',temp
   read_ncvar, dir+'tczyc\climate\tczyco.pfclann.nc','temp_mm_uo',temp_cont
   read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfsdann.nc','temp_mm_uo',temp_sd
   temp(where(temp gt 1000.)) = -999.
   temp_anom=temp-temp_cont
   temp_anom(where(temp gt 1000.)) = -999.
   nlon=n_elements(lon)
   nlat=n_elements(lat)
   ids = where(lon ge 180)
   lon(ids) = lon(ids)-360.
   lon = shift(lon,nlon/2)
   temp = shift(temp,nlon/2,0)
   temp_anom = shift(temp_anom,nlon/2,0)
   temp_cont = shift(temp_cont,nlon/2,0)
   temp_sd = shift(temp_sd,nlon/2,0)
endelse



;Write ocean only points to file
;It turns out it doesn't want the lat/lon values as the position in the array is all that matters.
model_dir='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Model_data\'+expt+path_sep()
file_mkdir, model_dir

get_lun, unit_cosine
openw, unit_cosine, model_dir+'cosine.txt'
get_lun, unit_cosine2
openw, unit_cosine2, model_dir+'cosine2.txt'
get_lun, unit_abs
openw, unit_abs, model_dir+'abs.txt'
get_lun, unit_latav
openw, unit_latav, model_dir+'latav.txt'
get_lun, unit_globav
openw, unit_globav, model_dir+'globav.txt'
if (expt ne 'HadISST') then begin
   get_lun, unit_anom
   openw, unit_anom, model_dir+'anom.txt'
endif
get_lun, unit_lsm
openw, unit_lsm, model_dir+'mask.txt'



mask=make_array(nlon,nlat,value=1,/integer)
mask(where(temp lt -500.)) = 0
latav = average_zonal(temp, lat, mask)
globav = average_global(temp, lon, lat, mask)
globav = replicate( globav, nlat)
; Estimate bounds of cosine curves from latav
hi=max(latav,/NaN)
lo=min(latav,/NaN)
cosine = (cos(lat*!pi/180)*(hi-lo))+lo
cosine2 = (cos(lat*!pi/180)*cos(lat*!pi/180)*(hi-lo))+lo
if (expt ne 'HadISST') then cont = average_zonal(temp_cont, lat, mask)


filename = model_dir+'latmean'+leg
plot_open, pngorps, filename
make_plot, 'Zonal Average Temperature: '+expt, 'Lat','Temperature (degrees C)',  [-90,90],45,[lo,hi],0, 16, 3.*leg_plot
oplot,lat,latav,color=fsc_color("Green"),thick=2
oplot,lat,globav,color=fsc_color("Red"),thick=2
oplot,lat,cosine,color=fsc_color("Royal Blue"),thick=2
oplot,lat,cosine2,color=fsc_color("Black"),thick=2
if (expt ne 'HadISST') then oplot,lat,cont,color=fsc_color("Sienna"),thick=2
if (leg_plot) then begin
   xyouts, 0.8,0.84,'Lat. Av.', color=fsc_color("Green"), charsize=2.5*charscale, charthick=1.5*linescale, font=1, /normal
   xyouts, 0.8,0.8 ,'Global Av.', color=fsc_color("Red"), charsize=2.5*charscale, charthick=1.5*linescale, font=1, /normal
   xyouts, 0.8,0.76,'Cosine', color=fsc_color("Royal Blue"), charsize=2.5*charscale, charthick=1.5*linescale, font=1, /normal
   xyouts, 0.8,0.72,'Cosine!u2!n', color=fsc_color("Black"), charsize=2.5*charscale, charthick=1.5*linescale, font=1, /normal
   if (expt ne 'HadISST') then xyouts, 0.8,0.68,'PI', color=fsc_color("Sienna"), charsize=2.5*charscale, charthick=1.5*linescale, font=1, /normal
   xyouts, 0.8,0.6,'Min: '+strcompress(string(lo,format='(f5.2)'),/remove_all), color=fsc_color("Black"), charsize=2.5*charscale, charthick=1.5*linescale, font=1, /normal
   xyouts, 0.8,0.56 ,'Max: '+strcompress(string(hi,format='(f5.2)'),/remove_all), color=fsc_color("Black"), charsize=2.5*charscale, charthick=1.5*linescale, font=1, /normal
endif   
plot_close, pngorps, filename


; This bit makes the latitudianlly averaged vectors into a global array for easy subtraction
cosine2=reform(cosine2,1,nlat)
cosine2=rebin(cosine2,nlon,nlat)
temp_cosine2=temp-cosine2
temp_cosine2(where(temp gt 100.)) = -999

cosine=reform(cosine,1,nlat)
cosine=rebin(cosine,nlon,nlat)
temp_cosine=temp-cosine
temp_cosine(where(temp gt 100.)) = -999

latav=reform(latav,1,nlat)
latav=rebin(latav,nlon,nlat)
temp_latav=temp-latav
temp_latav(where(temp gt 100.)) = -999

globav=reform(globav,1,nlat)
globav=rebin(globav,nlon,nlat)
temp_globav=temp-globav
temp_globav(where(temp gt 100.)) = -999

; Write out files
for j = 0, n_elements(lat)-1 do begin
   for i = 0, n_elements(lon)-1 do begin
      n=n_elements(lon)*j+i+1
      if(temp(i,j) gt -500.) then begin
         printf, unit_cosine2, temp_cosine2(i,j)
         printf, unit_lsm, n
      endif
   end
end

for j = 0, n_elements(lat)-1 do begin
   for i = 0, n_elements(lon)-1 do begin
      n=n_elements(lon)*j+i+1
      if(temp(i,j)gt -500.) then begin
         printf, unit_cosine, temp_cosine(i,j)
      endif
   end
end
for j = 0, n_elements(lat)-1 do begin
   for i = 0, n_elements(lon)-1 do begin
      n=n_elements(lon)*j+i+1
      if(temp(i,j)gt -500.) then begin
         printf, unit_abs, temp(i,j)
      endif
   end
end
for j = 0, n_elements(lat)-1 do begin
   for i = 0, n_elements(lon)-1 do begin
      n=n_elements(lon)*j+i+1
      if(temp(i,j)gt -500.) then begin
         printf, unit_latav, temp_latav(i,j)
      endif
   end
end
for j = 0, n_elements(lat)-1 do begin
   for i = 0, n_elements(lon)-1 do begin
      n=n_elements(lon)*j+i+1
      if(temp(i,j)gt -500.) then begin
         printf, unit_globav, temp_globav(i,j)
      endif
   end
end
if (expt ne 'HadISST') then begin
   for j = 0, n_elements(lat)-1 do begin
      for i = 0, n_elements(lon)-1 do begin
         n=n_elements(lon)*j+i+1
         if(temp(i,j)gt -500.) then begin
            printf, unit_anom, temp_anom(i,j)
         endif
      end
   end
endif
close, unit_cosine
free_lun, unit_cosine
close, unit_cosine2
free_lun, unit_cosine2
close, unit_abs
free_lun, unit_abs
close, unit_latav
free_lun, unit_latav
close, unit_globav
free_lun, unit_globav
if (expt ne 'HadISST') then begin
   close, unit_anom
   free_lun, unit_anom
endif   
close, unit_lsm
free_lun, unit_lsm

print, n_elements(lon), min(lon), max(lon), n_elements(lat), min(lat), max(lat)

; Plot model data
get_levels, 'SST',  levels, ncol, fmt, ticknames
bin_data, temp, temp_out, levels
map_colour_table, 99, ncol, 0, 0, "white"
filename = model_dir+'abs'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_out, lon, lat, ncol, 'Model SST: Abs!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


get_levels, 'SST_del',  levels, ncol, fmt, ticknames
bin_data, temp_cosine2, temp_cosine2_out, levels
bin_data, temp_cosine, temp_cosine_out, levels
bin_data, temp_latav, temp_latav_out, levels
bin_data, temp_anom, temp_anom_out, levels
ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
filename = model_dir+'cosine2'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_cosine2_out, lon, lat, ncol, 'Model SST Anomalies: cosine!u2!d!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse

filename = model_dir+'cosine'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_cosine_out, lon, lat, ncol, 'Model SST Anomalies: cosine!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse

filename = model_dir+'latav'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_latav_out, lon, lat, ncol, 'Model SST Anomalies: Lat av.!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse

if (expt ne 'HadISST') then begin
   filename = model_dir+'anom'+leg
   if (pngorps eq 'ps') then begin
      device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
   endif
   map_discontinuous, temp_anom_out, lon, lat, ncol, 'Model SST Anomalies: PI!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
   ;     Save plot
   if (pngorps eq 'png') then begin
      write_png, filename+'.png', tvrd(true=1)
      wdelete
   endif else begin
      device, /close
   endelse
endif   

get_levels, 'SST_del2',  levels, ncol, fmt, ticknames
bin_data, temp_globav, temp_globav_out, levels
ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
filename = model_dir+'globav'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_globav_out, lon, lat, ncol, 'Model SST Anomalies: globav!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Create pseudo-observations at PRISM data locations
obs_dir = 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Observation_data\'+expt+path_sep()
file_mkdir, obs_dir

; Read PRISM3 locations
get_lun, unit
openr, unit, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Observation_data\PRISM3+_lon_lat.txt'
data=read_ascii('C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Observation_data\PRISM3+_lon_lat.txt',data_Start=1)
close, unit
free_lun, unit

lon_prism =reform(data.field1(0,*))
lat_prism =reform(data.field1(1,*))
lon_prism_UM =reform(data.field1(2,*))
lat_prism_UM =reform(data.field1(3,*))
lon_prism_Had =reform(data.field1(4,*))
lat_prism_Had =reform(data.field1(5,*))

if (expt eq 'HadISST') then begin
   x=lon_prism_Had
   y=lat_prism_Had
endif else begin
   print, "UM shifted data points"
   x=lon_prism_UM
   y=lat_prism_UM
endelse

print, "Extracting pseudo-observations from data'
dig_points, x, y, lon, lat, temp_cosine2, -999., prism_points_cosine2, ''
dig_points, x, y, lon, lat, temp_cosine, -999., prism_points_cosine, ''
dig_points, x, y, lon, lat, temp, -999., prism_points_abs, ''
dig_points, x, y, lon, lat, temp_latav, -999., prism_points_latav, ''
dig_points, x, y, lon, lat, temp_globav, -999., prism_points_globav, ''
if (expt ne 'HadISST') then dig_points, x, y, lon, lat, temp_anom, -999., prism_points_anom, ''
dig_points, x, y, lon, lat, temp_sd, -999., prism_points_sd, ''


; Plot vs whole data field and just as points with no context
get_levels, 'SST',  levels, ncol, fmt, ticknames
map_colour_table, 99, ncol, 0, 0, "white"
filename = obs_dir+'abs_prism'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_out, lon, lat, ncol, 'Model SST: Abs!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
bin_data, prism_points_abs, prism_points_abs_out, levels
not_nans=where(prism_points_abs gt -500.)
nans    =where(prism_points_abs lt -500.)
for i=0,n_elements(not_nans)-1 do begin
   xx=lon_prism(not_nans(i)) & yy=lat_prism(not_nans(i))
   polyfill,circle(xx,yy,2,1),/data,color=prism_points_abs_out(not_nans(i))
   plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor   
if (nans(0) ne -1) then begin
   for i=0,n_elements(nans)-1 do begin
      xx=lon_prism(nans(i)) & yy=lat_prism(nans(i))
      plots,square(xx,yy,2,1),/data,color=fsc_color("Royal Blue"), thick=2
      print, "NaN point ",xx, yy
   endfor   
endif
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Plot points only
filename = obs_dir+'abs_prism_points'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_points,    prism_points_abs_out, lon_prism, lat_prism, ncol, 'Model SST: Abs!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


get_levels, 'SST_del',  levels, ncol, fmt, ticknames
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
filename = obs_dir+'cosine2_prism'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_cosine2_out, lon, lat, ncol, 'Model SST Anomalies: cosine!u2!d!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
bin_data, prism_points_cosine2, prism_points_cosine2_out, levels
not_nans=where(prism_points_cosine2 gt -500.)
nans    =where(prism_points_cosine2 lt -500.)
for i=0,n_elements(not_nans)-1 do begin
   xx=lon_prism(not_nans(i)) & yy=lat_prism(not_nans(i))
   polyfill,circle(xx,yy,2,1),/data,color=prism_points_cosine2_out(not_nans(i))
   plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor  
if (nans(0) ne -1) then begin
   for i=0,n_elements(nans)-1 do begin
      xx=lon_prism(nans(i)) & yy=lat_prism(nans(i))
      plots,square(xx,yy,2,1),/data,color=fsc_color("Royal Blue"), thick=2
   endfor  
endif 
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Plot points only
filename = obs_dir+'cosine2_prism_points'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_points,    prism_points_cosine2_out, lon_prism, lat_prism, ncol, 'Model SST: cosine!u2!d!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


filename = obs_dir+'cosine_prism'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_cosine_out, lon, lat, ncol, 'Model SST Anomalies: cosine!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
bin_data, prism_points_cosine, prism_points_cosine_out, levels
not_nans=where(prism_points_cosine gt -500.)
nans    =where(prism_points_cosine lt -500.)
for i=0,n_elements(not_nans)-1 do begin
   xx=lon_prism(not_nans(i)) & yy=lat_prism(not_nans(i))
   polyfill,circle(xx,yy,2,1),/data,color=prism_points_cosine_out(not_nans(i))
   plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor   
if (nans(0) ne -1) then begin
   for i=0,n_elements(nans)-1 do begin
      xx=lon_prism(nans(i)) & yy=lat_prism(nans(i))
      plots,square(xx,yy,2,1),/data,color=fsc_color("Royal Blue"), thick=2
   endfor   
endif
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Plot points only
filename = obs_dir+'cosine_prism_points'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_points,    prism_points_cosine_out, lon_prism, lat_prism, ncol, 'Model SST: cosine!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


filename = obs_dir+'latav_prism'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_latav_out, lon, lat, ncol, 'Model SST Anomalies: Lat av.!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
bin_data, prism_points_latav, prism_points_latav_out, levels
not_nans=where(prism_points_latav gt -500.)
nans    =where(prism_points_latav lt -500.)
for i=0,n_elements(not_nans)-1 do begin
   xx=lon_prism(not_nans(i)) & yy=lat_prism(not_nans(i))
   polyfill,circle(xx,yy,2,1),/data,color=prism_points_latav_out(not_nans(i))
   plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor   
if (nans(0) ne -1) then begin
   for i=0,n_elements(nans)-1 do begin
      xx=lon_prism(nans(i)) & yy=lat_prism(nans(i))
      plots,square(xx,yy,2,1),/data,color=fsc_color("Royal Blue"), thick=2
   endfor   
endif
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Plot points only
filename = obs_dir+'latav_prism_points'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_points,    prism_points_latav_out, lon_prism, lat_prism, ncol, 'Model SST: Lat av.!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse

if (expt ne 'HadISST') then begin   
   filename = obs_dir+'anom_prism'+leg
   if (pngorps eq 'ps') then begin
      device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
   endif
   map_discontinuous, temp_anom_out, lon, lat, ncol, 'Model SST Anomalies: PI!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
   bin_data, prism_points_anom, prism_points_anom_out, levels
   not_nans=where(prism_points_anom gt -500.)
   nans    =where(prism_points_anom lt -500.)
   for i=0,n_elements(not_nans)-1 do begin
      xx=lon_prism(not_nans(i)) & yy=lat_prism(not_nans(i))
      polyfill,circle(xx,yy,2,1),/data,color=prism_points_anom_out(not_nans(i))
      plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
   endfor   
   if (nans(0) ne -1) then begin
      for i=0,n_elements(nans)-1 do begin
         xx=lon_prism(nans(i)) & yy=lat_prism(nans(i))
         plots,square(xx,yy,2,1),/data,color=fsc_color("Royal Blue"), thick=2
      endfor   
   endif
   ;     Save plot
   if (pngorps eq 'png') then begin
      write_png, filename+'.png', tvrd(true=1)
      wdelete
   endif else begin
      device, /close
   endelse
   
   
   ; Plot points only
   filename = obs_dir+'anom_prism_points'+leg
   if (pngorps eq 'ps') then begin
      device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
   endif
   map_points,    prism_points_anom_out, lon_prism, lat_prism, ncol, 'Model SST: PI!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
   ;     Save plot
   if (pngorps eq 'png') then begin
      write_png, filename+'.png', tvrd(true=1)
      wdelete
   endif else begin
      device, /close
   endelse
endif

get_levels, 'SST_del2',  levels, ncol, fmt, ticknames
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
filename = obs_dir+'globav_prism'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_globav_out, lon, lat, ncol, 'Model SST Anomalies: globav!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
bin_data, prism_points_globav, prism_points_globav_out, levels
not_nans=where(prism_points_globav gt -500.)
nans    =where(prism_points_globav lt -500.)
for i=0,n_elements(not_nans)-1 do begin
   xx=lon_prism(not_nans(i)) & yy=lat_prism(not_nans(i))
   polyfill,circle(xx,yy,2,1),/data,color=prism_points_globav_out(not_nans(i))
   plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
endfor   
if (nans(0) ne -1) then begin
   for i=0,n_elements(nans)-1 do begin
      xx=lon_prism(nans(i)) & yy=lat_prism(nans(i))
      plots,square(xx,yy,2,1),/data,color=fsc_color("Royal Blue"), thick=2
   endfor   
endif
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Plot points only
filename = obs_dir+'globav_prism_points'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_points,    prism_points_globav_out, lon_prism, lat_prism, ncol, 'Model SST: globav!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Some points have 0 SD, especially where there is permanent sea ice, 
; but this will break the algorithm when it calculates 1/sigma2, so set to small value.

limit=0.1
ids = where(prism_points_sd lt limit)
if (ids(0) ne -1) then begin
   print, "SDs reset to ",limit
   for i=0,n_elements(ids)-1 do begin
      print, lon(ids(i)),lat(ids(i)),prism_points_sd(ids(i))
   endfor
   prism_points_sd(ids)=limit
endif

get_levels, 'SST_SD3',  levels, ncol, fmt, ticknames
temp_sd(where(temp_sd gt 1000.)) = -999
bin_data, temp_sd, temp_sd_out, levels
bin_data, prism_points_sd, prism_points_sd_out, levels
map_colour_table, 99, ncol, 0, 0, "white"
filename = obs_dir+'sd_prism'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_discontinuous, temp_sd_out, lon, lat, ncol, 'Model SD!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot
not_nans=where(prism_points_sd gt -500.)
nans    =where(prism_points_sd lt -500.)
for i=0,n_elements(not_nans)-1 do begin
   x=lon_prism(not_nans(i)) & y=lat_prism(not_nans(i))
   polyfill,circle(x,y,2,1),/data,color=prism_points_sd_out(not_nans(i))
   plots,   circle(x,y,2,1),/data,color=fsc_color("Black")
endfor 
if (ids(0) ne -1) then begin
   for i=0,n_elements(ids)-1 do begin
      x=lon_prism(ids(i)) & y=lat_prism(ids(i))
      plots,circle(x,y,2,1),/data,color=fsc_color("Yellow"), thick=2
   endfor  
endif  
if (nans(0) ne -1) then begin
   for i=0,n_elements(nans)-1 do begin
      x=lon_prism(nans(i)) & y=lat_prism(nans(i))
      plots,square(x,y,2,1),/data,color=fsc_color("Royal Blue"), thick=2
   endfor  
endif 
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


; Plot points only
filename = obs_dir+'sd_prism_points'+leg
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps',/encapsulated, /color, bits_per_pixel=8
endif
map_points,    prism_points_sd_out, lon_prism, lat_prism, ncol, 'Model SST: SD!c'+expt, 'lon', 'lat', ticknames, 0, 0, legend=leg_plot, /map
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse


;Write out to file

get_lun, unit
openw, unit, obs_dir+'prism_cosine2.txt'
   printf, unit, "x      y      z      std"
for i=0,n_elements(prism_points_cosine2)-1 do begin
   ;              lon,         lat,          anomaly,         uncertainty
   printf, unit, lon_prism(i), lat_prism(i), prism_points_cosine2(i), prism_points_sd(i)
endfor
close, unit
free_lun, unit

get_lun, unit
openw, unit, obs_dir+'prism_cosine.txt'
   printf, unit, "x      y      z      std"
for i=0,n_elements(prism_points_cosine)-1 do begin
   ;              lon,         lat,          anomaly,         uncertainty
   printf, unit, lon_prism(i), lat_prism(i), prism_points_cosine(i), prism_points_sd(i)
endfor
close, unit
free_lun, unit

get_lun, unit
openw, unit, obs_dir+'prism_latav.txt'
   printf, unit, "x      y      z      std"
for i=0,n_elements(prism_points_latav)-1 do begin
   ;              lon,         lat,          anomaly,         uncertainty
   printf, unit, lon_prism(i), lat_prism(i), prism_points_latav(i), prism_points_sd(i)
endfor
close, unit
free_lun, unit

get_lun, unit
openw, unit, obs_dir+'prism_globav.txt'
   printf, unit, "x      y      z      std"
for i=0,n_elements(prism_points_globav)-1 do begin
   ;              lon,         lat,          anomaly,         uncertainty
   printf, unit, lon_prism(i), lat_prism(i), prism_points_globav(i), prism_points_sd(i)
endfor
close, unit
free_lun, unit

get_lun, unit
openw, unit, obs_dir+'prism_abs.txt'
   printf, unit, "x      y      z      std"
for i=0,n_elements(prism_points_abs)-1 do begin
   ;              lon,         lat,          anomaly,         uncertainty
   printf, unit, lon_prism(i), lat_prism(i), prism_points_abs(i), prism_points_sd(i)
endfor
close, unit
free_lun, unit

if (expt ne 'HadISST') then begin
   get_lun, unit
   openw, unit, obs_dir+'prism_anom.txt'
      printf, unit, "x      y      z      std"
   for i=0,n_elements(prism_points_anom)-1 do begin
      ;              lon,         lat,          anomaly,         uncertainty
      printf, unit, lon_prism(i), lat_prism(i), prism_points_anom(i), prism_points_sd(i)
   endfor
   close, unit
   free_lun, unit
endif

return
end