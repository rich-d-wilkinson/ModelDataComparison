pro plot_krig, sst_filename, pngorps

sst_file=file_basename(sst_filename,'.txt')

;Read PRISM data points
read_PRISM3_SST_points, lat_prism, lon_prism, lon_shift, temp_plio, conf_plio, hadisst_modern

; Initialise arrays

lon_fit  = indgen(96)*3.75-180
lat_fit  = indgen(72)*2.5-90
std_fit  = make_array(96,72,/float, value=-999.0)
sst_fit  = make_array(96,72,/float, value=-999.0)

; Read my_fit from text file
get_lun, unit
openr, unit, "C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\output\"+sst_filename

dummy=' '
readf, unit, dummy, format='(A0)'

string=' '
while ~eof(unit) do begin
   readf,unit,string, format='(a0)'
   data = strsplit(string,'[]= ',/extract)
   i=fix(data(1))/15.+48
   j=fix(data(2))/10.+36
   ;lon_fit(i)=data(1)/4.+1.875
   ;lat_fit(j)=data(2)/4.+1.25
   std=float(data(4))
   sst=float(data(5))
   std_fit(i,j)=std
   sst_fit(i,j)=sst   
;   print, "Stuff"
endwhile

close, unit
free_lun, unit

; Plot fit output

; Setup plotting
if (pngorps eq 'png') then begin
   set_plot, 'win'
   device, decomposed=0
endif else begin
   set_plot, 'ps'
endelse
plot_setup, pngorps, scale, charscale, linescale
device, set_font='helvetica', /tt_font

; Set plotting levels
get_levels, 'SST_del', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, sst_fit, sst_fit_out, levels
bin_data, temp_plio-hadisst_modern, temp_plio_out, levels

;map_colour_table, 99, ncol, 0, 0, "white"
ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'

; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\output\'+sst_file
title='SST - Kriged Field!c'+sst_filename
xtitle='lon'
ytitle='lat'
extra_leg=3.5
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, sst_fit_out, lon_fit, lat_fit, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map
for i=0,n_elements(lon_prism)-1 do begin
   xx=lon_prism(i) & yy=lat_prism(i)
   case conf_plio(i) of
   4: begin
         polyfill,circle(xx,yy,3,1),/data,color=temp_plio_out(i)
         plots,   circle(xx,yy,3,1),/data,color=fsc_color("Black")
      end
   3: begin
         polyfill,square(xx,yy,3,1),/data,color=temp_plio_out(i)
         plots,   square(xx,yy,3,1),/data,color=fsc_color("Black")
      end
   2: begin
         polyfill,triangle(xx,yy,3,1),/data,color=temp_plio_out(i)
         plots,   triangle(xx,yy,3,1),/data,color=fsc_color("Black")
      end
   1: begin
         polyfill,diamond(xx,yy,3,1),/data,color=temp_plio_out(i)
         plots,   diamond(xx,yy,3,1),/data,color=fsc_color("Black")
      end
   endcase
endfor
         xyouts, 0.87, 0.78,'CONFIDENCE',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   circle(0.88,0.72,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.71,'Very high',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   square(0.88,0.66,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.65,'High',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   triangle(0.88,0.60,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.59,'Medium',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   diamond(0.88,0.54,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.53,'Low',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal

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
bin_data, std_fit, std_fit_out, levels
map_colour_table, 99, ncol, 0, 0, "white"

; Setup plotting
if (pngorps eq 'png') then begin
   set_plot, 'win'
   device, decomposed=0
endif else begin
   set_plot, 'ps'
endelse
plot_setup, pngorps, scale, charscale, linescale
device, set_font='helvetica', /tt_font

; FOR EACH PLOT REQUIRED
; CALL SPEICALIST PLOTTING SUBROUTINE
; Set plot attributes
;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\output\'+sst_file+'_SD'
title='SST - SD From Kriging!c'+sst_filename
xtitle='lon'
ytitle='lat'
extra_leg=3.5
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

map_discontinuous, std_fit_out, lon_fit, lat_fit, ncol,         title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend, /map

for i=0,n_elements(lon_prism)-1 do begin
   xx=lon_prism(i) & yy=lat_prism(i)
   case conf_plio(i) of
   4: begin
         plots,   circle(xx,yy,2,1),/data,color=fsc_color("Black")
      end
   3: begin
         plots,   square(xx,yy,2,1),/data,color=fsc_color("Black")
      end
   2: begin
         plots,   triangle(xx,yy,2,1),/data,color=fsc_color("Black")
      end
   1: begin
         plots,   diamond(xx,yy,2,1),/data,color=fsc_color("Black")
      end
   endcase
endfor
         xyouts, 0.87, 0.78,'CONFIDENCE',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   circle(0.88,0.72,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.71,'Very high',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   square(0.88,0.66,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.65,'High',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   triangle(0.88,0.60,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.59,'Medium',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal
         plots,   diamond(0.88,0.54,0.005,2),/normal,color=fsc_color("Black"),thick=2.*linescale
         xyouts, 0.9, 0.53,'Low',color=fsc_color("Black"), charsize=charscale*2.,charthick=2.*linescale, /normal

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