pro plot_matrix, analysis, pngorps

; Setup plotting
if (pngorps eq 'png') then begin
   set_plot, 'win'
   device, decomposed=0
endif else begin
   set_plot, 'ps'
endelse
plot_setup, pngorps, scale, charscale, linescale
device, set_font='helvetica', /tt_font

; Read results matrix
get_lun, unit
openr, unit, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\'+analysis+'\Results_matrix.txt'
data=read_ascii('C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\Results_matrix.txt', data_Start=1)
close, unit
free_lun, unit

results = data.field1[1:8,0:7]
x=findgen(8)+1
co2    = [280, 315, 350, 375, 405, 475, 560, 1000]

; Set plotting levels
get_levels, 'krig', levels, ncol, format, ticknames
; Bin data for each dataset if gridded data required
bin_data, results, results_out, levels

ncol=n_elements(levels)
map_colour_table, 99, ncol, 0, 0, "white"


;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\Results_matrix'
title='Results Matrix'
xtitle='UM Model'
ytitle='Fit Field'
extra_leg=2
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

;imdisp, results_out, /noscale, color=fsc_color('Black'), position=[0.1,0.1,0.9,0.9], /normal

;map_discontinuous,  results_out, x, x, ncol, title, xtitle, ytitle, ticknames, extra_leg, 0,  /legend
plot_discontinuous, results_out, x, x, ncol, title, xtitle, ytitle, ticknames, extra_leg, 0, /legend
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse
;     END PLOT

;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\Results_matrix_co2'
title='Results Matrix'
xtitle='UM Model'
ytitle='Fit Field'
extra_leg=2
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

contour_fjb,    results, co2, co2, ncol, levels, title, xtitle, ytitle, ticknames, 0, 0, /legend
;     Save plot
if (pngorps eq 'png') then begin
   write_png, filename+'.png', tvrd(true=1)
   wdelete
endif else begin
   device, /close
endelse
;     END PLOT


;Create filename stem for plot
filename='C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\Output\Results_matrix_ln'
title='Results Matrix'
xtitle='UM Model'
ytitle='Fit Field'
extra_leg=2
if (pngorps eq 'ps') then begin
   device, file=filename+'.eps', /encapsulated, /color, bits_per_pixel=8
endif   

contour_fjb,    results, alog(co2), alog(co2), ncol, levels, title, xtitle, ytitle, ticknames, 0, 0, /legend
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
