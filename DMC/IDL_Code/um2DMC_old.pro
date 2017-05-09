pro um2krig, expt

; Read SST data for expt
dir = 'O:\Documents\11_Projects\03_PhD_Pliocene_CO2\HadCM3_jobs\webpages\'
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','longitude',lon
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','latitude',lat
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfclann.nc','temp_mm_uo',temp
read_ncvar, dir+expt+path_sep()+'climate'+path_sep()+expt+'o.pfsdann.nc','temp_mm_uo',temp_sd
read_ncvar, dir+'tczyc'+path_sep()+'climate'+path_sep()+'tczyco.pfclann.nc','temp_mm_uo',temp_cont

;Write ocean only points to file
;It turns out it doesn't want the lat/lon values as the position in the array is all that matters.
get_lun, unit
openw, unit, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'.txt'
get_lun, unit_anom
openw, unit_anom, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'_anom.txt'
get_lun, unit_lsm
openw, unit_lsm, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\lsm.txt'

ids = where(lon ge 180)
lon(ids) = lon(ids)-360.
lon = shift(lon,144)
temp = shift(temp,144,0)
temp_cont = shift(temp_cont,144,0)
temp_sd = shift(temp_sd,144,0)


for j = 0, n_elements(lat)-1 do begin
   for i = 0, n_elements(lon)-1 do begin
      n=n_elements(lon)*j+i+1
      if(temp(i,j)lt 1000.) then begin
         printf, unit, temp(i,j)
         printf, unit_anom, temp(i,j)-temp_cont(i,j)
         printf, unit_lsm, n
      endif
   end
end

close, unit
free_lun, unit
close, unit_anom
free_lun, unit_anom
close, unit_lsm
free_lun, unit_lsm

print, n_elements(lon), min(lon), max(lon), n_elements(lat), min(lat), max(lat)

set_plot, 'win'
device, decomposed=0
get_levels, 'SST',  levels, ncol, fmt, ticknames
bin_data, temp, temp_out, levels
map_colour_table, 99, ncol, 0, 0, "white"
map_discontinuous, temp_out, lon, lat, ncol, 'Model Data', 'lon', 'lat', ticknames, 0, 0, /legend, /map
write_png, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'.png', tvrd(true=1)
wdelete

get_levels, 'SST_del',  levels, ncol, fmt, ticknames
bin_data, temp-temp_cont, temp_anom_out, levels
ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
map_discontinuous, temp_anom_out, lon, lat, ncol, 'Model Data Anomalies', 'lon', 'lat', ticknames, 0, 0, /legend, /map
write_png, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'_anom.png', tvrd(true=1)
wdelete

; Interpolate onto PRISM data locations
; ; Read PRISM3 locations
get_lun, unit
openr, unit, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\observation_data\PRISM3_lon_lat.txt'

data=read_ascii('C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\observation_data\PRISM3_lon_lat.txt')

close, unit
free_lun, unit

lon_prism =reform(data.field1(0,*))
lat_prism =reform(data.field1(1,*))
lon_prism_shift =reform(data.field1(2,*))
lat_prism_shift =reform(data.field1(3,*))

x = VALUE_LOCATE(lon-0.625,lon_prism_shift) 
y = VALUE_LOCATE(lat-0.625,lat_prism_shift)
prism_points=temp[x,y]-temp_cont(x,y)
prism_points_sd=temp_sd[x,y]
;print, "SST interpolated at PRISM3 data locations "
;print, temp[x,y]
;print, "SST cont interpolated at PRISM3 data locations "
;print, temp_cont(x,y)
;print, "Anomalies interpolated at PRISM3 data locations "
;print, prism_points
;print, "SDs interpolated at PRISM3 data locations "
;print, prism_points_sd
print, "NaNs"
nans=where(temp[x,y] gt 100.)
not_nans=where(temp[x,y] le 100.)
if (nans(0) ne -1) then begin
    print, "NaN indices"
    print, x(nans)
    print, y(nans)
    for i=0,n_elements(nans)-1 do begin
      print, strcompress(string(i)), lon(x(nans(i))),lat(y(nans(i))), lon_prism_shift(nans(i)),lat_prism_shift(nans(i))
      print, temp(x(nans(i))-1:x(nans(i))+1,y(nans(i))-1:y(nans(i))+1)
    endfor   
endif

; Check anomaly locations and values, find NaNs
map_discontinuous, temp_anom_out, lon, lat, ncol, 'Model Anomalies', 'lon', 'lat', ticknames, 0, 0, /legend
bin_data, prism_points, prism_points_out, levels
for i=0,n_elements(not_nans)-1 do begin
   x=lon_prism(not_nans(i)) & y=lat_prism(not_nans(i))
   polyfill,circle(x,y,2,1),/data,color=prism_points_out(not_nans(i))
   plots,   circle(x,y,2,1),/data,color=fsc_color("Black")
endfor   
for i=0,n_elements(nans)-1 do begin
   x=lon_prism(nans(i)) & y=lat_prism(nans(i))
   plots,square(x,y,2,1),/data,color=fsc_color("Royal Blue"), thick=2
endfor   
write_png, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\observation_data\'+expt+'_anom_prism.png', tvrd(true=1)
wdelete

get_levels, 'SST_SD',  levels, ncol, fmt, ticknames
temp_sd(where(temp_sd gt 1000.)) = -999
bin_data, temp_sd, temp_sd_out, levels
bin_data, prism_points_sd, prism_points_sd_out, levels
map_colour_table, 99, ncol, 0, 0, "white"
map_discontinuous, temp_sd_out, lon, lat, ncol, 'Model SD', 'lon', 'lat', ticknames, 0, 0, /legend
for i=0,n_elements(not_nans)-1 do begin
   x=lon_prism(not_nans(i)) & y=lat_prism(not_nans(i))
   polyfill,circle(x,y,2,1),/data,color=prism_points_sd_out(not_nans(i))
   plots,   circle(x,y,2,1),/data,color=fsc_color("Black")
endfor   
for i=0,n_elements(nans)-1 do begin
   x=lon_prism(nans(i)) & y=lat_prism(nans(i))
   plots,square(x,y,2,1),/data,color=fsc_color("Royal Blue"), thick=2
endfor   
write_png, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\observation_data\'+expt+'_sd_prism.png', tvrd(true=1)
wdelete

;Write out to file

get_lun, unit
openw, unit, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\observation_data\'+expt+'_prism.txt'
for i=0,n_elements(prism_points)-1 do begin
   ;              lon,         lat,          anomaly,         uncertainty
   printf, unit, lon_prism(i), lat_prism(i), prism_points(i), prism_points_sd(i)
endfor
close, unit
free_lun, unit

;---------------------------------------------------------------------------------

;Regrid onto atmosphere resolution
; Shift lon and temp a further 1 cell around to align with whole atmosphere cell
print, "lon array before shifting: ", lon[0:5]
lon = shift(lon,1)
lon[0] = lon[0]-360.
print, "Shifted lon array before rebinning: ", lon[0:5]
lon_at = rebin(lon,96)
print, "Shifted lon array after rebinning: ", lon_at[0:5]
; Shift temp arrays in same manner longitudinally
temp = shift(temp,1,0)
temp_cont = shift(temp_cont,1,0)

; Remove first and last rows from lat direction as these are not included in atmosphere cells
print, "lat array before rebinning: ", lat[0:5]
lat_at = rebin(lat[1:142],71)
print, "lat array after rebinning: ", lat_at[0:5]
; Rebin temp arrays in same manner
temp_at = rebin(temp[*,1:142],96,71)
temp_cont_at = rebin(temp_cont[*,1:142],96,71)

get_levels, 'SST',  levels, ncol, fmt, ticknames
bin_data, temp_at, temp_at_out, levels
map_colour_table, 99, ncol, 0, 0, "white"
map_discontinuous, temp_at_out, lon_at, lat_at, ncol, 'Reduced Resolution', 'lon', 'lat', ticknames, 0, 0, /legend, /map
write_png, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'_atm.png', tvrd(true=1)
wdelete

get_levels, 'SST_del',  levels, ncol, fmt, ticknames
bin_data, temp_at-temp_cont_at, temp_anom_at_out, levels
ncol=n_elements(levels)
ncol_side=ncol/2
one_two_tone_MO, 2, 'blue', ncol_side, 'red', ncol_side, 'white'
map_discontinuous, temp_anom_at_out, lon_at, lat_at, ncol, 'Reduced Resolution Anomalies', 'lon', 'lat', ticknames, 0, 0, /legend, /map
write_png, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'_anom_atm.png', tvrd(true=1)
wdelete

;Write ocean only points to file
;It turns out it doesn't want the lat/lon values as the position in the array is all that matters.
get_lun, unit
openw, unit, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'_atm.txt'
get_lun, unit_anom
openw, unit_anom, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\'+expt+'_anom_atm.txt'
get_lun, unit_lsm
openw, unit_lsm, 'C:\Users\fb7897\Dropbox\Work\Paleo_obs_tests\Fran\model_data\lsm_atm.txt'

for j = 0, n_elements(lat_at)-1 do begin
   for i = 0, n_elements(lon_at)-1 do begin
      n=n_elements(lon_at)*j+i+1
      if(temp_at(i,j)lt 1000.) then begin
         printf, unit, temp_at(i,j)
         printf, unit_anom, temp_at(i,j)-temp_cont_at(i,j)
         printf, unit_lsm, n
      endif
   end
end

close, unit
free_lun, unit
close, unit_anom
free_lun, unit_anom
close, unit_lsm
free_lun, unit_lsm

print, n_elements(lon_at), min(lon_at), max(lon_at), n_elements(lat_at), min(lat_at), max(lat_at)

return
end