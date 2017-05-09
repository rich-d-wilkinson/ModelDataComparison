pro read_DMC_output, out_dir, expt, lon, lat, temp, SD, mask

; Read gridded output of fitted field
; ----------------------------------
; FE 2 grid conversion runs from min to max lon and lat found in output and splits evenly into 'n=detail' points
; Find min and max and n_elements for lon/lat arrays, then can convert to indices to create mean/std arrays

; Read mean and std from text file
out_filename = out_dir+expt+'_grid_output.txt'
get_lun, unit
openr, unit, out_filename
data=read_ascii(out_filename, missing_value=-999, data_start=1)
close, unit
free_lun, unit

lon_fit  = reform(data.field01(1,*)*180./!pi)
lat_fit  = reform(data.field01(2,*)*180./!pi)
temp_fit = reform(data.field01(3,*))
SD_fit   = reform(data.field01(6,*))
mask_fit = reform(data.field01(9,*))

; Find unique values of lon and lat
lon = lon_fit(UNIQ(lon_fit, SORT(lon_fit)))
lat = lat_fit(UNIQ(lat_fit, SORT(lat_fit)))
; Convert lon_fit and lat_fit to indices for assigning temp and SD to arrays
lon_ind= round((lon_fit-min(lon_fit))/(lon[1]-lon[0]))
lat_ind= round((lat_fit-min(lat_fit))/(lat[1]-lat[0]))

; Generate temp array and populate
temp=make_array(n_elements(lon),n_elements(lat),/float,value=-999)
temp(lon_ind,lat_ind)=temp_fit
SD=make_array(n_elements(lon),n_elements(lat),/float,value=-999)
SD(lon_ind,lat_ind)=SD_fit
mask=make_array(n_elements(lon),n_elements(lat),/float,value=-999)
mask(lon_ind,lat_ind)=mask_fit
return
end
