pro um2krig_all
expts = ['tdgth', 'tczyi', 'tdgtj', 'tczyj', 'tdgtg', 'tczyk', 'tdgtk', 'tdgti']

for i=0, n_elements(expts)-1 do begin
   um2krig, expts[i], 'CO2', /model, /obs
   um2krig, expts[i], 'CO2', /model, /obs, /anom
endfor

return
end
