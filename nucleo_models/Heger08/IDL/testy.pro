;restore h08 25m nucleo.idl
;409 is index of lowest mass zone for He/C
;903 is wind (so heC-env=409-902 inclusive)

restore,'C:\Users\Larry\SNmodels\Heger08\s25A@nucleo-mar2012.idl'
totdata=data

numd=n_elements(totdata)
zonedata=totdata(0)
tn=tag_names(totdata)
dm=totdata.mass-shift(totdata.mass,1)
w=lindgen(903-409)+409
i=0  
    tmass=total(dm(w))
    print,zones(i),min(dm(w)),max(dm(w)),n_elements(w)
    
    zonedata(i).zone=i+1
    zonedata(i).zonename='He/C through env'
    zonedata(i).mass=tmass
    
    for j=2,n_elements(tn)-3 do begin
    
      dum=execute('var=totdata(w).'+tn(j))
        tot=total(var*dm(w))
      mfrac=tot/tmass
     ;    dum=execute('zonedata(i).'+tn(j)+'=tot')
      dum=execute('zonedata(i).'+tn(j)+'=mfrac')
    endfor

end
