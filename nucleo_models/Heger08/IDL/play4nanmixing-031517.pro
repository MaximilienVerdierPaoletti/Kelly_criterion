si28=d.si28
si29=d.si29
si30=d.si30
ti48=d.ti48+d.cr48+d.v48
ti49=d.ti49+d.v49+d.cr49
ti50=d.ti50


num=1e3
f=maken(0,1,num)
ind1=1
ind2=5
ms1=d[ind1].mass
ms2=d[ind2].mass
m28=f*si28[ind1]/ms1+(1-f)*si28[ind2]/ms2
m29=f*si29[ind1]/ms1+(1-f)*si29[ind2]/ms2
m30=f*si30[ind1]/ms1+(1-f)*si30[ind2]/ms2
m48=f*ti48[ind1]/ms1+(1-f)*ti48[ind2]/ms2
m49=f*ti49[ind1]/ms1+(1-f)*ti49[ind2]/ms2
m50=f*ti50[ind1]/ms1+(1-f)*ti50[ind2]/ms2
m50b=f*ti50[ind1]/ms1+(1-f)*ti49[ind2]/ms2*stdrat(50)/stdrat(49)
md29=rat2del(m29/m28*28./29,stdrat(29))
md30=rat2del(m30/m28*28./30,stdrat(30))
md49=rat2del(m49/m48*48./49,stdrat(49))
md50=rat2del(m50/m48*48./50,stdrat(50))
md50b=rat2del(m50b/m48*48./50,stdrat(50))



a1=ti49[1]/ms1
a2=ti49[5]/ms2
b1=ti48[1]/ms1
b2=ti48[5]/ms2
r0=stdrat(49)

f1=f
f2=1-f
del=1e3*((f1*a1+f2*a2)/(f1*b1+f2*b2)/r0*48/49.-1)  ;d49
d1=1e3*((f1*a1)/(f1*b1+f2*b2)/r0*48/49.-1)
d2=1e3*((f2*a2)/(f1*b1+f2*b2)/r0*48/49.-1)


d49star=md49-md50-1000
d49starb=md49-md50b-1000

end
