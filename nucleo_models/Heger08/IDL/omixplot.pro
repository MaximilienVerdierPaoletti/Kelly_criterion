ff=findfile('*@nucleo.idl',count=cnt)
psinit,0,file='snomix0610.ps'

for i=0,cnt-1 do begin
  restore,ff(i)
o16=zonedata.o16
o17=zonedata.o17
o18=zonedata.o18
o16=totdata.o16
o17=totdata.o17
o18=totdata.o18

;r18=zonedata.o18/zonedata.o16*16./18
;r17=zonedata.o17/zonedata.o16*16./17
r18=totdata.o18/totdata.o16*16./18
r17=totdata.o17/totdata.o16*16./17
wf=where(finite(r18))


f=maken(0,1,1e3)
;i1=5
;i2=6 ;He/N
wf=where(finite(r18))
i1=where(r18 eq max(r18(wf)))
wf=where(finite(r17))
i2=where(r17 eq max(r17(wf)))
i1=i1(0)
i2=i2(0)
m16=f*o16(i1)+(1.-f)*o16(i2)
m17=f*o17(i1)+(1.-f)*o17(i2)
m18=f*o18(i1)+(1.-f)*o18(i2)


mr18=m18/m16*16./18
mr17=m17/m16*16./17

plot_oo,r18,r17,psym=4,xrange=[1e-12,100],yrange=[1e-7,1],title=ff(i)
solarax,'o2'
oplot,mr18,mr17
plots,.08,2e-3,psym=sym('star')

oplot,totdata.o18/totdata.o16*16/18.,totdata.o17/totdata.o16*16/17.,psym=-5


erase
endfor
psterm



end
