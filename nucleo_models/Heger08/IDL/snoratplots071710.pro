tek_color
restore,'D:\Heger08\oxdbdata.sav'
ff=findfile('*@nucleo.idl',count=cnt)
psinit,0,file='snoratprofs.ps'
device,/color
tek_color
first=1
for i=0,cnt-1 do begin
pos2=strpos(ff(i),'A@')
ms=fix(strmid(ff(i),1,pos2-1))
if ms lt 50 then begin ;only consider low SN masses

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

plot_io,totdata.mass,totdata.h1,yrange=[1e-6,1e2],xrange=[2,.4*ms],title=ff(i)
oplot,totdata.mass,r18,color=2
oplot,totdata.mass,r17,color=3
oplot,totdata.mass,totdata.he4,color=4
oplot,totdata.mass,totdata.c12,color=5




f=maken(0,1,1e3)
;i1=5
;i2=6 ;He/N
wf=where(finite(r18))
i1=where(r18 eq max(r18(wf)))
ver,totdata(i1).mass
wf2=where(finite(r17))
i2=where(r17 eq max(r17(wf2)))
ver,totdata(i2).mass


endif
erase
endfor

psterm



end
