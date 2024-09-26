restore,'D:\Heger08\IDL\oxdbdata.sav'
ff=findfile('*@nucleo.idl',count=cnt)
psinit,0,file='snomix070710.ps'
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



f=maken(0,1,1e3)
;i1=5
;i2=6 ;He/N
wf=where(finite(r18))
i1=where(r18 eq max(r18(wf)))
wf2=where(finite(r17))
i2=where(r17 eq max(r17(wf2)))
;print,ms,max(r17(wf2)),max(r18(wf))

i1=i1(0)
i2=i2(0)
m16=f*o16(i1)+(1.-f)*o16(i2)
m17=f*o17(i1)+(1.-f)*o17(i2)
m18=f*o18(i1)+(1.-f)*o18(i2)
m12=f*totdata(i1).c12+(1.-f)*totdata(i2).c12
m13=f*totdata(i1).c13+(1.-f)*totdata(i2).c13

mr18=m18/m16*16./18
mr17=m17/m16*16./17
mc_o=(m12+m13)/(m16+m18)*16/12.


if first then begin
 plot_oo,r18,r17,psym=4,xrange=[1e-6,1000],/xstyl,yrange=[1e-6,.1],/nod,$
 xtitle=rlab(18,16),ytit=rlab(17,16),charsi=2
solarax,'o2'
plots,.1,3e-3,psym=sym('star')
first=0
endif


oplot,mr18,mr17,color=i+1
xyouts,mr18(999),mr17(999),stc(ms),color=i+1
wo=where(mc_o lt 1)
if wo(0) ne -1 then oplot,mr18(wo),mr17(wo),color=i+1,thick=3

wc=where(mc_o ge 1)

;print,ms,mr17(0),mr18(0),mr17(999),mr18(999),mr18(wc(0))
 
;oplot,totdata.o18/totdata.o16*16/18.,totdata.o17/totdata.o16*16/17.,psym=-5

endif
endfor
oplot,dw.o18_16,dw.o17_16,psym=4,color=7,thick=1.8
oplot,dw2.o18_16,dw2.o17_16,psym=sym('cir',color=4),thick=1.9

psterm



end
