ff=findfile('*@nucleo.idl',count=cnt)
psinit,0,file='snomix0610.ps',/full
!p.multi=[0,1,2]

for i=0,cnt-1 do begin
  restore,ff(i)
si28=zonedata.si28
si29=zonedata.si29
si30=zonedata.si30
si28=totdata.si28
si29=totdata.si29
si30=totdata.si30

r29=totdata.si29/totdata.si28*28/29.
r30=totdata.si30/totdata.si28*28/30.



;wf=where(finite(r18))


f=maken(0,1,1e3)
;i1=5
;i2=6 ;He/N
;wf=where(finite(r18))
;i1=where(r18 eq max(r18(wf)))
;wf=where(finite(r17))
;i2=where(r17 eq max(r17(wf)))
;i1=i1(0)
;i2=i2(0)
;m16=f*o16(i1)+(1.-f)*o16(i2)
;m17=f*o17(i1)+(1.-f)*o17(i2)
;m18=f*o18(i1)+(1.-f)*o18(i2);


;mr18=m18/m16*16./18
;mr17=m17/m16*16./17
d29=rat2del(r29,stdrat(29))
d30=rat2del(r30,stdrat(30))

plot_oo,r30/stdrat(30),r29/stdrat(29),psym=4,title=ff(i);,xrange=[1e-12,100],yrange=[1e-7,1],title=ff(i)
ver,1,lines=1
hor,1,lines=1

plots,[1e-4,1],[.341,1]

if i mod 2 then erase
endfor
!p.multi=0
psterm



end
