import mocsy

# define aliases to use flxco2 functions easily
flxco2 = mocsy.gasx.flxco2
vars = mocsy.mvars.vars

# print vars.__doc__

# compute fluxes and carbonate chemistry
co2flux,co2ex,dpco2,ph,pco2,fco2,co2,hco3,co3,omegaa,omegac,betad,rhosw,p,tempis = flxco2(temp=20,sal=35,alk=2300*1.028e-3,dic=2000*1.028e-3,sil=0,phos=0,kw660=1,xco2=400,patm=1,dz1=10,optt='Tinsitu')
