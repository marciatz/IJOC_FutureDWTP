import sys
import numpy as np
import xarray as xr
import DWTP_auxiliary as optb
import matplotlib.pyplot as plt

area=sys.argv[1]

match area:
    case 'global':
        xxi,xxf,yyi,yyf=[0,360,50,-50]
    case 'SAm':
        xxi,xxf,yyi,yyf=[275.,330.,15.0,-60.0]
    case 'SAfr':
        xxi,xxf,yyi,yyf=[10.0,55.0,0.0,-35.0]
    case 'EAfr':
        xxi,xxf,yyi,yyf=[25.0,55.0,15.0,-10.0]
    case 'Afr':
        xxi,xxf,yyi,yyf=[-25,57.0,40.0,-46.0]
    case _:
        print('Uknown domain')
        sys.exit()

maindir='/gws/ssde/j25a/firstrains/OnsetPeriod/data/'
pltdir='/gws/ssde/j25a/firstrains/OnsetPeriod/plots/SensTests/'

if(area == 'SAm'):
    mdls=['CHIRPS-0','UKMO-control','UKMO-rcp85','UKMO-n512','UKMO-n512rcp85','h216-hist','h216-ssp585']
if(area in ['Afr','SAfr']):
    mdls=['CHIRPS-0','CP4-CP4','CP4-CP4.fc','CP4-GA25','CP4-GA25.fc','h216-hist','h216-ssp585']

ccol=['k','#2255AA','#2255AA','#006644','#006644','#882288','#882288']
mdlCol=dict(zip(mdls,ccol))
aa='Afr' if area=='SAfr' else area

#==========
#TEST 1 - CLIMATOLOGICAL ENVELOPE - LENGTH OF THE BUTTERWORTH FILTER
tst='BWFlt'
hstcol=['#6F4C9B', '#4D8AC6', '#60AB9E','#A6BE54','#E49C39','#DF4828','#721E17']

aveMdl={};sdMdl={}
for mdl in mdls:
    dset=mdl.split('-')[0]
    prd=mdl.split('-')[1]
    ndays=365 if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG'] else 360

    if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG']: 
        if prd == '0':
            envfile=maindir+dset+'/'+aa+"_OnsPerClimEnvelope_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+aa+"_"+prd+"_OnsPerClimEnvelope_"+tst+"_025dg.aw.nc"
    else: 
        if ((dset == 'bam') & (prd=='all')):
            pass

        if dset in ['UKMO','SAAG']:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsPerClimEnvelope_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsPerClimEnvelope_"+tst+".nc"

    clmEnv=xr.open_dataarray(envfile)
    if(area == 'SAfr'):
        if(dset in ['h216','h96']): 
            clmEnv=clmEnv.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
        else:
            clmEnv=clmEnv.sel(lat=slice(yyi,yyf),lon=slice(xxi,xxf))

    tmp=xr.where(clmEnv.sel(sea=0,Stats='FirstDay')>120,\
                 clmEnv.sel(sea=0,Stats='FirstDay'),clmEnv.sel(sea=0,Stats='FirstDay')+ndays)
    
    bse=tmp.sel(BWFlt=60)
    flt=tmp.stack(z=['lat','lon']).dropna("z")
    dff=(tmp-bse).stack(z=['lat','lon']).dropna("z")


    #Histograms comparing the baseline and the sensitivity test
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,bw in enumerate(clmEnv.BWFlt.values):
        hst=np.histogram(flt.sel(BWFlt=bw).values,bins=np.arange(1,370,10))
        plt.step(hst[1][:-1],100*(hst[0]/flt.z.size),hstcol[ii],linewidth=.75,where='post')
    hst=np.histogram(flt.sel(BWFlt=60).values,bins=np.arange(1,370,10))
    plt.step(hst[1][:-1],100*(hst[0]/flt.z.size),'k',linewidth=1,where='post')

    plt.xlim(0,370)
    plt.xticks(fontsize=8); plt.xlabel('DOY',fontsize=8)
    plt.ylim(0,20)
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

    outfile=pltdir+'Hist_ClmEnvBWFlt-FirstDay_'+dset+'-'+prd+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()

    #Histogram with the difference btw the baseline and the sensitivity test
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,bw in enumerate(clmEnv.BWFlt.values):
        if bw == 60: continue
        hst=np.histogram(dff.sel(BWFlt=bw).values,bins=np.arange(-40,41,1))
        plt.step(hst[1][:-1],100*(hst[0]/dff.z.size),hstcol[ii],linewidth=.75,where='post')
    plt.plot([0,0],[0,20],'k--',linewidth=.3)

    plt.xlim(-40,40)
    plt.xticks(fontsize=8); plt.xlabel('$\Delta$ days',fontsize=8)
    plt.ylim(0,20)
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

    outfile=pltdir+'HistDiff_ClmEnvBWFlt-FirstDay_'+dset+'-'+prd+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()

    #Mean change per filter length
    aveMdl[mdl]=dff.mean(dim='z')
    sdMdl[mdl]=dff.std(dim='z')

    clmEnv.close()

fig=plt.figure(figsize=(4,3))
ax1=fig.add_axes([.12,.13,.86,.84])
lst=['-','-','--','-','--','-','--']
for ii,mdl in enumerate(mdls):
    plt.plot(aveMdl[mdl].BWFlt.values,aveMdl[mdl].values,color=mdlCol[mdl],linewidth=0.5,linestyle=lst[ii])

plt.plot(aveMdl['CHIRPS-0'].BWFlt.values,\
             np.repeat(0,aveMdl['CHIRPS-0'].BWFlt.size),'k--',linewidth=.3)
plt.xticks(aveMdl['CHIRPS-0'].BWFlt.values,fontsize=8)
plt.xlabel('Butterworth filter length (days)',fontsize=8)
plt.ylim(-30,30)
plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Mean $\Delta$ days',fontsize=8)

outfile=pltdir+'MeanDiff_ClmEnvBWFlt-FirstDay_'+area+'.png'
plt.savefig(outfile, format='png', dpi=300)
plt.close()


#=======
#Test 2 - CLIMATOLOGICAL ENVELOPE - THRESHOLD FOR THE ONSET DEFINITION
tst='ShrYr'
hstcol=['#6F4C9B', '#4D8AC6', '#A6BE54','#E49C39','#721E17']

aveMdl={};sdMdl={}
for mdl in mdls:
    dset=mdl.split('-')[0]
    prd=mdl.split('-')[1]
    ndays=365 if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG'] else 360

    if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG']: 
        if prd == '0':
            envfile=maindir+dset+'/'+aa+"_OnsPerClimEnvelope_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+aa+"_"+prd+"_OnsPerClimEnvelope_"+tst+"_025dg.aw.nc"
    else: 
        if ((dset == 'bam') & (prd=='all')):
            pass

        if dset in ['UKMO','SAAG']:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsPerClimEnvelope_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsPerClimEnvelope_"+tst+".nc"

    clmEnv=xr.open_dataarray(envfile)
    if(area == 'SAfr'):
        if(dset in ['h216','h96']): 
            clmEnv=clmEnv.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
        else:
            clmEnv=clmEnv.sel(lat=slice(yyi,yyf),lon=slice(xxi,xxf))


    tmp=xr.where(clmEnv.sel(sea=0,Stats='FirstDay')>120,\
                 clmEnv.sel(sea=0,Stats='FirstDay'),clmEnv.sel(sea=0,Stats='FirstDay')+ndays)
    
    bse=tmp.sel(ShrYr=50)
    flt=tmp.stack(z=['lat','lon']).dropna("z")
    dff=(tmp-bse).stack(z=['lat','lon']).dropna("z")

    #Histograms comparing the baseline and the sensitivity test
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,ff in enumerate(clmEnv.ShrYr.values):
        hst=np.histogram(flt.sel(ShrYr=ff).values,bins=np.arange(1,370,10))
        plt.step(hst[1][:-1]+5,100*(hst[0]/flt.z.size),hstcol[ii],linewidth=.75,where='post')
    hst=np.histogram(flt.sel(ShrYr=50).values,bins=np.arange(1,370,10))
    plt.step(hst[1][:-1],100*(hst[0]/flt.z.size),'k',linewidth=1,where='post')

    plt.xlim(0,370)
    plt.xticks(fontsize=8); plt.xlabel('DOY',fontsize=8)
    plt.ylim(0,20)
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

    outfile=pltdir+'Hist_ClmEnvShrYr-FirstDay_'+dset+'-'+prd+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()

    #Histogram with the difference btw the baseline and the sensitivity test
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,ff in enumerate(clmEnv.ShrYr.values):
        if ff == 50: continue
        hst=np.histogram(dff.sel(ShrYr=ff).values,bins=np.arange(-40,41,1))
        plt.step(hst[1][:-1],100*(hst[0]/dff.z.size),hstcol[ii],linewidth=.75,where='post')
    plt.plot([0,0],[0,30],'k--',linewidth=.3)

    plt.xlim(-40,40)
    plt.xticks(fontsize=8); plt.xlabel('$\Delta$ days',fontsize=8)
    plt.ylim(0,30)
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

    outfile=pltdir+'HistDiff_ClmEnvShrYr-FirstDay_'+dset+'-'+prd+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()

    #Mean change per filter length
    aveMdl[mdl]=dff.mean(dim='z')
    sdMdl[mdl]=dff.std(dim='z')

    clmEnv.close()

fig=plt.figure(figsize=(4,3))
ax1=fig.add_axes([.12,.13,.86,.84])
lst=['-','-','--','-','--','-','--']
for ii,mdl in enumerate(mdls):
    plt.plot(aveMdl[mdl].ShrYr.values,aveMdl[mdl].values,color=mdlCol[mdl],linewidth=0.5,linestyle=lst[ii])
lst=['-','-','--','-','--','-','--']


plt.plot(aveMdl['CHIRPS-0'].ShrYr.values,\
             np.repeat(0,aveMdl['CHIRPS-0'].ShrYr.size),'k--',linewidth=.3)
plt.xticks(aveMdl['CHIRPS-0'].ShrYr.values,fontsize=8)
plt.xlabel('Percentage of Agreement (%)',fontsize=8)
plt.ylim(-20,20)
plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Mean $\Delta$ days',fontsize=8)

outfile=pltdir+'MeanDiff_ClmEnvShrYr-FirstDay_'+area+'.png'
plt.savefig(outfile, format='png', dpi=300)
plt.close()

#=======
#TEST 3 - DWTP CHARACTERISTICS - CLIMATOLOGICAL ENVELOPE

tst='ClmEnv'
hstcol=['#6F4C9B', '#4D8AC6', '#60AB9E','#A6BE54','#E49C39','#DF4828','#721E17']

aveMdl={};sdMdl={}
for mdl in mdls:
    aveMdl[mdl]={};sdMdl[mdl]={}
    dset=mdl.split('-')[0]
    prd=mdl.split('-')[1]
    ndays=365 if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG'] else 360

    if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG']: 
        if prd == '0':
            envfile=maindir+dset+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+aa+"_"+prd+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
    else: 
        if ((dset == 'bam') & (prd=='all')):
            pass

        if dset in ['UKMO','SAAG']:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsetPeriod.V2_"+tst+".nc"

    ons=xr.open_dataarray(envfile)
    if(area == 'SAfr'):
        if(dset in ['h216','h96']): 
            ons=ons.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
        else:
            ons=ons.sel(lat=slice(yyi,yyf),lon=slice(xxi,xxf))

    xlm={'FirstDay':[0,370,10],'Length':[20,140,5],'OnsIntensity':[-10,30,1],
         'PrecipRate':[0,14,0.5],'FrcDryDays':[0,100,5]}
    ylm={'FirstDay':15,'Length':12,'OnsIntensity':15,
         'PrecipRate':12,'FrcDryDays':15}
    xdflm={'FirstDay':[-100,100,10],'Length':[-20,20,1],'OnsIntensity':[-5,5,0.2],
           'PrecipRate':[-2,2,0.1],'FrcDryDays':[-20,20,1]}
    ydflm={'FirstDay':90,'Length':50,'OnsIntensity':80,
           'PrecipRate':40,'FrcDryDays':40}
    lbl={'FirstDay':'DOY','Length':'Length (days)','OnsIntensity':'Intensity ($x10^{-2} mm.day^{-2}$)',
         'PrecipRate':'Precip Rate (mm/day)','FrcDryDays':'Dry Days (%)'}
    lbdf={'FirstDay':'$\Delta$ days','Length':'$\Delta$ days',\
          'OnsIntensity':'$\Delta$ Intensity ($x10^{-2} mm.day^{-2}$)',
          'PrecipRate':'$\Delta$ Precip Rate (mm/day)','FrcDryDays':'$\Delta$ Dry Days (%)'}
   
    tmp=ons.sel(sea=0)

    tmp.loc[{'nstats':'FirstDay'}]=xr.where(tmp.sel(nstats='FirstDay')>120,\
                                            tmp.sel(nstats='FirstDay'),\
                                            tmp.sel(nstats='FirstDay')+ndays)
    tmp=xr.where(np.isinf(tmp),np.nan,tmp)
    
    bse=tmp.sel(ClmEnv=0)
    dff=(tmp-bse).stack(z=['lat','lon','years']).dropna("z")
    
    #Mean change per filter length
    aveMdl[mdl]=dff.mean(dim='z')
    sdMdl[mdl]=dff.std(dim='z')

    ons.close()

    for nvar in ['FirstDay','Length','OnsIntensity','PrecipRate','FrcDryDays']:

        flt=tmp.sel(nstats=nvar).stack(z=['lat','lon','years']).dropna("z")
        dff2=dff.sel(nstats=nvar).copy()
        if(nvar == 'OnsIntensity'): flt=flt*100; dff2=dff2*100

        #Histograms comparing the baseline and the sensitivity test
        fig=plt.figure(figsize=(4,3))
        ax1=fig.add_axes([.12,.13,.86,.84])
        for ii,ff in enumerate(flt.ClmEnv.values):
            hst=np.histogram(flt.sel(ClmEnv=ff).values,\
                             bins=np.arange(xlm[nvar][0],xlm[nvar][1]+xlm[nvar][2],xlm[nvar][2]))
            plt.step(hst[1][:-1],100*(hst[0]/flt.z.size),hstcol[ii],linewidth=.75,where='post')
        hst=np.histogram(flt.sel(ClmEnv=0).values,\
                         bins=np.arange(xlm[nvar][0],xlm[nvar][1]+xlm[nvar][2],xlm[nvar][2]))
        plt.step(hst[1][:-1],100*(hst[0]/flt.z.size),'k',linewidth=0.75,where='post')
    
        plt.xlim(xlm[nvar][0],xlm[nvar][1])
        plt.xticks(fontsize=8); plt.xlabel(lbl[nvar],fontsize=8)
        plt.ylim(0,ylm[nvar])
        plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

        outfile=pltdir+'Hist_DWTPClmEnv-'+nvar+'_'+dset+'-'+prd+'_'+area+'.png'
        plt.savefig(outfile, format='png', dpi=300)
        plt.close()

        #Histogram with the difference btw the baseline and the sensitivity test
        fig=plt.figure(figsize=(4,3))
        ax1=fig.add_axes([.12,.13,.86,.84])
        for ii,ff in enumerate(flt.ClmEnv.values):
            if ff == 0: continue
            hst=np.histogram(dff2.sel(ClmEnv=ff).values,\
                             bins=np.arange(xdflm[nvar][0],xdflm[nvar][1]+xdflm[nvar][2],xdflm[nvar][2]))
            plt.step(hst[1][:-1],100*(hst[0]/dff2.z.size),hstcol[ii],linewidth=.75,where='post')
        plt.plot([0,0],[0,ydflm[nvar]],'k--',linewidth=.3)

        plt.xlim(xdflm[nvar][0],xdflm[nvar][1])
        plt.xticks(fontsize=8); plt.xlabel(lbl[nvar],fontsize=8)
        plt.ylim(0,ydflm[nvar])
        plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

        outfile=pltdir+'HistDiff_DWTPClmEnv-'+nvar+'_'+dset+'-'+prd+'_'+area+'.png'
        plt.savefig(outfile, format='png', dpi=300)
        plt.close()

ylm={'FirstDay':(-8,8),'Length':(-5,5),'OnsIntensity':(-0.5,0.5),
     'PrecipRate':(-0.5,0.5),'FrcDryDays':(-4,4)}
lbl={'FirstDay':'$\Delta$ First Day (days)','Length':'$\Delta$ Length (days)',\
     'OnsIntensity':'$\Delta$ Intensity ($x10^{-2} mm.day^{-2}$)',
     'PrecipRate':'$\Delta$ Precip Rate (mm/day)','FrcDryDays':'$\Delta$ Dry Days (%)'}
lst=['-','-','--','-','--','-','--']

for nvar in ['FirstDay','Length','OnsIntensity','PrecipRate','FrcDryDays']:
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,mdl in enumerate(mdls):
        tmp2=aveMdl[mdl].sel(nstats=nvar)
        if(nvar == 'OnsIntensity'): tmp2=tmp2*100
        plt.plot(tmp2.ClmEnv.values,tmp2.values,color=mdlCol[mdl],linewidth=0.5,linestyle=lst[ii])

    plt.plot(aveMdl['CHIRPS-0'].ClmEnv.values,\
             np.repeat(0,aveMdl['CHIRPS-0'].ClmEnv.size),'k--',linewidth=.3)
    plt.xticks(aveMdl['CHIRPS-0'].ClmEnv.values,fontsize=8)
    plt.xlabel('Clim Envelope ($\Delta$days)',fontsize=8)
    plt.ylim(ylm[nvar][0],ylm[nvar][1])
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel(lbl[nvar],fontsize=8)

    outfile=pltdir+'MeanDiff_DWTPClmEnv-'+nvar+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()

#Changes in FUT-PD
if area == 'SAm':
    cmpmdls={'PD':['h216-hist','UKMO-control','UKMO-n512'],
             'Fut':['h216-ssp585','UKMO-rcp85','UKMO-n512rcp85']}

if area in ['SAfr','Afr']:
    cmpmdls={'PD':['h216-hist','CP4-CP4','CP4-GA25'],
             'Fut':['h216-ssp585','CP4-CP4.fc','CP4-GA25.fc']}
AveDiff={}
for pdm,fut in zip(cmpmdls['PD'],cmpmdls['Fut']):

    dset=pdm.split('-')[0]
    ndays=360

    AveDiff[fut]={}

    if dset in ['UKMO','SAAG']:
        pdfile=maindir+dset+'/'+pdm.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        futfile=maindir+dset+'/'+fut.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        
    else:
        pdfile=maindir+dset+'/'+pdm.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+".nc"
        futfile=maindir+dset+'/'+fut.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+".nc"

    with xr.open_dataarray(pdfile) as inp:
        onsPD=inp.sel(sea=0)
    del inp
    with xr.open_dataarray(futfile) as inp:
        onsFUT=inp.sel(sea=0)
    del inp
    
    if(area == 'SAfr'):
        if(dset in ['h216','h96']): 
            onsPD=onsPD.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
            onsFUT=onsFUT.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
        else:
            onsFUT=onsFUT.sel(lat=slice(yyi,yyf),lon=slice(xxi,xxf))

    onsPD.loc[{'nstats':'FirstDay'}]=xr.where(onsPD.sel(nstats='FirstDay')>120,\
                                              onsPD.sel(nstats='FirstDay'),\
                                              onsPD.sel(nstats='FirstDay')+ndays)
    onsPD=xr.where(np.isinf(onsPD),np.nan,onsPD)
    
    onsFUT.loc[{'nstats':'FirstDay'}]=xr.where(onsFUT.sel(nstats='FirstDay')>120,\
                                               onsFUT.sel(nstats='FirstDay'),\
                                               onsFUT.sel(nstats='FirstDay')+ndays)
    onsFUT=xr.where(np.isinf(onsFUT),np.nan,onsFUT)

    avePD=onsPD.mean(dim='years')
    avePD.loc[{'nstats':'FirstDay'}]=optb.CircularMedian(onsPD.sel(nstats='FirstDay'),ndays)

    aveFUT=onsFUT.mean(dim='years')
    aveFUT.loc[{'nstats':'FirstDay'}]=optb.CircularMedian(onsFUT.sel(nstats='FirstDay'),ndays)
   
    dff=aveFUT-avePD

    dff=(dff-dff.sel(ClmEnv=0)).stack(z=['lat','lon']).dropna("z")
    
    #Mean change per filter length
    AveDiff[fut]=dff.mean(dim='z')



ylm={'FirstDay':(-2,2),'Length':(-1,1),'OnsIntensity':(-0.2,0.2),
     'PrecipRate':(-0.2,0.2),'FrcDryDays':(-1,1)}
lbl={'FirstDay':'$\Delta$ First Day (days)','Length':'$\Delta$ Length (days)',\
     'OnsIntensity':'$\Delta$ Intensity ($x10^{-2} mm.day^{-2}$)',
     'PrecipRate':'$\Delta$ Precip Rate (mm/day)','FrcDryDays':'$\Delta$ Dry Days (%)'}

for nvar in ['FirstDay','Length','OnsIntensity','PrecipRate','FrcDryDays']:
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,mdl in enumerate(AveDiff.keys()):
        tmp2=AveDiff[mdl].sel(nstats=nvar)
        if(nvar == 'OnsIntensity'): tmp2=tmp2*100
        plt.plot(tmp2.ClmEnv.values,tmp2.values,color=mdlCol[mdl],linewidth=0.5)

    plt.plot(AveDiff[cmpmdls['Fut'][0]].ClmEnv.values,\
             np.repeat(0,AveDiff[cmpmdls['Fut'][0]].ClmEnv.size),'k--',linewidth=.3)
    plt.xticks(AveDiff[cmpmdls['Fut'][0]].ClmEnv.values,fontsize=8)
    plt.xlabel('Clim Envelope ($\Delta$days)',fontsize=8)
    plt.ylim(ylm[nvar][0],ylm[nvar][1])
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel(lbl[nvar],fontsize=8)

    outfile=pltdir+'MeanDiff_DWTPClmEnv-FUT-PD_'+nvar+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()



#==========
#TEST 4 - DWTP SENSITIVITY TO THE DRY SPELL DEFINITION

tst='DrySpl'
hstcol=['#6F4C9B', '#4D8AC6', '#60AB9E','#A6BE54','#E49C39','#DF4828','#721E17']

aveMdl={};sdMdl={}
for mdl in mdls:
    aveMdl[mdl]={};sdMdl[mdl]={}
    dset=mdl.split('-')[0]
    prd=mdl.split('-')[1]
    ndays=365 if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG'] else 360

    if dset in ['CHIRPS','Br-DWGD','ERA5','IMERG']: 
        if prd == '0':
            envfile=maindir+dset+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+aa+"_"+prd+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
    else: 
        if ((dset == 'bam') & (prd=='all')):
            pass

        if dset in ['UKMO','SAAG']:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        else:
            envfile=maindir+dset+'/'+prd+'/'+aa+"_OnsetPeriod.V2_"+tst+".nc"

    ons=xr.open_dataarray(envfile)
    if(area == 'SAfr'):
        if(dset in ['h216','h96']): 
            ons=ons.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
        else:
            ons=ons.sel(lat=slice(yyi,yyf),lon=slice(xxi,xxf))
    tmp=ons.sel(sea=0)

    tmp.loc[{'nstats':'FirstDay'}]=xr.where(tmp.sel(nstats='FirstDay')>120,\
                                            tmp.sel(nstats='FirstDay'),\
                                            tmp.sel(nstats='FirstDay')+ndays)
    tmp=xr.where(np.isinf(tmp),np.nan,tmp)
    
    bse=tmp.sel(DrySpl=10)
    dff=(tmp-bse).stack(z=['lat','lon','years']).dropna("z")
    
    #Mean change per filter length
    aveMdl[mdl]=dff.mean(dim='z')
    sdMdl[mdl]=dff.std(dim='z')

    ons.close()

    xlm={'FirstDay':[0,370,10],'Length':[20,140,5],'OnsIntensity':[-10,30,1],
         'PrecipRate':[0,14,0.5],'FrcDryDays':[0,100,5]}
    ylm={'FirstDay':15,'Length':12,'OnsIntensity':15,
         'PrecipRate':12,'FrcDryDays':15}
    xdflm={'FirstDay':[-100,100,10],'Length':[-20,20,1],'OnsIntensity':[-5,5,0.2],
           'PrecipRate':[-2,2,0.1],'FrcDryDays':[-20,20,1]}
    ydflm={'FirstDay':90,'Length':50,'OnsIntensity':80,
           'PrecipRate':40,'FrcDryDays':40}
    lbl={'FirstDay':'DOY','Length':'Length (days)','OnsIntensity':'Intensity ($x10^{-2} mm.day^{-2}$)',
         'PrecipRate':'Precip Rate (mm/day)','FrcDryDays':'Dry Days (%)'}
    lbdf={'FirstDay':'$\Delta$ days','Length':'$\Delta$ days',\
          'OnsIntensity':'$\Delta$ Intensity ($x10^{-2} mm.day^{-2}$)',
          'PrecipRate':'$\Delta$ Precip Rate (mm/day)','FrcDryDays':'$\Delta$ Dry Days (%)'}
    
    for nvar in ['FirstDay','Length','OnsIntensity','PrecipRate','FrcDryDays']:

        flt=tmp.sel(nstats=nvar).stack(z=['lat','lon','years']).dropna("z")
        dff2=dff.sel(nstats=nvar).copy()
        if(nvar == 'OnsIntensity'): flt=flt*100; dff2=dff2*100

        #Histograms comparing the baseline and the sensitivity test
        fig=plt.figure(figsize=(4,3))
        ax1=fig.add_axes([.12,.13,.86,.84])
        for ii,ff in enumerate(flt.DrySpl.values):
            hst=np.histogram(flt.sel(DrySpl=ff).values,\
                             bins=np.arange(xlm[nvar][0],xlm[nvar][1]+xlm[nvar][2],xlm[nvar][2]))
            plt.step(hst[1][:-1],100*(hst[0]/flt.z.size),hstcol[ii],linewidth=.75,where='post')
        hst=np.histogram(flt.sel(DrySpl=10).values,\
                         bins=np.arange(xlm[nvar][0],xlm[nvar][1]+xlm[nvar][2],xlm[nvar][2]))
        plt.step(hst[1][:-1],100*(hst[0]/flt.z.size),'k',linewidth=0.75,where='post')
    
        plt.xlim(xlm[nvar][0],xlm[nvar][1])
        plt.xticks(fontsize=8); plt.xlabel(lbl[nvar],fontsize=8)
        plt.ylim(0,ylm[nvar])
        plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

        outfile=pltdir+'Hist_DWTPDrySpl-'+nvar+'_'+dset+'-'+prd+'_'+area+'.png'
        plt.savefig(outfile, format='png', dpi=300)
        plt.close()

        #Histogram with the difference btw the baseline and the sensitivity test
        fig=plt.figure(figsize=(4,3))
        ax1=fig.add_axes([.12,.13,.86,.84])
        for ii,ff in enumerate(flt.DrySpl.values):
            if ff == 10: continue
            hst=np.histogram(dff2.sel(DrySpl=ff).values,\
                             bins=np.arange(xdflm[nvar][0],xdflm[nvar][1]+xdflm[nvar][2],xdflm[nvar][2]))
            plt.step(hst[1][:-1],100*(hst[0]/dff2.z.size),hstcol[ii],linewidth=.75,where='post')
        plt.plot([0,0],[0,ydflm[nvar]],'k--',linewidth=.3)

        plt.xlim(xdflm[nvar][0],xdflm[nvar][1])
        plt.xticks(fontsize=8); plt.xlabel(lbl[nvar],fontsize=8)
        plt.ylim(0,ydflm[nvar])
        plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel('Frequency (%)',fontsize=8)

        outfile=pltdir+'HistDiff_DWTPDrySpl-'+nvar+'_'+dset+'-'+prd+'_'+area+'.png'
        plt.savefig(outfile, format='png', dpi=300)
        plt.close()


ylm={'FirstDay':(-5,5),'Length':(-4,4),'OnsIntensity':(-0.5,0.5),
     'PrecipRate':(-0.5,0.5),'FrcDryDays':(-4,4)}
lbl={'FirstDay':'$\Delta$ First Day (days)','Length':'$\Delta$ Length (days)',\
     'OnsIntensity':'$\Delta$ Intensity ($x10^{-2} mm.day^{-2}$)',
     'PrecipRate':'$\Delta$ Precip Rate (mm/day)','FrcDryDays':'$\Delta$ Dry Days (%)'}

lst=['-','-','--','-','--','-','--']

for nvar in ['FirstDay','Length','OnsIntensity','PrecipRate','FrcDryDays']:
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,mdl in enumerate(mdls):
        tmp2=aveMdl[mdl].sel(nstats=nvar)
        if(nvar == 'OnsIntensity'): tmp2=tmp2*100
        plt.plot(tmp2.DrySpl.values,tmp2.values,color=mdlCol[mdl],linewidth=0.5,linestyle=lst[ii])

    plt.plot(aveMdl['CHIRPS-0'].DrySpl.values,\
             np.repeat(0,aveMdl['CHIRPS-0'].DrySpl.size),'k--',linewidth=.3)
    plt.xticks(aveMdl['CHIRPS-0'].DrySpl.values,fontsize=8)
    plt.xlabel('Dry Spell Threshold (%)',fontsize=8)
    plt.ylim(ylm[nvar][0],ylm[nvar][1])
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel(lbl[nvar],fontsize=8)

    outfile=pltdir+'MeanDiff_DWTPDrySpl-'+nvar+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()




AveDiff={}
for pdm,fut in zip(cmpmdls['PD'],cmpmdls['Fut']):

    dset=pdm.split('-')[0]
    ndays=360

    AveDiff[fut]={}

    if dset in ['UKMO','SAAG']:
        pdfile=maindir+dset+'/'+pdm.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        futfile=maindir+dset+'/'+fut.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+"_025dg.aw.nc"
        
    else:
        pdfile=maindir+dset+'/'+pdm.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+".nc"
        futfile=maindir+dset+'/'+fut.split('-')[1]+'/'+aa+"_OnsetPeriod.V2_"+tst+".nc"

    with xr.open_dataarray(pdfile) as inp:
        onsPD=inp.sel(sea=0)
    del inp
    with xr.open_dataarray(futfile) as inp:
        onsFUT=inp.sel(sea=0)
    del inp
    
    if(area == 'SAfr'):
        if(dset in ['h216','h96']): 
            onsPD=onsPD.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
            onsFUT=onsFUT.sel(lat=slice(yyf,yyi),lon=slice(xxi,xxf))
        else:
            onsFUT=onsFUT.sel(lat=slice(yyi,yyf),lon=slice(xxi,xxf))

    onsPD.loc[{'nstats':'FirstDay'}]=xr.where(onsPD.sel(nstats='FirstDay')>120,\
                                              onsPD.sel(nstats='FirstDay'),\
                                              onsPD.sel(nstats='FirstDay')+ndays)
    onsPD=xr.where(np.isinf(onsPD),np.nan,onsPD)
    
    onsFUT.loc[{'nstats':'FirstDay'}]=xr.where(onsFUT.sel(nstats='FirstDay')>120,\
                                               onsFUT.sel(nstats='FirstDay'),\
                                               onsFUT.sel(nstats='FirstDay')+ndays)
    onsFUT=xr.where(np.isinf(onsFUT),np.nan,onsFUT)

    avePD=onsPD.mean(dim='years')
    avePD.loc[{'nstats':'FirstDay'}]=optb.CircularMedian(onsPD.sel(nstats='FirstDay'),ndays)

    aveFUT=onsFUT.mean(dim='years')
    aveFUT.loc[{'nstats':'FirstDay'}]=optb.CircularMedian(onsFUT.sel(nstats='FirstDay'),ndays)
   
    dff=aveFUT-avePD

    dff=(dff-dff.sel(DrySpl=10)).stack(z=['lat','lon']).dropna("z")
    
    #Mean change per filter length
    AveDiff[fut]=dff.mean(dim='z')



ylm={'FirstDay':(-3,3),'Length':(-1,1),'OnsIntensity':(-0.2,0.2),
     'PrecipRate':(-0.1,0.1),'FrcDryDays':(-0.2,0.2)}
lbl={'FirstDay':'$\Delta$ First Day (days)','Length':'$\Delta$ Length (days)',\
     'OnsIntensity':'$\Delta$ Intensity ($x10^{-2} mm.day^{-2}$)',
     'PrecipRate':'$\Delta$ Precip Rate (mm/day)','FrcDryDays':'$\Delta$ Dry Days (%)'}

for nvar in ['FirstDay','Length','OnsIntensity','PrecipRate','FrcDryDays']:
    fig=plt.figure(figsize=(4,3))
    ax1=fig.add_axes([.12,.13,.86,.84])
    for ii,mdl in enumerate(AveDiff.keys()):
        tmp2=AveDiff[mdl].sel(nstats=nvar)
        if(nvar == 'OnsIntensity'): tmp2=tmp2*100
        plt.plot(tmp2.DrySpl.values,tmp2.values,color=mdlCol[mdl],linewidth=0.5)

    plt.plot(AveDiff[cmpmdls['Fut'][0]].DrySpl.values,\
             np.repeat(0,AveDiff[cmpmdls['Fut'][0]].DrySpl.size),'k--',linewidth=.3)
    plt.xticks(AveDiff[cmpmdls['Fut'][0]].DrySpl.values,fontsize=8)
    plt.xlabel('Dry Spell Threshold (%)',fontsize=8)
    plt.ylim(ylm[nvar][0],ylm[nvar][1])
    plt.yticks(fontsize=8,rotation=90,va='center'); plt.ylabel(lbl[nvar],fontsize=8)

    outfile=pltdir+'MeanDiff_DWTPDrySpl-FUT-PD_'+nvar+'_'+area+'.png'
    plt.savefig(outfile, format='png', dpi=300)
    plt.close()
