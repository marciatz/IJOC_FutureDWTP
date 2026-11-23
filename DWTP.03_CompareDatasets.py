import sys
import numpy as np
import DWTP_auxiliary as optb
import DWTP_PlotFunctions as myplt
import xarray as xr
from scipy.ndimage import gaussian_filter as gsf
import itertools as it
import scipy.stats as sps

area=sys.argv[1]

match area:
    case 'SAm':
        mdls=['CHIRPS-0','CHIRPS-UKMO','h216-hist','h216-amip','h216-ssp585',\
              'hg96-hist','hg96-amip','hg96-ssp585','bam-m1','bam-m2','bam-m3','bam-m4',\
              'bam-all','UKMO-control','UKMO-n512','UKMO-rcp85','UKMO-n512rcp85',\
              'UKMO-hindcast','UKMO-ERAI25km']
        comps={'ObsMdl':{'obs':[['CHIRPS-0'],['CHIRPS-UKMO']],\
                         'mdl':[['h216-hist','h216-amip','hg96-hist','hg96-amip',\
                                 'bam-m1','bam-m2','bam-m3','bam-m4','bam-all'],\
                                ['UKMO-control','UKMO-hindcast','UKMO-n512','UKMO-ERAI25km']]},
               'FutPD':{'PD':[['h216-hist'],['hg96-hist'],['UKMO-control'],['UKMO-n512']],
                        'Fut':[['h216-ssp585'],['hg96-ssp585'],['UKMO-rcp85'],['UKMO-n512rcp85']]},
               'MdlBC':{'Mdl':[['UKMO-control'],['UKMO-hindcast'],['UKMO-rcp85']],
                        'BC':[['UKMO-n512'],['UKMO-ERAI25km'],['UKMO-n512rcp85']]}
            }
        sub={
             'WAmz':[[290,300,300,290,290],[-14,-14,-4.5,-4.5,-14]],\
             'EAmz':[[300,310,310,300,300],[-15.5,-15.5,-5.5,-5.5,-15.5]],\
             'EBr':[[310,318.5,318.5,310,310],[-18.,-18.,-10,-10,-18]],
             'HMArea':[[290,-12.5],[320,-12.5],[np.nan,7.5],'lon']
             }
             
        seqcol=['xkcd:azure','xkcd:minty green','xkcd:light peach','xkcd:cherry red']
    case aa if aa in ['SAfr','Afr']:
        mdls=['CHIRPS-0','CHIRPS-CP4','CP4-CP4','CP4-CP4.fc','CP4-RA25','CP4-RA25.fc',\
              'CP4-GA25','CP4-GA25.fc','h216-hist','h216-amip','h216-ssp585',\
              'hg96-hist','hg96-amip','hg96-ssp585','bam-m1','bam-m2','bam-m3',\
              'bam-m4','bam-all']
        comps={'ObsMdl':{'obs':[['CHIRPS-0'],['CHIRPS-CP4']],\
                         'mdl':[['bam-m1','bam-m2','bam-m3','bam-m4','bam-all',\
                                 'h216-hist','h216-amip','hg96-hist','hg96-amip'],\
                                ['CP4-CP4','CP4-RA25','CP4-GA25']]},
               'FutPD':{'PD':[['h216-hist'],['hg96-hist'],['CP4-CP4'],['CP4-RA25'],['CP4-GA25']],
                        'Fut':[['h216-ssp585'],['hg96-ssp585'],['CP4-CP4.fc'],['CP4-RA25.fc'],['CP4-GA25.fc']]},
               'MdlBC':{'Mdl':[['CP4-CP4'],['CP4-CP4.fc'],['CP4-RA25'],['CP4-RA25.fc']],
                        'BC':[['CP4-RA25'],['CP4-RA25.fc'],['CP4-GA25'],['CP4-GA25.fc']]}
            }
        sub={'ZBS':[[22,32,32,22,22],[-18,-18,-12,-12,-18]],
             'LPP':[[24,32,32,24,24],[-26,-26,-20,-20,-26]],
             'SAf':[[23,28,28,23,23],[-33,-33,-27,-27,-33]],
             'HMArea':[[18,-5],[30.,-30.],[np.nan,5],'lat']}
        seqcol=['xkcd:azure','xkcd:minty green','xkcd:light peach','xkcd:cherry red']
        if aa == 'SAfr':
            xxi,xxf,yyi,yyf=[10.0,55.0,0.0,-35.0]
        else: 
            xxi,xxf,yyi,yyf=[-18.0,52.0,15.0,-35.0]

    case _:
        print('Area not defined')
        sys.exit()

areacol=dict(zip(sub.keys(),seqcol))

maindir='/gws/ssde/j25a/firstrains/OnsetPeriod/data/'
pltdir='/gws/ssde/j25a/firstrains/OnsetPeriod/plots/'

# Load the data
onsets={}
for mdl in mdls:
    print(mdl)
    onsets[mdl]={}

    ndays=365 if(mdl.split('-')[0] in ['CHIRPS','bam']) else 360
    dd=[0,31,59,90,120,151,181,212,243,273,304,334] if ndays > 360 else list(range(0,ndays,30))
    frdd=np.asarray(dd)/float(ndays)
    dtxt=["1/"+str(mm) for mm in range(1,13)]

    #loading onset data
    indir=maindir+mdl.split('-')[0]+'/'
    prd=mdl.split('-')[1]
    if mdl.split('-')[0] in ['CHIRPS','Br-DWGD']:
        infile=indir+area+'_OnsetPeriod.V2_025dg.aw.nc' if prd=='0' else\
               indir+area+'_'+prd+'_OnsetPeriod.V2_025dg.aw.nc'
        b19file=indir+area+'_OnsetDemise_Bombardi2019.nc' if prd=='0' else\
                indir+area+'_'+prd+'_OnsetDemise_Bombardi2019_025dg.aw.nc'
        if(area == 'SAfr'):
            infile=indir+'Afr_OnsetPeriod.V2_025dg.aw.nc' if prd=='0' else\
                   indir+'Afr_'+prd+'_OnsetPeriod.V2_025dg.aw.nc'
            b19file=indir+'Afr_OnsetDemise_Bombardi2019.nc' if prd=='0' else\
                    indir+'Afr'+'_'+prd+'_OnsetDemise_Bombardi2019_025dg.aw.nc'
    else:
        if (mdl == 'bam-all'):
            pass
        else:
            indir=indir+prd+'/'

        if mdl.split('-')[0] in ['UKMO']:
            infile=indir+area+'_OnsetPeriod.V2_025dg.aw.nc'
            b19file=indir+area+'_OnsetDemise_Bombardi2019_025dg.aw.nc'
        else: 
            if area == 'SAfr': 
                infile=indir+'Afr_OnsetPeriod.V2.nc'
                b19file=indir+'Afr_OnsetDemise_Bombardi2019.nc'
            else: 
                infile=indir+area+'_OnsetPeriod.V2.nc'
                b19file=indir+area+'_OnsetDemise_Bombardi2019.nc'
        
    ons=xr.load_dataarray(infile)
    b19=xr.load_dataset(b19file).onset 
    b19=b19.rename({'latitude':'lat','longitude':'lon','year':'years'})

    if (ons.lat.diff(dim='lat') > 0).all():
        ons=ons.sel(lat=slice(None,None,-1)) #reversing latitudes
        b19=b19.sel(lat=slice(None,None,-1)) #reversing latitudes
        
    if(area in ['Afr','SAfr']):
        ons=ons.sel(lon=slice(xxi,xxf),lat=slice(yyi,yyf))
        b19=b19.sel(lon=slice(xxi,xxf),lat=slice(yyi,yyf))

    #Nr of valid data in each grid point
    nrvld=ons.sel(nstats='FirstDay').notnull(keep_attrs=True).sum(dim='years')
    nrvld=nrvld.where(nrvld > 0,np.nan)

    #Nr of rainy seasons
    nrrainy=xr.where(nrvld > 0,1,0.).sum(dim='sea')
    nrrainy=nrrainy.where(nrrainy > 0,np.nan)

    avetmp=ons.sel(nstats='FirstDay').median(dim='years',skipna=True)
    avetmp=avetmp.where(avetmp < 365,avetmp-365)

    if(area == 'SAm'):
        tmp=xr.where((avetmp.sel(sea=0) > avetmp.sel(sea=1)) | \
                     (avetmp.sel(sea=1).isnull()),0,1).stack(z=('lat','lon'))
    else:
        tmp=xr.where(avetmp.sel(sea=1) > avetmp.sel(sea=0),1,0).stack(z=('lat','lon'))

    pxl=tmp.where(tmp == 1,drop=True).coords['z'].values
    
    tmp=ons.copy(deep=True)
    for yy,xx in pxl:
        tmp.loc[{'lat':yy,'lon':xx}]=ons.sel(lat=yy,lon=xx).roll(sea=1)
    
    ons=tmp
    del tmp

    #Preparing the onset date stats
    b19=b19.expand_dims({'nstats':1,'sea':2}).transpose('lat','lon','years','sea','nstats') #expanding to match the years dimension
    b19=b19.assign_coords({'nstats':['Onset'],'sea':[0,1],'years':ons.years.values}).copy()
    b19.loc[{'sea':1}]=np.nan
    b19=b19%ndays

    #concatenating b19 to the onset data array
    ons=xr.concat([ons.sel(nstats=['FirstDay']),b19,ons.isel(nstats=slice(1,None))],dim='nstats')
    b19.close()
    del b19

    #loading CHIRPS mask to remove ocean - only for h216 and h96
    if mdl.split('-')[0] in ['h216','hg96','bam']:
        with xr.open_dataarray(maindir+'mask_CHIRPS_0.25dg.nc') as inp:
            msk=inp.sel(longitude=slice(ons.lon[0].values,ons.lon[-1].values),\
                        latitude=slice(ons.lat[0].values,ons.lat[-1].values)).isel(time=0)
        del inp
        #regrid msk to model's grid
        msk=msk.rename({'latitude':'lat','longitude':'lon'})
        msk=msk.interp_like(ons.isel(years=0,sea=0,nstats=0),method='nearest')  
        ons=ons*msk
        del msk


    #Median and average onset period stats
    ave=ons.mean(dim='years',skipna=True)
    ave.loc[{'nstats':['FirstDay','Onset']}]=optb.CircularMedian(ons.sel(nstats=['FirstDay','Onset']),ndays)

    onsets[mdl]['NrYears']=nrvld
    onsets[mdl]['NrRainy']=nrrainy
    onsets[mdl]['MeanStats']=ave
    onsets[mdl]['VarStats']=ons.var(dim='years',skipna=True)
    onsets[mdl]['VarStats'].loc[{'nstats':['FirstDay','Onset']}]=optb.CircularVar(ons.sel(nstats=['FirstDay','Onset']),ndays)
 
    #Variable comparing B19 onset date and first day
    tmpons=ons.sel(nstats=['Onset','FirstDay'])
    tmpons=tmpons.where(tmpons > 180,tmpons+365)

    diffons=xr.where(((tmpons.sel(sea=0,nstats='Onset')>=tmpons.sel(sea=0,nstats='FirstDay')) &\
                      (tmpons.sel(sea=0,nstats='Onset')<=(tmpons.sel(sea=0,nstats='FirstDay')+\
                                             ons.sel(sea=0,nstats='Length')))),0,np.nan)
    diffons=diffons.where(tmpons.sel(sea=0,nstats='Onset') > tmpons.sel(sea=0,nstats='FirstDay'),-1) 
    diffons=diffons.where(tmpons.sel(sea=0,nstats='Onset') < (\
                          tmpons.sel(sea=0,nstats='FirstDay')+ons.sel(sea=0,nstats='Length')),1) 
    diffons=diffons.where(ave.sel(sea=0,nstats='FirstDay').notnull(),np.nan) #Masking ocean out
    diffons=diffons.where(ave.sel(nstats='Onset',sea=0).notnull(),np.nan) #Masking areas without B19

    ons.close()

    #Plotting maps
    #Nr of dry-to-wet transition seasons
    indat={
        'shd':{'lat':ons.lat.values,'lon':ons.lon.values,'ave':nrrainy.values},
        'cbar':{'pst':[.25,0.75],'txt':['One','Two']},
    }

    outf=pltdir+'Map_NrRainy_'+area+'_'+mdl
    myplt.Map(area,indat,'xkcd:charcoal grey','tab10',(1,3,1,1),'Nr Rainy Seasons',outf)

    #Percentage of years with B19 onset date outside the DWTP: negative values indicate that the 
    #B19 onset date is before the start of the DWTP, while positive values indicate that the 
    #B19 onset date is after the end of the DWTP
    tmp=gsf(diffons.mean('years').values,2)
    tx1='Ons > DWTP: '+str(np.round(((diffons > 0).mean())*100,1).item())+'%'
    ps1=(304,-60) if area == 'SAm' else (40,-35)
    tx2='Ons < DWTP: '+str(np.round(((diffons < 0).mean())*100,1).item())+'%'
    ps2=(304,-63) if area == 'SAm' else (40,-37)

    indat={
        'shd':{'lat':ons.lat.values,'lon':ons.lon.values,'ave':(100*diffons.mean('years')).values},
        'ctl':{'lat':ons.lat.values,'lon':ons.lon.values,'ave':100*tmp,\
               'lvl':[-100,-75,-50,-25,-10,-5,5,10,25,50,75,100],'col':'navy'},
        'txt':{'lat':ps1[1],'lon':ps1[0],'txt':tx1},
        'txt2':{'lat':ps2[1],'lon':ps2[0],'txt':tx2},
    }

    outf=pltdir+'Map_OnsetFirstDay_'+area+'_'+mdl
    myplt.Map(area,indat,'xkcd:charcoal grey','cmo.balance',(-100,100,5,25),'Percentage of Years',outf)

    #preparing the coords to plot the HM area
    subs=sub.copy()
    if('HMArea' in subs.keys()):
        aa='HMArea'
        dims=[(abs(ons[subs[aa][3]][2]-ons[subs[aa][3]][1]).item()),subs[aa][2][1]]
        hv,mskhv=myplt.Hovmoller(subs[aa][0],subs[aa][1],ons.sel(sea=0,nstats='FirstDay'),dims,\
                                 subs[aa][3],circular=ndays)
      
        yy,xx=(mskhv == 1).nonzero()
        if((subs[aa][0][1]==subs[aa][1][1]) | \
           (subs[aa][0][0]==subs[aa][1][0])):
            lstx=[xx[0],xx[-1],xx[-1],xx[0],xx[0]]
            lsty=[yy[-1],yy[-1],yy[0],yy[0],yy[-1]]
        else:
            lstx=[np.min(xx),xx[0],np.max(xx),xx[-1],np.min(xx)]
            lsty=[yy[(xx==np.min(xx)).nonzero()[0][0]],yy[0],yy[(xx==np.max(xx)).nonzero()[0][-1]],\
                yy[-1],yy[(xx==np.min(xx)).nonzero()[0][0]]]

        subs[aa]=[list(ons.lon[lstx].values),list(ons.lat[lsty].values)]

    #Maps with median onset
    for ss in range(int(np.nanmax(nrrainy))):

        #Map with other statistics for the onset periods - Median
        for vv in ave.nstats.values:
            if((vv == 'Onset') & (ss==1)): 
                continue  #no B19 Onset Date in season 2 so skipping to next variable

            shd,shdave,unit,ctl=myplt.MapParam(vv,ndays)

            lvls=(shdave[0],shdave[1],shdave[3])
            caption=vv+' '+unit
            tmp=ave.sel(nstats=vv,sea=ss).values
            if(vv in ['Intensity','OnsIntensity']): tmp=tmp*10

            if(vv in ['FirstDay','Onset']):
                caption='Date'
                #Map with the First day of the season - Median
                indat={
                    'shd':{'lat':ave.lat.values,'lon':ave.lon.values,'ave':tmp},
                    'cbar':{'pst':frdd,'txt':dtxt},

                }
                if (vv == 'FirstDay'):
                    indat['ctl']={'lat':ave.lat.values,'lon':ave.lon.values,'ave':nrrainy.values,\
                                  'lvl':[4],'col':'xkcd:copper','extra':[1]}
            else:
                gf=3 if mdl.split('-')[0] in ('CHIRPS','UKMO','CP4') else 1
                if(mdl.split('-')[0] in ['hg96']): gf=0

                smth=np.where(np.isnan(tmp),0,tmp)
                smth=gsf(smth,gf)
                smth=np.where(np.isnan(tmp),np.nan,smth)
                indat={
                    'shd':{'lat':ave.lat.values,'lon':ave.lon.values,'ave':tmp},
                    'ctl':{'lat':ave.lat.values,'lon':ave.lon.values,\
                           'lvl':lvls,'col':ctl,'ave':smth}
                }

            for ia,aa in enumerate(subs):
                if(ia == 0):
                    indat['poly']={'poly':subs[aa],'ccol':areacol[aa]}
                else:
                    indat['poly'+str(ia)]={'poly':subs[aa],'ccol':areacol[aa]}

            outf=pltdir+'Map_OnsPeriod_Median_'+vv+'_S'+str(ss+1)+'_'+area+'_'+mdl
            myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdave,caption,outf)


#Compare obs and models
obslists=comps['ObsMdl']['obs'] 
mdllists=comps['ObsMdl']['mdl']

for ii in range(len(obslists)):
    for obs,mdl in it.product(obslists[ii],mdllists[ii]):

        ndmdl=365 if(mdl.split('-')[0] in ['CHIRPS','Br-DWGD','bam']) else 360
        ndobs=365
         
        if mdl.split('-')[0] in ('CHIRPS','UKMO','CP4'):
            gf=3
        elif mdl.split('-')[0] in('h216','bam'): gf=1
        else: gf=0

        #Difference in the number of rainy seasons
        obsnry=onsets[obs]['NrRainy'].copy(deep=True)
        mdlnry=onsets[mdl]['NrRainy'].copy(deep=True)

        if(mdl.split('-')[0] in ['h216','hg96','bam']):
            obsnry=obsnry.interp_like(mdlnry,method='nearest')  
            #note that here only 1 and 2 are valid values, so nearest is a better method
            
        #Ocean mask - obs to mask out the ocean
        msk=xr.where(obsnry.isnull(),np.nan,1)

        nrdiff=(mdlnry-obsnry)*msk

        lab='$\Delta$ Nr Dry-to-Wet Transitions'
        shd,shdave,unt,ctl,lvl=myplt.MapParam('NrRainy',ndmdl,ptype='diff')

        indat={
            'shd':{'lat':nrdiff.lat.values,'lon':nrdiff.lon.values,'ave':nrdiff.values},
            'cbar':{'pst':[.165,.5,.825],'txt':['-1','0','+1']},
        }
        outf=pltdir+'Map_DiffOnset_NrRainy_'+mdl+'-'+obs.split('-')[0]+'_'+area
        myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdave,lab,outf)

        for ss in range(int(np.nanmax(mdlnry))):
            for vv in onsets[obs]['MeanStats'].nstats.values:
                
                if((vv == 'Onset') & (ss==1)): 
                    continue  #no B19 Onset Date in season 2 so skipping to next variable

                obstmp=onsets[obs]['MeanStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                obsvar=onsets[obs]['VarStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                obsNyr=onsets[obs]['NrYears'].sel(sea=ss).copy(deep=True)
                mdltmp=onsets[mdl]['MeanStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                mdlvar=onsets[mdl]['VarStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                mdlNyr=onsets[mdl]['NrYears'].sel(sea=ss).copy(deep=True)

                if(mdl.split('-')[0] in ['h216','hg96','bam']):
                    obstmp=obstmp.interp_like(mdltmp,method='nearest')  
                    obsvar=obsvar.interp_like(mdlvar,method='nearest')  
                    obsNyr=obsNyr.interp_like(mdlNyr,method='nearest')  


                #Difference in the onset period stats            
                if(vv in ['FirstDay','Onset']):
                    #Transforming the dates to 365 day calendar. 
                    if(ndmdl != ndobs): mdltmp=mdltmp*(ndobs/ndmdl)
                    ndiff=(((mdltmp-obstmp+ndobs/2) % ndobs)*msk)-ndobs/2
                else:
                    ndiff=(mdltmp-obstmp)*msk
                    if(vv in ['Intensity','OnsIntensity']): ndiff=ndiff*10
                        
                #Significance test
                dfnum=((obsvar/obsNyr)+(mdlvar/mdlNyr))**2
                dfden=(obsvar**2/((obsNyr-1)*obsNyr**2))+\
                      (mdlvar**2/((mdlNyr-1)*mdlNyr**2))
                df=np.divide(dfnum,dfden)
                t=(ndiff)/np.sqrt((obsvar/obsNyr)+(mdlvar/mdlNyr))
                sgnf = sps.distributions.t.sf(np.abs(t), df) * 2

                #Plotting difference maps
                shd,shdave,unt,ctl,lvl=myplt.MapParam(vv,ndobs,ptype='diff')
                lab='$\Delta$ '+vv+' '+unt 

                if(vv in ['FirstDay','Onset']):
                    ave=np.empty_like(ndiff)
                else:
                    ave=xr.where(np.isnan(obstmp),0,obstmp)
                    ave=gsf(ave,gf) 
                    #using the gf from the mdl because for lower resolution simulations, 
                    #the obs is regrided to the mdl gridspace
                    ave=xr.where(np.isnan(obstmp),np.nan,ave)

                indat={
                    'shd':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':ndiff.values},
                    'ctl':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,\
                           'ave':(ave*msk).values,'lvl':lvl,'col':ctl},
                    'sct':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':sgnf},
                }

                outf=pltdir+'Map_DiffOnset_'+vv+'_S'+str(ss+1)+'_'+mdl+'-'+obs.split('-')[0]+'_'+area
                myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdave,lab,outf,whiteout=True)

                #Plotting differences as percentages
                if(vv in ['FirstDay','Onset']): continue

                percdiff=100*np.divide(ndiff,obstmp)*msk  
                if(vv in ['Intensity','OnsIntensity']): percdiff=100*np.divide(ndiff,obstmp*10)*msk

                shdprc=(-100,100,25,50)
                if(vv in ['FrcWetDays','LengthWetSpells','MaxWetSpell','FrcDryDays']):
                    shdprc=(-50,50,5,10)
                if(vv == 'Length'):shdprc=(-80,80,10,20)

                lab='$\Delta$ '+vv+' (%)'
                indat={
                    'shd':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,\
                            'ave':percdiff.values},
                    'ctl':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'lvl':lvl,\
                        'col':ctl,'ave':(ave*msk).values},
                    'sct':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':sgnf},
                }
                outf=pltdir+'Map_DiffOnsetPerc_'+vv+'_S'+str(ss+1)+'_'+mdl+'-'+obs.split('-')[0]+'_'+area
                myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdprc,lab,outf,whiteout=True)

#Compare future and present day
pdlists=comps['FutPD']['PD'] 
futlists=comps['FutPD']['Fut']

for ii in range(len(pdlists)):
    for pd,fut in it.product(pdlists[ii],futlists[ii]):
        
        ndays=365 if(fut.split('-')[0] in ['CHIRPS','Br-DWGD','bam']) else 360
        if pd.split('-')[1] in ('control','n512','CP4','RA25','GA25'):
            gf=2
        elif pd.split('-')[0] in('h216'): gf=1
        else: gf=0
        #Difference in the number of rainy seasons
        pdnry=onsets[pd]['NrRainy'].copy(deep=True)
        futnry=onsets[fut]['NrRainy'].copy(deep=True)

        #CP4-GA25 extends further south than GA25.fc
        if(pd.split('-')[1]=='GA25'): 
            pdnry=pdnry.interp_like(futnry)

        #Interpolating mask:
        #Using obs data with similar area coverage to mask out the ocean
        if (pd.split('-')[0] in ['UKMO','CP4']):
            mskOcean=onsets['CHIRPS-'+pd.split('-')[0]]['NrRainy'].copy(deep=True)
        else:
            mskOcean=onsets['CHIRPS-0']['NrRainy'].copy(deep=True)
        mskOcean=xr.where(mskOcean.isnull(),np.nan,1)
        msk=mskOcean.interp_like(pdnry,method='nearest')

        nrdiff=(futnry-pdnry)*msk

        lab='$\Delta$ Nr Dry-to-Wet Transitions'
        shd,shdave,unt,ctl,lvl=myplt.MapParam('NrRainy',ndays,ptype='diff')

        indat={
            'shd':{'lat':nrdiff.lat.values,'lon':nrdiff.lon.values,'ave':nrdiff.values},
            'cbar':{'pst':[.165,.5,.825],'txt':['-1','0','+1']},
        }

        tagout=pd.split('-')[1]+'_FUT-PD_' if fut.split('-')[0] in ('CP4','UKMO') else \
               fut.split('-')[0]+'_FUT-PD_'
        outf=pltdir+'Map_DiffOnset_NrRainy_'+tagout+area
        myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdave,lab,outf)

        for ss in range(int(np.nanmax(pdnry))):
            for vv in onsets[pd]['MeanStats'].nstats.values:
                
                if((vv == 'Onset') & (ss==1)): 
                    continue  #no B19 Onset Date in season 2 so skipping to next variable

                pdtmp=onsets[pd]['MeanStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                pdvar=onsets[pd]['VarStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                pdNyr=onsets[pd]['NrYears'].sel(sea=ss).copy(deep=True)
                futtmp=onsets[fut]['MeanStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                futvar=onsets[fut]['VarStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                futNyr=onsets[fut]['NrYears'].sel(sea=ss).copy(deep=True)

                #CP4-GA25 extends further south than GA25.fc
                if(pd.split('-')[1]=='GA25'): 
                    pdtmp=pdtmp.interp_like(futtmp)
                    pdvar=pdvar.interp_like(futvar)
                    pdNyr=pdNyr.interp_like(futNyr)


                #Difference in the onset period stats
                if(vv in ['FirstDay','Onset']):
                    ndiff=(((futtmp-pdtmp+ndays/2) % ndays)*msk)-ndays/2
                else:
                    ndiff=(futtmp-pdtmp)*msk
                    if(vv in ['Intensity','OnsIntensity']): ndiff=ndiff*10
                        
                #Significance test
                dfnum=((pdvar/pdNyr)+(futvar/futNyr))**2
                dfden=(pdvar**2/((pdNyr-1)*pdNyr**2))+\
                      (futvar**2/((futNyr-1)*futNyr**2))
                df=np.divide(dfnum,dfden)
                t=(ndiff)/np.sqrt((pdvar/pdNyr)+(futvar/futNyr))
                sgnf = sps.distributions.t.sf(np.abs(t), df) * 2

                #Plotting difference maps
                shd,shdave,unt,ctl,lvl=myplt.MapParam(vv,ndays,ptype='diff')
                lab='$\Delta$ '+vv+' '+unt 

                if(vv in ['FirstDay','Onset']):
                    ave=np.empty_like(ndiff); ave[:]=np.nan
                else:
                    ave=xr.where(np.isnan(pdtmp),0,pdtmp)
                    ave=gsf(ave,gf)
                    ave=xr.where(np.isnan(pdtmp),np.nan,ave)
                    if(vv in ['Intensity','OnsIntensity']): ave=ave*10


                indat={
                    'shd':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':ndiff.values},
                    'ctl':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,\
                           'ave':(ave*msk).values,'lvl':lvl,'col':ctl},
                    'sct':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':sgnf},
                }

                outf=pltdir+'Map_DiffOnset_'+vv+'_S'+str(ss+1)+'_'+tagout+area
                myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdave,lab,outf,whiteout=True)

                #Plotting differences as percentages
                if(vv in ['FirstDay','Onset']): continue

                percdiff=100*np.divide(ndiff,pdtmp)*msk  
                if(vv in ['Intensity','OnsIntensity']): percdiff=100*np.divide(ndiff,pdtmp*10)*msk

                shdprc=(-100,100,25,50)
                if(vv in ['FrcWetDays','LengthWetSpells','MaxWetSpell','FrcDryDays']):
                    shdprc=(-50,50,5,10)
                if(vv == 'Length'):shdprc=(-80,80,10,20)

                lab='$\Delta$ '+vv+' (%)'
                indat={
                    'shd':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,\
                            'ave':percdiff.values},
                    'ctl':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'lvl':lvl,\
                        'col':ctl,'ave':(ave*msk).values},
                    'sct':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':sgnf},
                }
                outf=pltdir+'Map_DiffOnsetPerc_'+vv+'_S'+str(ss+1)+'_'+tagout+area
                myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdprc,lab,outf,whiteout=True)

#Compare mdl and BC
mdllists=comps['MdlBC']['Mdl'] 
BClists=comps['MdlBC']['BC']

for ii in range(len(mdllists)):
    for mdl,bc in it.product(mdllists[ii],BClists[ii]):

        ndays=365 if(mdl.split('-')[0] in ['CHIRPS','Br-DWGD','bam']) else 360
        gf=2
        
        #Difference in the number of rainy seasons
        mdlnry=onsets[mdl]['NrRainy'].copy(deep=True)
        bcnry=onsets[bc]['NrRainy'].copy(deep=True)

        #Interpolating mask:
        #Using obs data with similar area coverage to mask out the ocean
        mskOcean=onsets['CHIRPS-'+mdl.split('-')[0]]['NrRainy'].copy(deep=True)
        mskOcean=xr.where(mskOcean.isnull(),np.nan,1)
        msk=mskOcean.interp_like(mdlnry,method='nearest')

        nrdiff=(mdlnry-bcnry)*msk

        lab='$\Delta$ Nr Dry-to-Wet Transitions'
        shd,shdave,unt,ctl,lvl=myplt.MapParam('NrRainy',ndays,ptype='diff')

        indat={
            'shd':{'lat':nrdiff.lat.values,'lon':nrdiff.lon.values,'ave':nrdiff.values},
            'cbar':{'pst':[.165,.5,.825],'txt':['-1','0','+1']},
        }
        outf=pltdir+'Map_DiffOnset_NrRainy_'+mdl+'-'+bc.split('-')[1]+'_'+area
        myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdave,lab,outf)

        for ss in range(int(np.nanmax(mdlnry))):
            for vv in onsets[mdl]['MeanStats'].nstats.values:
                
                if((vv == 'Onset') & (ss==1)): 
                    continue  #no B19 Onset Date in season 2 so skipping to next variable

                mdltmp=onsets[mdl]['MeanStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                mdlvar=onsets[mdl]['VarStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                mdlNyr=onsets[mdl]['NrYears'].sel(sea=ss).copy(deep=True)
                bctmp=onsets[bc]['MeanStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                bcvar=onsets[bc]['VarStats'].sel(nstats=vv,sea=ss).copy(deep=True)
                bcNyr=onsets[bc]['NrYears'].sel(sea=ss).copy(deep=True)

                #Difference in the onset period stats            
                if(vv in ['FirstDay','Onset']):
                    ndiff=(((mdltmp-bctmp+ndays/2) % ndays)*msk)-ndays/2
                else:
                    ndiff=(mdltmp-bctmp)*msk
                    if(vv in ['Intensity','OnsIntensity']): ndiff=ndiff*10
                        
                #Significance test
                dfnum=((bcvar/bcNyr)+(mdlvar/mdlNyr))**2
                dfden=(bcvar**2/((bcNyr-1)*bcNyr**2))+\
                      (mdlvar**2/((mdlNyr-1)*mdlNyr**2))
                df=np.divide(dfnum,dfden)
                t=(ndiff)/np.sqrt((bcvar/bcNyr)+(mdlvar/mdlNyr))
                sgnf = sps.distributions.t.sf(np.abs(t), df) * 2

                #Plotting difference maps
                shd,shdave,unt,ctl,lvl=myplt.MapParam(vv,ndays,ptype='diff')
                lab='$\Delta$ '+vv+' '+unt 

                if(vv in ['FirstDay','Onset']):
                    ave=xr.where(np.isnan(ndiff),0,ndiff)
                    ave=gsf(ave,1)
                    ave=xr.where(np.isnan(ndiff),np.nan,ave)
                else:
                    ave=xr.where(np.isnan(bctmp),0,bctmp)
                    ave=gsf(ave,1)
                    ave=xr.where(np.isnan(bctmp),np.nan,ave)

                indat={
                    'shd':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':ndiff.values},
                    'ctl':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,\
                           'ave':(ave*msk).values,'lvl':lvl,'col':ctl},
                    'sct':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':sgnf},
                }

                outf=pltdir+'Map_DiffOnset_'+vv+'_S'+str(ss+1)+'_'+mdl+'-'+bc.split('-')[1]+'_'+area
                myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdave,lab,outf,whiteout=True)

                #Plotting differences as percentages
                if(vv in ['FirstDay','Onset']): continue

                percdiff=100*np.divide(ndiff,bctmp)*msk  
                if(vv in ['Intensity','OnsIntensity']): percdiff=100*np.divide(ndiff,bctmp*10)*msk

                shdprc=(-100,100,25,50)
                if(vv in ['FrcWetDays','LengthWetSpells','MaxWetSpell','FrcDryDays']):
                    shdprc=(-50,50,5,10)
                if(vv == 'Length'):shdprc=(-80,80,5,20)

                lab='$\Delta$ '+vv+' (%)'
                indat={
                    'shd':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,\
                            'ave':percdiff.values},
                    'ctl':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'lvl':lvl,\
                        'col':ctl,'ave':(ave*msk).values},
                    'sct':{'lat':ndiff.lat.values,'lon':ndiff.lon.values,'ave':sgnf},
                }
                outf=pltdir+'Map_DiffOnsetPerc_'+vv+'_S'+str(ss+1)+'_'+mdl+'-'+bc.split('-')[1]+'_'+area
                myplt.Map(area,indat,'xkcd:charcoal grey',shd,shdprc,lab,outf,whiteout=True)
