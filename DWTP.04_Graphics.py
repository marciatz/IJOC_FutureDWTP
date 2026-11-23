import sys
import numpy as np
import DWTP_auxiliary as optb
import DWTP_PlotFunctions as myplt
import xarray as xr
import cmocean as cmo
import csv

area=sys.argv[1]

match area:
    case 'SAm':
        mdls=['CHIRPS-0','CHIRPS-UKMO','h216-hist','h216-amip','h216-ssp585',\
              'hg96-hist','hg96-amip','hg96-ssp585','bam-m1','bam-m2','bam-m3','bam-m4',\
              'bam-all','UKMO-control','UKMO-n512','UKMO-rcp85','UKMO-n512rcp85',\
              'UKMO-hindcast','UKMO-ERAI25km']
        sub={
             'WAmz':[[290,300,300,290,290],[-14,-14,-4.5,-4.5,-14]],\
             'EAmz':[[300,310,310,300,300],[-15.5,-15.5,-5.5,-5.5,-15.5]],\
             'EBr':[[310,318.5,318.5,310,310],[-18.,-18.,-10,-10,-18]],\
             'HMArea':[[290,-12.5],[320,-12.5],[np.nan,7.5],'lon']}

        seqcol=['xkcd:azure','xkcd:minty green','xkcd:orange','magenta','xkcd:cherry red']
    case aa if aa in ['Afr','SAfr']:
        mdls=['CHIRPS-0','CHIRPS-CP4','h216-hist','h216-amip','h216-ssp585',\
              'hg96-hist','hg96-amip','hg96-ssp585','bam-m1','bam-m2','bam-m3','bam-m4',\
              'bam-all','CP4-CP4','CP4-CP4.fc','CP4-RA25','CP4-RA25.fc',\
              'CP4-GA25','CP4-GA25.fc']
        sub={'ZBS':[[22,32,32,22,22],[-18,-18,-12,-12,-18]],
             'LPP':[[24,32,32,24,24],[-26,-26,-20,-20,-26]],
             'SAf':[[23,28,28,23,23],[-33,-33,-27,-27,-33]],
             'HMArea':[[18,-5],[30.,-30.],[np.nan,5],'lat']}
        seqcol=['xkcd:azure','xkcd:minty green','xkcd:orange','magenta','yellow','xkcd:cherry red']
        if aa == 'SAfr':
            xxi,xxf,yyi,yyf=[10.0,55.0,0.0,-35.0]
        else: 
            xxi,xxf,yyi,yyf=[-18.0,52.0,15.0,-35.0]
    case _:
        print('Area not defined')
        sys.exit()

#Defining HPC used
maindir='/gws/ssde/j25a/firstrains/OnsetPeriod/data/'
pltdir='/gws/ssde/j25a/firstrains/OnsetPeriod/plots/'

# Load the data
onsets={}; trends={}
for mdl in mdls:
    onsets[mdl]={} 

    ndays=365 if(mdl.split('-')[0] in ['CHIRPS','bam']) else 360
    cal='noleap' if(mdl.split('-')[0] in ['CHIRPS','bam']) else '360_day'
    #loading onset data
    indir=maindir+mdl.split('-')[0]+'/'
    prd=mdl.split('-')[1]
    if mdl.split('-')[0] == 'CHIRPS':
        infile=indir+area+'_OnsetPeriod.V2_025dg.aw.nc' if prd=='0' else\
               indir+area+'_'+prd+'_OnsetPeriod.V2_025dg.aw.nc'
        ppfile=indir+area+'_BW60Precip_025dg.aw.nc' if prd=='0' else\
               indir+area+'_'+prd+'_BW60Precip_025dg.aw.nc'
        b19file=indir+area+'_OnsetDemise_Bombardi2019.nc' if prd=='0' else\
                indir+area+'_'+prd+'_OnsetDemise_Bombardi2019_025dg.aw.nc'
        if(area == 'SAfr'):
            infile=indir+'Afr_OnsetPeriod.V2_025dg.aw.nc' if prd=='0' else\
                   indir+'Afr_'+prd+'_OnsetPeriod.V2_025dg.aw.nc'
            ppfile=indir+'Afr_BW60Precip_025dg.aw.nc' if prd=='0' else\
                   indir+'Afr_'+prd+'_BW60Precip_025dg.aw.nc'
            b19file=indir+'Afr_OnsetDemise_Bombardi2019.nc' if prd=='0' else\
                    indir+'Afr_'+prd+'_OnsetDemise_Bombardi2019_025dg.aw.nc'
    else:
        if (mdl == 'bam-all'):
            pass
        else:
            indir=indir+prd+'/'

        if mdl.split('-')[0] in ['UKMO']:
            infile=indir+area+'_OnsetPeriod.V2_025dg.aw.nc'
            ppfile=indir+area+'_BW60Precip_025dg.aw.nc'
            b19file=indir+area+'_OnsetDemise_Bombardi2019_025dg.aw.nc'
        else: 
            if area == 'SAfr': 
                infile=indir+'Afr_OnsetPeriod.V2.nc'
                ppfile=indir+'Afr_BW60Precip.nc'
                b19file=indir+'Afr_OnsetDemise_Bombardi2019.nc'

            else: 
                infile=indir+area+'_OnsetPeriod.V2.nc'
                ppfile=indir+area+'_BW60Precip.nc'
                b19file=indir+area+'_OnsetDemise_Bombardi2019.nc'
        
    ons=xr.load_dataarray(infile)
    bwpp=xr.load_dataarray(ppfile)
    b19=xr.load_dataset(b19file).onset #.isel(years=slice(0,-1))
    b19=b19.rename({'latitude':'lat','longitude':'lon','year':'years'})#.isel(years=slice(None,-1))

    if (ons.lat.diff(dim='lat') > 0).all(): 
        ons=ons.sel(lat=slice(None,None,-1)) #reversing latitudes
        bwpp=bwpp.sel(lat=slice(None,None,-1)) #reversing latitudes
        b19=b19.sel(lat=slice(None,None,-1)) #reversing latitudes

    if(area in ['Afr','SAfr']):
        ons=ons.sel(lon=slice(xxi,xxf),lat=slice(yyi,yyf))
        bwpp=bwpp.sel(lon=slice(xxi,xxf),lat=slice(yyi,yyf))
        b19=b19.sel(lon=slice(xxi,xxf),lat=slice(yyi,yyf))

    #Sorting onsets
    avetmp=ons.sel(nstats='FirstDay').median(dim='years',skipna=True)
    avetmp=avetmp.where(avetmp < 365,avetmp-365)

    if(area == 'SAm'):
        tmp=xr.where((avetmp.sel(sea=0) > avetmp.sel(sea=1)) | \
                     (avetmp.sel(sea=1).isnull()),0,1).stack(z=('lat','lon'))
    else:
        tmp=xr.where((avetmp.sel(sea=1) < avetmp.sel(sea=0)) | \
                     (avetmp.sel(sea=1).isnull()),0,1).stack(z=('lat','lon'))

    pxl=tmp.where(tmp == 1,drop=True).coords['z'].values
    
    tmp=ons.copy(deep=True)
    for yy,xx in pxl:
        tmp.loc[{'lat':yy,'lon':xx}]=ons.sel(lat=yy,lon=xx).roll(sea=1)
    
    ons=tmp
    del tmp

    msk=xr.where((ons.sel(nstats='FirstDay').notnull()).sum(dim='years') > 0.75*ons.years.size,1,np.nan)

    #b19 data
    b19=b19.expand_dims({'nstats':1,'sea':2}).transpose('lat','lon','years','sea','nstats') #expanding to match the years dimension
    b19=b19.assign_coords({'nstats':['Onset'],'sea':[0,1],'years':ons.years.values}).copy()
    #here b19 is the same for both seasons
    b19=b19%ndays

    #concatenating b19 to the onset data array
    ons=xr.concat([ons.sel(nstats=['FirstDay']),b19,ons.isel(nstats=slice(1,None))],dim='nstats')
    b19.close()
    del b19

    #Mean over areas:
    for aa in sub:
        print(mdl,aa)

        if(aa == 'HMArea'):
            #Hovmoller plot for FirstDay and OnsetDate, so different averaging
            dims=[abs((ons[sub[aa][3]][2]-ons[sub[aa][3]][1]).item()),sub[aa][2][1]]

            #First Day
            fd,mskhv=myplt.Hovmoller(sub[aa][0],sub[aa][1],ons.sel(sea=0,nstats='FirstDay'),dims,\
                                        sub[aa][3],circular=ndays)
            #CP4 has some data in year 11 (finishes in feb of year 11) so removing it)
            vyr=[yr for yr in fd.years.values if fd.sel(years=yr).notnull().any()]
            fd=fd.sel(years=vyr)

            avearea=optb.CircularMedian(fd,ndays)

            #B19 onset day
            b19hv,mskb19=myplt.Hovmoller(sub[aa][0],sub[aa][1],ons.sel(sea=0,nstats='Onset'),dims,\
                                        sub[aa][3],circular=ndays)
            b19area=optb.CircularMedian(b19hv,ndays)
            
            shd,shdave,unit,ctl=myplt.MapParam('FirstDay',ndays)

            shd=cmo.tools.crop_by_percent(shd,20,which='min',N=None)

            outf=pltdir+'Graphs_HovmollerFirstDay_'+mdl+'_'+area
            myplt.GraphHovmoller(fd.T,avearea.values,shd,(200,340),'',sub[aa][3],outf,extra=b19area)
        
        else:
            tmparea=optb.poly2mask(ons,sub[aa])

            #Cropping data to the area bbox
            xmin,xmax=np.min(sub[aa][0]),np.max(sub[aa][0])
            ymin,ymax=np.min(sub[aa][1]),np.max(sub[aa][1])
            
            tmparea=tmparea.sel(lat=slice(ymax,ymin),lon=slice(xmin,xmax)) 

            ave=optb.CircularMedian(tmparea.sel(nstats='FirstDay'),ndays)

            #Number of existing pixels in the area
            mskarea=optb.poly2mask(msk,sub[aa])
            
            #separating the seasons (if any second season in the area)
            if(ave.sel(sea=1).notnull().any()):
                if(area == 'SAm'):
                    tmp=xr.where((ave.sel(sea=0) > ave.sel(sea=1)) | \
                                (ave.sel(sea=1).isnull()),0,1).stack(z=('lat','lon'))
                else:   
                    tmp=xr.where(ave.sel(sea=1) > ave.sel(sea=0),1,0).stack(z=('lat','lon'))

                pxl=tmp.where(tmp == 1,drop=True).coords['z'].values
                
                tmp=tmparea.copy(deep=True)
                for yy,xx in pxl:
                    tmp.loc[{'lat':yy,'lon':xx}]=tmparea.sel(lat=yy,lon=xx).roll(sea=1)
                
                tmparea=tmp
                del tmp

                #Recalculating average for each season
                ave=optb.CircularMedian(tmparea.sel(nstats='FirstDay'),ndays)

                #Refining: checking if pixels are closer to the area average of te season
                areaave=optb.CircularMedian(ave,ndays,dim=['lat','lon'])

                #Circular difference:
                fwd=(ave.sel(sea=0).expand_dims(sea=areaave.sea)-areaave) % ndays
                bwd=(areaave-ave.sel(sea=0).expand_dims(sea=areaave.sea)) % ndays
                avediff=np.minimum(fwd, bwd)

                tmp=xr.where(avediff.sel(sea=1) < avediff.sel(sea=0),1,0).stack(z=('lat','lon'))
                pxl=tmp.where(tmp == 1,drop=True).coords['z'].values

                tmp=tmparea.copy(deep=True)
                for yy,xx in pxl:
                    tmp.loc[{'lat':yy,'lon':xx}]=tmparea.sel(lat=yy,lon=xx).roll(sea=1)
                    
                tmparea=tmp
                del tmp
                
                #Recalculating average for each season
                ave=optb.CircularMedian(tmparea.sel(nstats='FirstDay'),ndays)
                ##This could be plotted as a map, focusing only on the bbox of the area polygon
        
            #Selecting only seasons with more than 25% of the area valid pixels over the area
            mxpx=mskarea.sum(dim=['lat','lon']).values  # Nr vld pxls in each season
            totpx=xr.where(mskarea.sum(dim='sea') > 0,1,0).sum(dim=['lat','lon']).item() 
            #total nr of pixels in the area
            
            vld=(mxpx >= 0.25*totpx).nonzero()[0] #Selecting only seasons with data in more than 25% of the area
            
            mskarea=mskarea.sel(sea=vld)
            tmparea=tmparea.sel(sea=vld)

            mxpx=mskarea.sum(dim=['lat','lon']).values
            tmp=xr.where(xr.where(tmparea.sel(nstats='FirstDay').notnull(),1,0).\
                                                    sum(dim=['lat','lon']) \
                                            >= (mxpx*.5),1,np.nan)
            #**Statistics over the years, averaged over the area
            onsarea=tmparea.mean(dim=['lat','lon'],skipna=True)*tmp
            onsarea.loc[{'nstats':['FirstDay','Onset']}]=\
                        optb.CircularMedian(tmparea.sel(nstats=['FirstDay','Onset']),ndays,dim=['lat','lon'])*tmp
            #CP4 has some data in year 11 (finishes in feb of year 11) so removing it)
            vyr=[yr for yr in onsarea.years.values if onsarea.sel(years=yr).notnull().any()]
            onsarea=onsarea.sel(years=vyr)
            tmparea=tmparea.sel(years=vyr)

            avearea=onsarea.quantile([.05,.25,.5,.75,.95],dim=['years'])
            avearea.loc[{'nstats':['FirstDay','Onset']}] = \
                    optb.CircularQuantile(onsarea.sel(nstats=['FirstDay','Onset']),ndays,[.05,.25,.5,.75,.95])

            onsets[mdl][aa]=avearea

            #Array with fraction of area in a transition period.
            trsmsk=xr.DataArray(coords=(xr.date_range(start=str(tmparea.years[0].item())+"-01-01",\
                                                      periods=ndays*tmparea.years.size,freq="D",calendar=cal),\
                                tmparea.lat.values,tmparea.lon.values),\
                                dims=('time','lat','lon'))
            
            mskarea=xr.where(mskarea.sel(lat=slice(ymax,ymin),lon=slice(xmin,xmax)).sum(dim='sea') >= 1,1,np.nan)
            tmp=mskarea.stack(z=('lat','lon'))
            pxl=tmp.where(tmp == 1,drop=True).coords['z'].values    
            for yy,xx in pxl:
                tmp=xr.DataArray(0,coords=((list(range(trsmsk.time.size))),),dims=('time',))

                for yr in tmparea.years.values:
                    di=(tmparea.sel(lon=xx,lat=yy,nstats='FirstDay',sea=0,years=yr).item()+\
                        (yr-tmparea.years[0].values)*ndays)
                    df=di+tmparea.sel(lon=xx,lat=yy,nstats='Length',sea=0,years=yr).item()
                    if(np.isfinite(di)): tmp.loc[{'time':slice(int(di)-1,int(df)-1)}]=1


                if(tmparea.sea.size > 1):
                    for yr in tmparea.years.values:
                        di=(tmparea.sel(lon=xx,lat=yy,nstats='FirstDay',sea=1,years=yr).item()+\
                            (yr-tmparea.years[0].values)*ndays)
                        df=di+tmparea.sel(lon=xx,lat=yy,nstats='Length',sea=1,years=yr).item()
                        if(np.isfinite(di)): tmp.loc[{'time':slice(int(di)-1,int(df)-1)}]=2

                trsmsk.loc[{'lon':xx,'lat':yy}]=tmp.values
            
            #Reshape it for yrs and days
            trsmsk=trsmsk.assign_coords(years=("time",trsmsk.time.dt.year.values),
                                        days=("time",trsmsk.time.dt.dayofyear.values)).\
                                        set_index(time=["years","days"]).unstack("time")
            
            #**Fraction of the area in transition period
            frc=100*(trsmsk.sum(dim=('lon','lat'))/mskarea.sum())

            #mask with the intensity of thetransition period in the area
            varm=xr.zeros_like(frc)
            for yr in varm.years.values:
                di=onsarea.sel(years=yr,nstats='FirstDay',sea=0).item()
                if np.isnan(di): 
                    df=di
                    varm.loc[{'years':yr}]=np.nan
                else:
                    df=di+onsarea.sel(years=yr,nstats='Length',sea=0).item()
                    varm.loc[{'years':yr,'days':slice(int(di)-1,int(df)-1)}]=\
                        onsarea.sel(years=yr,nstats='OnsIntensity',sea=0)
            varm=varm.where(varm > 0,np.nan)

            #***Average filtered precipitation
            avepp=(bwpp.sel(lat=slice(ymax,ymin),lon=slice(xmin,xmax))*mskarea).mean(dim=['lat','lon'],skipna=True)
            clm=avepp.groupby(avepp.time.dt.dayofyear).mean('time')

            outf=pltdir+'Graphs_FrcAreaFltPP_'+mdl+'_'+area+'_'+aa
            myplt.GraphOnsetPeriods(frc.values,avearea.sea.size,onsarea.years[0].item(),\
                                      avearea.sel(quantile=0.5,nstats=['FirstDay','Length','Onset']).values,\
                                      clm,(0,100,5,10),'Blues',mdl+' - '+aa,'Frc Area (%)',outf,\
                                      b19=onsarea.sel(nstats='Onset',sea=0).values)

            outf=pltdir+'Graphs_OnsIntFltPP_'+mdl+'_'+area+'_'+aa
            myplt.GraphOnsetPeriods(varm.values*10,avearea.sea.size,onsarea.years.values[0],\
                                      avearea.sel(quantile=0.5,nstats=['FirstDay','Length','Onset']).values,\
                                      clm,(0,1.6,.2,.4),'BuPu',mdl+' - '+aa,'Ons Intensity $(x10^{-1} mm.day^{-2})$',outf,\
                                      b19=onsarea.sel(nstats='Onset',sea=0).values)
    ons.close()
    bwpp.close()

  
#Saving the area averaged data
outfile=maindir+'OnsPeriod.V2_AreaAve_'+area+'.npy'
np.save(outfile,onsets,allow_pickle=True)

#boxplot with quantiles over all areas considering selected models
match area:
    case 'SAm':

        mdllst={'All':{'obs':['CHIRPS-UKMO'],
                        'mdl':['UKMO-control','UKMO-rcp85','UKMO-n512','UKMO-n512rcp85',\
                               'h216-hist','h216-ssp585','hg96-hist','hg96-ssp585'],
                        'mcol':['#66CCEE','#66CCEE',"#4477AA","#4477AA",
                                "#EE6677","#EE6677","#AA3377","#AA3377"],
                        'mname':['CHIRPS','CPM-PD','CPM-FUT','GCM-PD','GCM-FUT',\
                                 'HadGEM3-n216-PD','HadGEM3-n216-FUT','HadGEM3-n96-PD',\
                                 'HadGEM3-n96-FUT']},
            }

    case aa if aa in ['Afr','SAfr']:
         mdllst={'All':{'obs':['CHIRPS-CP4'],
                        'mdl':['CP4-CP4','CP4-CP4.fc','CP4-GA25','CP4-GA25.fc',\
                               'h216-hist','h216-ssp585','hg96-hist','hg96-ssp585'],
                        'mcol':['#66CCEE','#66CCEE',"#4477AA","#4477AA",
                                "#EE6677","#EE6677","#AA3377","#AA3377"],
                        'mname':['CHIRPS','CPM-PD','CPM-FUT','GCM-PD','GCM-FUT',\
                                 'HadGEM3-n216-PD','HadGEM3-n216-FUT','HadGEM3-n96-PD',\
                                 'HadGEM3-n96-FUT']},
            }
    case _:
        print('Area not defined for boxplot')
        sys.exit()

vname={'FirstDay':'First Day','Length':'Length','OnsIntensity':'Onset Intensity',\
       'Onset':'Onset Date','PrecipRate':'Precipitation Rate','FrcDryDays':'Fraction Dry Days',
       'TotPrecip':'Total Precipitation','PrecipInt':'Precipitation Intensity',\
       'FrcWetDays':'Fraction Wet Days','LengthWetSpells':'Wet Spells Length',\
       'MaxWetSpell':'Max Wet Spells','LengthDrySpells':'Dry Spells Length',\
       'MaxDrySpell':'Max Dry Spells'}

for scn in mdllst:

    obs=mdllst[scn]['obs']
    mdl=mdllst[scn]['mdl']
    mcol=mdllst[scn]['mcol']

    tmp={}
    selareas=['WAmz','EAmz','EBr'] if area == 'SAm' else ['ZBS','LPP','SAf']
    for dset in obs+mdl:
        tmp[dset]=xr.concat([onsets[dset][aa] for aa in selareas],\
                             dim=xr.DataArray(selareas,dims='area',name='area'))

    ave=xr.concat([tmp[dset] for dset in tmp],dim=xr.DataArray(obs+mdl,dims='mdl',name='mdl'))

    del tmp

    #renaming models
    ave=ave.assign_coords(mdl=mdllst[scn]['mname'])
    obs=[ave.mdl[0].item()]
    mdl=list(ave.mdl[1:].values)

    for vv in ave.nstats.values:
        tmp=ave.sel(nstats=vv,sea=0)
        shd,shdave,unit,ctl=myplt.MapParam(vv,365)

        if(vv in ['FirstDay','Onset']):
            ##Need to transfrom 360 day values to 365 days
            slc=[xx for xx in mdl if not 'bam' in xx]
            tmp.loc[{'mdl':slc}]=tmp.sel(mdl=slc)*(365/360)

            tmp=tmp.where(tmp > 100,tmp+365)
        
        if(vv == 'OnsIntensity'): tmp=tmp*10

        outf=pltdir+'Graphs_CompareAreas_'+scn+'_'+area+'_'+vv
        myplt.LineBox(tmp,(obs,mdl),'area',vname[vv]+' '+unit,vv,mcol,outf)
 
    #Graphic with bars marking the median DWTP in each model
    outf=pltdir+'Graphs_FDLenOnset_'+scn+'_'+area

    indat=ave.sel(sea=0,nstats=['FirstDay','Onset','Length'])
    dset=(obs[0],mdl)
    myplt.FDLengthOns(indat,dset,np.repeat("#899eaa",indat.area.size),\
                      mcol,outf)


    #Saving main variables as a table with format "median [25th;75th]"
    outfile=maindir+area+'_MedianCharacteristic_'+scn+'.csv'
    nvars=['FirstDay','Length','Onset','OnsIntensity','PrecipRate','FrcDryDays']
    ave.loc[{'nstats':'OnsIntensity'}]=ave.sel(nstats='OnsIntensity')*100
    with open(outfile,mode='w') as out_file:
        outf=csv.writer(out_file,delimiter=',',quoting=csv.QUOTE_NONE)
        for aa in selareas:
            row=[aa]+nvars
            outf.writerow(row)

            rnd=[0,1,0,1,1,1]
            for ii in obs+mdl:
                mdn=[str(np.round(ave.sel(mdl=ii,sea=0,nstats=nvars[jj],area=aa,quantile=[0.5]).item(),rnd[jj])) for jj in range(6)]
                lwr=[str(np.round(ave.sel(mdl=ii,sea=0,nstats=nvars[jj],area=aa,quantile=[0.25]).item(),rnd[jj])) for jj in range(6)]
                upr=[str(np.round(ave.sel(mdl=ii,sea=0,nstats=nvars[jj],area=aa,quantile=[0.75]).item(),rnd[jj])) for jj in range(6)]
                mrged=list(zip(mdn,lwr,upr))
                row=[ii]+[f"{a} [{b}; {c}]" for a, b, c in mrged]
                outf.writerow(row)
