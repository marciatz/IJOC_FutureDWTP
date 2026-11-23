#=================================================================================
#This code tests the sensitivity of the climatological envelope to
#different lengths of the butterworth filter
#08/10/2025 - Marcia Zilli
#=================================================================================
import sys, glob
import numpy as np
import DWTP_auxiliary as optb
import xarray as xr
import concurrent.futures
import multiprocessing
from timeit import default_timer as timer

def PxlProcess(xy,pp):

    clmEnv={}
    for ff in bwfilt: 
        minlen=30 #minimum length of the onset period envelope
        deltat=pp.time.dt.dayofyear[0].item()  #==1 for most datasets except CP4, RA25 that starts in March
        fwin=ff
        vyr=np.unique(pp.time.dt.year.values)

        flt=optb.butterworth(pp,fwin)  #filtering

        grd=flt.differentiate('time',datetime_unit='D')
        rny=xr.where(((flt > 0) & (grd > 0)),1,0)
        
        #Climatological envelope
        fd,lns=optb.FindOnsetPeriods(rny.values,minlen) 
        fd=[ii+deltat-1 for ii in fd] 
        #shifting to count from 1 Jan since some datasets start in March
        #final values are in days from Jan 1st YR0, starting from 0.

        if(len(fd) > 0):

            frst,lnt,w0=optb.FindNrRainy(fd,lns,len(vyr),ndays,mxnr=mxrainy)
            #Check the intensity of the precipitation in each envelope if more than one transition
            if(len(frst) > 1):
                tmp=np.zeros((len(vyr)*ndays)); tmp[:]=np.nan
                tmp[deltat-1:len(flt)+deltat-1]=flt 
                tmp=np.nanmean(np.reshape(tmp,(len(vyr),ndays)),0)
                fave=np.asarray([np.nanmean(tmp[frst[ss]:frst[ss]+lnt[ss]]) \
                                    for ss in range(len(frst))])
                pst=[ii for ii in range(len(frst)) if fave[ii] >= (np.max(fave)*.35)] 
                #keep envelopes that have at least 20% of the maximum mean precipitation
                if(len(pst) > 0):
                    frst=frst[pst]; lnt=lnt[pst]
            if(len(frst) == 0):
                frst,lnt,w0=np.nan,np.nan,np.nan

        else:
            print(xy,'no onset periods found')
            frst,lnt,w0=np.nan

        clmEnv[ff]=[frst+1,lnt,w0+1]
    return xy,clmEnv
        

def workerXY(wrk):
    class_data=[]
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(PxlProcess,pxl[xy],precip.sel(lon=pxl[xy][1],\
                                              lat=pxl[xy][0])) for xy in xylst[wrk]]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]

    class_data.append([wrk,results])
    
    del results
    print(wrk,'worker done')

    return class_data

#===================================================================================
#Reading daily precip data
dset=sys.argv[1]  #Br-DWGD ERA5 CHIRPS bam hg96 h216, UKMO (SAM), SAAG (SAM), CP4 (SAfr)
prd=sys.argv[2]   #Simulation: 'amip' 'hist' 'ssp585' for UKMO-GCM models
                  #            'm1' 'm2' 'm3' 'm4' for BAM; all for ensemble mean
                  #            'control' 'n512' 'hindcast' 'ERAI25km' 'rcp85' 'n512rcp85' for UKMO-CPM
                  #            '' for SAAG-WRF CPM
                  #            'CP4', 'RA25', 'GA25', 'CP4.fc', 'RA25.fc', 'GA25.fc'  for CP4 SAfr CPM
                  #            '' for observations or UKMO, SAAG, or CP4 if comparing to CPM
area=sys.argv[3]

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

#Defining HPC used
maindir='/gws/ssde/j25a/firstrains/OnsetPeriod/data/'


yrs={'IMERG':(1998,2025),'CHIRPS':(1981,2025), 
     'bam':(1981,2010),'hg96':(1979,2014),'h216':(1979,2014),
     'UKMO':(1998,2007),'SAAG':(2000,2021),'CP4':(1997,2006),
     } 
yr0,yrf=yrs[dset]; nyrs=yrf-yr0+1

nstats=['FirstDay','Length','OnsIntensity','TotPrecip','PrecipRate',\
        'PrecipInt','FrcWetDays','LengthWetSpells','MaxWetSpell',\
        'FrcDryDays','LengthDrySpells','MaxDrySpell']
 
#Reading daily precipitation values
match dset:
    case 'CHIRPS':
        indir=maindir+dset+'/'+dset+'_'
        infile=indir+'precip.day.total_1981-2025.025dg.aw.nc' 

        if area=='global':
            infile=indir+'precip.day.total_1981-2025.05dg.bl.nc'

        if prd == '0':   #using all data available
            with xr.open_dataset(infile) as inp:
                if(area == 'Afr'):
                    inp=inp.assign_coords(longitude=(((inp.longitude+180)%360)-180)).sortby('longitude')
                precip=inp.total_precipitation.sel(longitude=slice(xxi,xxf),latitude=slice(yyi,yyf))
            del inp
        else: #using same period as pr dataset (CPM comparison)
            yr0,yrf=yrs[prd]; nyrs=yrf-yr0+1
            if(prd == 'UKMO'):
                xxi,xxf=(276.5,328.5); yyi,yyf=(13.5,-38.5)  #UKMO CPM area over SAm
            elif(prd == 'SAAG'):
                xxi,xxf=(261.,343.0); yyi,yyf=(16,-58.)      #SAAG CPM area over SAm
            else:
                xxi,xxf=(-24.5,56.5); yyi,yyf=(39.75,-45.75)      #CP4 CPM area over Afr

            with xr.open_dataset(infile) as inp:
                if(area == 'Afr'):
                    inp=inp.assign_coords(longitude=(((inp.longitude+180)%360)-180)).sortby('longitude')

                precip=inp.total_precipitation.sel(longitude=slice(xxi,xxf),\
                                                   latitude=slice(yyi,yyf),\
                                                time=inp.time.dt.year.isin(range(yr0,yrf+1)))

            del inp
        
        #rename 'longitude' and 'latitude' to 'lon' and 'lat'
        precip=precip.rename({'longitude':'lon','latitude':'lat'})

        #Mask for CHIRPS
        mskfile=maindir+'mask_CHIRPS_0.25dg.nc'
        with xr.open_dataarray(mskfile) as inp:
            if(area == 'Afr'):
                inp=inp.assign_coords(longitude=(((inp.longitude+180)%360)-180)).sortby('longitude')
            msk=inp.sel(longitude=precip.lon.values,latitude=precip.lat.values).isel(time=0)
        del inp

        msk=msk.rename({'longitude':'lon','latitude':'lat'})


    case 'IMERG':
        indir=maindir+dset+'/'+dset+'_'    
        infile='/work/scratch-nopw2/franmorr/gpm_daily/imerg_1998-2025_daily.nc4'

        if prd == '0':   #using all data available within the configured IMERG span
            with xr.open_dataset(infile) as inp:
                if(area == 'Afr'):
                    inp=inp.assign_coords(lon=(((inp.lon+180)%360)-180)).sortby('lon')
                precip=inp.precipitation.sortby('lat',ascending=False).sel(
                    lon=slice(xxi,xxf),
                    lat=slice(yyi,yyf),
                    time=slice(f'{yr0}-01-01',f'{yrf}-12-31')
                )
            del inp
        else: #using same period as pr dataset (CPM comparison)
            yr0,yrf=yrs[prd]; nyrs=yrf-yr0+1
            if(prd == 'UKMO'):
                xxi,xxf=(276.5,328.5); yyi,yyf=(13.5,-38.5)  #UKMO CPM area over SAm
            elif(prd == 'SAAG'):
                xxi,xxf=(261.,343.0); yyi,yyf=(16,-58.)      #SAAG CPM area over SAm
            else:
                xxi,xxf=(-24.5,56.5); yyi,yyf=(39.75,-45.75)      #CP4 CPM area over Afr

            with xr.open_dataset(infile) as inp:
                if(area == 'Afr'):
                    inp=inp.assign_coords(lon=(((inp.lon+180)%360)-180)).sortby('lon')

                precip=inp.precipitation.sel(
                    lon=slice(xxi,xxf),
                    lat=slice(yyi,yyf),
                    time=slice(f'{yr0}-01-01',f'{yrf}-12-31')
                )

            del inp

        msk=xr.where(precip.notnull().sum('time') > 0,1,np.nan)
        yr0=precip.time.dt.year[0].item()
        yrf=precip.time.dt.year[-1].item(); nyrs=yrf-yr0+1
        
    case 'bam':
        indir=maindir+dset+'/'
        if prd == 'all':   #using ensemble mean
            filelst=glob.glob(indir+'*/GPOS1000PREC19810101122010012012C.fct.TQ0126L042*')
            with xr.open_mfdataset(filelst,combine='nested',concat_dim='ens') as inp:
                if(area == 'Afr'):
                    inp=inp.assign_coords(lon=(((inp.lon+180)%360)-180)).sortby('lon')
                precip=inp.prec.sel(lon=slice(xxi,xxf),lat=slice(yyf,yyi)).mean(dim='ens').compute()

        else:
            indir=indir+prd+'/'
            infile=indir+'GPOS1000PREC19810101122010012012C.fct.TQ0126L042_'+prd+'.nc'

            with xr.open_dataset(infile) as inp:
                if(area == 'Afr'):
                    inp=inp.assign_coords(lon=(((inp.lon+180)%360)-180)).sortby('lon')
                precip=inp.prec.sel(lon=slice(xxi,xxf),lat=slice(yyf,yyi))
            del inp

        #Creating its onw mask, considering sea areas
        msk=xr.where(precip.notnull().sum('time') > 0,1,np.nan)

    case d if d in ['h216','hg96']:
        indir=maindir+dset+'/'+prd+'/'

        if (prd == 'ssp585'): yr0=2070; yrf=2099; nyrs=yrf-yr0+1
        infile=indir+"hadgem3_n96_"+prd+".r1.precip.daily.CMIP6.nc" \
                if (dset == 'hg96') else indir+"hadgem3_n216_"+prd+".r1.precip.daily.CMIP6.nc"

        with xr.open_dataset(infile) as inp:
            if(area == 'Afr'):
                inp=inp.assign_coords(longitude=(((inp.longitude+180)%360)-180)).sortby('longitude')
            precip=inp.pr.sel(longitude=slice(xxi,xxf),latitude=slice(yyf,yyi),\
                              time=inp.pr.time.dt.year.isin(range(yr0,yrf+1))) \
                   if (prd in ['amip','hist']) else \
                   inp.precipitation_flux.sel(longitude=slice(xxi,xxf),latitude=slice(yyf,yyi),\
                              time=inp.precipitation_flux.time.dt.year.isin(range(yr0,yrf+1)))
        del inp

        precip=precip*86400 #converting from kg.m-2.s-1 to mm.day-1
        precip=precip.rename({'longitude':'lon','latitude':'lat'})
        
        #Creating its onw mask, considering sea areas
        msk=xr.where(precip.notnull().sum('time') > 0,1,np.nan)

    case 'UKMO':
        indir=maindir+dset+'/'+prd+'/'+prd+'_'

        infile=indir+'precip.day.total_'+str(yr0)+'-'+str(yrf)+'.025dg.aw.nc'

        with xr.load_dataset(infile) as inp:
            if(prd == 'n512rcp85'):
                xxi,xxf=(276.5,328.5); yyi,yyf=(13.5,-38.5)
                precip=inp.total_precipitation.sel(longitude=slice(xxi,xxf),latitude=slice(yyi,yyf))
            else:
                precip=inp.total_precipitation
        del inp

        #rename 'longitude' and 'latitude' to 'lon' and 'lat'
        precip=precip.rename({'longitude':'lon','latitude':'lat'})

        #Mask considering CHIRPS land-only area at era5 resloution
        mskfile=maindir+'mask_CHIRPS_0.25dg.nc'
        with xr.open_dataarray(mskfile) as inp:
            msk=inp.sel(longitude=precip.lon.values,latitude=precip.lat.values).isel(time=0)
        del inp

        msk=msk.rename({'longitude':'lon','latitude':'lat'})

    case 'SAAG':
        indir=maindir+dset+'/'+dset.lower()+'_'

        infile=indir+'precip.day.total_'+str(yr0)+'-'+str(yrf)+'.025dg.aw.nc'
        # if(dset == 'SAAG-WRF'): calendar='standard'

        with xr.open_dataset(infile) as inp:
            precip=inp.total_precipitation
        del inp
        #rename 'longitude' and 'latitude' to 'lon' and 'lat'
        precip=precip.rename({'longitude':'lon','latitude':'lat'})
        precip=xr.where(precip > 10e4,np.nan,precip)

        #Mask considering CHIRPS land-only area at era5 resloution
        mskfile=maindir+'mask_CHIRPS_0.25dg.nc'
        with xr.open_dataarray(mskfile) as inp:
            msk=inp.sel(longitude=precip.lon.values,latitude=precip.lat.values).isel(time=0)
        del inp

        msk=msk.rename({'longitude':'lon','latitude':'lat'})
        
    case 'CP4':
        indir=maindir+dset+'/'+prd+'/'+prd+'.'
        infile=indir+str(yr0)+'-'+str(yrf)+'.'+area+'.ERA5.nc'

        with xr.open_dataset(infile) as inp:
            if 'GA25' in prd:
                inp=inp.assign_coords(longitude=(((inp.longitude+180)%360)-180)).sortby('longitude')
                xxi,xxf=(-25.5,56.5); yyi,yyf=(39.75,-45.75)
                precip=inp.precipitation_flux.sel(longitude=slice(xxi,xxf),\
                                                    latitude=slice(yyi,yyf))
            else:
                precip=inp.pr

        precip=precip.rename({'longitude':'lon','latitude':'lat'})

        #Mask considering CHIRPS land-only area at era5 resloution
        #in CP4, the longitude writen as -180 to 180 because it bridges over the zero meridian
        with xr.open_dataarray(maindir+'mask_CHIRPS_0.25dg.nc') as inp:
            inp=inp.assign_coords(longitude=(((inp.longitude+180)%360)-180)).sortby('longitude')
            msk=inp.sel(longitude=precip.lon.values,latitude=precip.lat.values).isel(time=0)
        del inp

        msk=msk.rename({'longitude':'lon','latitude':'lat'})
    
        yr0=precip.time.dt.year[0].item()
        yrf=precip.time.dt.year[-1].item(); nyrs=yrf-yr0+1

    case _: 
        print('Dataset not found')

#Remove leap day
if precip.time.dt.calendar in ['noleap','360_day']:
    ndays=360 if precip.time.dt.calendar == '360_day' else 365
else:
    ndays=365
    if np.any((precip.time.dt.month==2) & (precip.time.dt.day ==29)):
        leap=precip.where((precip.time.dt.month == 2) &\
                            (precip.time.dt.day == 29),drop=True).coords['time'].values
        precip=precip.drop_sel(time=leap)
    
    #Converting the calendar enconding from datetime64 to cftime
    time=xr.date_range(start=str(yr0),periods=precip.time.size,\
                       freq='D',calendar='noleap',use_cftime=True)
    precip=precip.assign_coords(time=time)

mxrainy=2 #maximum number of climatological envelopes
tstname='BWFlt'; bwfilt=[30,40,50,60,70,80,90]

climEnv=xr.DataArray(coords=(precip.lat.values,precip.lon.values,\
                             list(range(mxrainy)),bwfilt,['FirstDay','Length']),\
                        dims=('lat','lon','sea',tstname,'Stats'))

tmp=(precip.isnull().sum(dim='time') > 10)

msk=xr.where(tmp,np.nan,msk)
tmp=xr.where(np.isnan(msk),0,1).stack(z=('lat','lon'))
pxl=tmp.where(tmp == 1,drop=True).coords['z'].values 

print(dset)

#Distributing tasks
mxwkrs=60
ntasks=int(len(pxl)/mxwkrs)
xylst=[[ii for ii in range(wrk*ntasks,wrk*ntasks+ntasks)] for wrk in range(mxwkrs)]
xtr=len(pxl) % mxwkrs
if(xtr > 0): 
    xtasks=len(pxl)-mxwkrs*ntasks
    for ii in range(xtasks): 
        xylst[ii % mxwkrs].append(mxwkrs*ntasks+ii)

t0=timer()
with multiprocessing.Pool(processes=mxwkrs) as pool:
    class_data = pool.map(workerXY,list(range(mxwkrs)))

t1=timer()
print('Refining onset periods: ',t1-t0)

#Retrieving results
for ii in range(mxwkrs):
    wrk=class_data[ii][0][0]
    tmp=class_data[ii][0][1]
    for jj in range(len(tmp)):
        yy,xx=tmp[jj][0]

        clm=tmp[jj][1]
        for kk in clm.keys():
            if(np.all(np.isnan(clm[kk][0]))):
                climEnv.loc[{'lat':yy,'lon':xx,'BWFlt':kk}]=np.nan
            else:
                climEnv.loc[{'lat':yy,'lon':xx,'sea':slice(0,len(clm[kk][0])),'BWFlt':kk,'Stats':'FirstDay'}]=clm[kk][0]
                climEnv.loc[{'lat':yy,'lon':xx,'sea':slice(0,len(clm[kk][1])),'BWFlt':kk,'Stats':'Length'}]=clm[kk][1]


#######################
#Saving the results


outdir=maindir+dset+'/'
if dset=='IMERG':
    outdir='/gws/ssde/j25a/firstrains/imerg_onset/'+dset

aa='' if area=='global' else area

if dset in ['CHIRPS','IMERG']: 
    if prd == '0':
        envfile=outdir+aa+"_OnsPerClimEnvelope_"+tstname+"_025dg.aw.nc"
    else:
        envfile=outdir+aa+"_"+prd+"_OnsPerClimEnvelope_"+tstname+"_025dg.aw.nc"
else: 
    if ((dset == 'bam') & (prd=='all')):
        pass
    else:
        outdir=outdir+prd+'/'

    if dset in ['UKMO','SAAG']:
        envfile=outdir+aa+"_OnsPerClimEnvelope_"+tstname+"_025dg.aw.nc"
    else:
        envfile=outdir+aa+"_OnsPerClimEnvelope_"+tstname+".nc"

#Climatological Envelope
clm=xr.Dataset({'clmEnv':climEnv},coords={'lat':climEnv.lat,'lon':climEnv.lon,\
                                          'sea':climEnv.sea,'BWFlt':climEnv.BWFlt,'Stats':climEnv.Stats})
print("Saving",envfile)
clm.to_netcdf(envfile,'w',encoding={'clmEnv':{"zlib":True}})

del clm