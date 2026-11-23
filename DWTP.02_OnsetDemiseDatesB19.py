#=========================================================================================
#This code calculates the onset and demise of the rainy season and is heavely based
#on the methodology described in Bombardi et al 2019 (doi: 10.1175/BAMS-D-18-0177.1)
#and related codes (https://github.com/rjbombardi/onset_demise_rainy_season/tree/master)
#22/05/2023 - Marcia Zilli
#25/06/2024 - Updated the code to improve is performance. Marcia Zilli
#=========================================================================================
import sys,glob
import numpy as np
import xarray as xr
import DWTP_auxiliary as optb
import concurrent.futures
import multiprocessing
import cftime
from timeit import default_timer as timer
from datetime import timedelta as dt

def PxlProcess(xy,pp):

    yy,xx=xy

    #Creating the harmonic variables    
    cycle,w0=optb.AnnualCycleB19(pp)

    rm=pp.mean(dim='time')

    #The accumulated precipitation is accounted from the onset to the demise from 
    #that water year; the dry season precipitation is accumulated from the demise
    #of the rainy season in the previous water year to the onset in this water year. 
    #That way, the order of the seasons are dry and then wet, even though the 
    #dry season starts before the beginning of the water year and finishes after. 
    yrs=np.unique(pp.time.dt.year.values)

    ons=xr.DataArray(np.nan, dims=['year'], coords={'year': yrs})
    dem=xr.full_like(ons,np.nan)
    length=xr.full_like(ons,np.nan)
    wet=xr.full_like(ons,np.nan)
    dry=xr.full_like(ons,np.nan)

    if(np.isnan(w0)):  #Rainy seasons not well defined, skipping grid point
        pass
    else:
        vyr=pp.sel(time=pp.time.dt.dayofyear==w0).time
        #Considering the year as the w0 year. If w0 is on day 330 (1 Dec) and the onset 
        #and demise are day 20 and 200 of the following year, the onset and demise are
        #assigned to the year of day 330, i.e., the water year.

        ppanom=pp-rm
        fd,dm=optb.SeasonOnsetDemiseB19(ppanom,w0)

        for ii in range(len(vyr)):
            #Locating the onset/demise closest to the begining of the water year (w0)
            dst=np.zeros(len(fd))
            for kk,jj in enumerate(fd):
                dst[kk]=(jj.item()-vyr[ii].item()).days if jj.notnull() else np.nan
            dst=np.where(dst > 0,dst,np.nan)  #Considering only the onsets that are after the start of the water year
            if  (dst < ndays).any(): #
                #Onset is consiredered if there is a candidate within one year of the beginning of the water year; otherwise the year is set as missing
        
                ps0=fd[np.nanargmin(dst)].item() 
                
                #onset of the rainy season for the water year ii
                ons.loc[{'year':vyr[ii].time.dt.year}]=ps0.dayofyr

            else:
                ps0=np.nan

            #Similar for demise
            dst=np.zeros(len(dm))
            for kk,jj in enumerate(dm):
                dst[kk]=(jj.item()-vyr[ii].item()).days if jj.notnull() else np.nan
            dst=np.where(dst > 0,dst,np.nan)  
            #Considering only the demises that are after the start of the water year
            if (dst < ndays).any():
                ps1=dm[np.nanargmin(dst)].item()

                if(ons.sel(year=vyr[ii].time.dt.year).notnull()):
                    if((ps1-ons.sel(year=vyr[ii].time.dt.year).time.item()).days > 0):
                        dem.loc[{'year':vyr[ii].time.dt.year}]=ps1.dayofyr
                        #demise of the rainy season for the water year ii
                else:
                    ps1=pp.isel(time=-1).time.item() if ii == len(vyr)-1 else xr.DataArray(np.nan)

            else:
                ps1=np.nan

            if ((isinstance(ps0,cftime.datetime)) & (isinstance(ps1,cftime.datetime))):
                length.loc[{'year':vyr[ii].time.dt.year}]=(ps1-ps0).days
                #total precipitation during the rainy season starting after w0
                wet.loc[{'year':vyr[ii].time.dt.year}]=\
                    pp.sel(time=slice(ps0,ps1)).sum().item()
                dry.loc[{'year':vyr[ii].time.dt.year}]=\
                    pp.sel(time=slice(vyr[ii],vyr[ii]+dt(days=ndays-1))).sum().item()-\
                    wet.sel(year=vyr[ii].time.dt.year)
            else:
                length.loc[{'year':vyr[ii].time.dt.year}]=np.nan
                wet.loc[{'year':vyr[ii].time.dt.year}]=np.nan
                dry.loc[{'year':vyr[ii].time.dt.year}]=np.nan


    return xy,(ons,dem,length,wet,dry)


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

#===================

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
        xxi,xxf,yyi,yyf=[0,360,90,-90]
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

yrs={'CHIRPS':(1981,2025),
     'bam':(1981,2010),'hg96':(1979,2014),'h216':(1979,2014),
     'UKMO':(1998,2007),'SAAG':(2000,2021),'CP4':(1997,2006)} 
yr0,yrf=yrs[dset]; nyrs=yrf-yr0+1

#Reading precipitation data and mask (masking ocean)
match dset:
    case 'CHIRPS':
        indir=maindir+dset+'/'+dset+'_'
        infile=indir+'precip.day.total_1981-2025.025dg.aw.nc'
        if area=='global':
            infile=indir+'precip.day.total_1981-2025.05dg.bl.nc'

        if prd == '0':   #using all data available
            with xr.load_dataset(infile) as inp:
                if(area == 'Afr'):
                    inp=inp.assign_coords(longitude=(((inp.longitude+180)%360)-180)).sortby('longitude')
                if(area == 'global'):
                    precip=inp.precip  #Read data with variable now called precip 
                    precip=precip.assign_coords(lon=(inp.lon % 360)).sortby("lon") #Fix longitude from -180W-180E to 0-360
                    precip=precip.reindex(lat=list(reversed(precip.lat))) #flipping the latitude so that it goes from N to S, as before
                    precip=precip.rename({'lon':'longitude','lat':'latitude'})  #Changing coord names to latitude and longitude
                else:
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

            with xr.load_dataset(infile) as inp:
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
        infile=indir+"hadgem3_n96_"+prd+".r1.precip.daily.CMIP6.nc" if (dset == 'hg96') else \
            indir+"hadgem3_n216_"+prd+".r1.precip.daily.CMIP6.nc"

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
    cal='360_day' if precip.time.dt.calendar == '360_day' else 'noleap'
    hrtime=xr.date_range(start=str(yr0)+'-'+str(precip.time.dt.month[0].item()).zfill(2)+'-01',\
                        periods=precip.time.size,freq='D',calendar=cal,use_cftime=True)
else:
    ndays=365
    if np.any((precip.time.dt.month==2) & (precip.time.dt.day ==29)):
        leap=precip.where((precip.time.dt.month == 2) &\
                            (precip.time.dt.day == 29),drop=True).coords['time'].values
        precip=precip.drop_sel(time=leap)
        precip=precip.convert_calendar('noleap') #convert to cftime noleap calendar

#Locating pixels with valid data (land only)
msk=xr.where(xr.where(precip.isnull(),0,1).sum(dim='time') >= .5*precip.time.size,1,np.nan) 
tmp=msk.stack(z=('lat','lon'))
pxl=tmp.where(tmp.notnull(),drop=True).coords['z'].values


#Parallelising the processing of the grid points
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
print('Calculating onset and demise dates: ',t1-t0)

#Retreiveing outputs    
onset=xr.DataArray(coords={'year':np.unique(precip.time.dt.year.values),\
                           'latitude':precip.lat.values,\
                           'longitude':precip.lon.values},\
                    dims=['year','latitude','longitude'])
demise=xr.full_like(onset,np.nan)
duration=xr.full_like(onset,np.nan)
wet=xr.full_like(onset,np.nan)
dry=xr.full_like(onset,np.nan)

#Retrieving results
for ii in range(mxwkrs):
    wrk=class_data[ii][0][0]
    tmp=class_data[ii][0][1]
    for jj in range(len(tmp)):
        yy,xx=tmp[jj][0]
        onset.loc[{'latitude':yy,'longitude':xx}]=tmp[jj][1][0]
        demise.loc[{'latitude':yy,'longitude':xx}]=tmp[jj][1][1]
        duration.loc[{'latitude':yy,'longitude':xx}]=tmp[jj][1][2]
        wet.loc[{'latitude':yy,'longitude':xx}]=tmp[jj][1][3]
        dry.loc[{'latitude':yy,'longitude':xx}]=tmp[jj][1][4]

#========================================================================
# Saving rainy and wet season characteristics
#========================================================================


outdir=maindir+dset+'/'

aa='' if area=='global' else area

if dset == 'CHIRPS': 
    if prd == '0':
        outfile=outdir+aa+"_OnsetDemise_Bombardi2019.nc" 
    else:
        outfile=outdir+aa+"_"+prd+"_OnsetDemise_Bombardi2019_025dg.aw.nc"
else: 
    if ((dset == 'bam') & (prd=='all')):
        pass
    else:
        outdir=outdir+prd+'/'

    if dset in ['UKMO','SAAG']:
        outfile=outdir+aa+"_OnsetDemise_Bombardi2019_025dg.aw.nc"
    else:
        outfile=outdir+aa+"_OnsetDemise_Bombardi2019.nc"

print('Saving results at ',outfile)


onsTime=xr.Dataset(data_vars={'onset':onset,'demise':demise,'duration':duration,\
                              'wet':wet,'dry':dry},\
                   attrs={'description':"Onset and demise of the rainy season, duration of the wet season, and total precipitation during the wet and dry seasons. Note that the onset and demise dates are in day of the year, but they might not have occurred in the same calendar year. For example, if the onset is on day 330 (1 Dec) and the demise is on day 20 of the following year, the onset and demise are assigned to the year of day 330, i.e., the water year."})
onsTime.to_netcdf(outfile,
                  encoding={"onset":{'zlib':True,'complevel':4},
                            "demise":{'zlib':True,'complevel':4},
                            "duration":{'zlib':True,'complevel':4},
                            "wet":{'zlib':True,'complevel':4},
                            "dry":{'zlib':True,'complevel':4}})
