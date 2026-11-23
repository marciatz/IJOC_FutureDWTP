import sys
import numpy as np
import xarray as xr
import itertools as it
from datetime import timedelta as dt
import scipy.fft as fft

#############################
#Main functions related to Onset period analysys
##############################

def FindOnsetPeriods(indat,minlen):
    '''
        Given an time series of 0's and 1's, this functions locates the 
        position and length of sequences of 1's with minimum number of entries
        defined by minlen. 
        
        INPUT:
            indat - time series of 0's and 1's
            minlen - minimum number of days in a sequence of 1's
        
        RETURNS:
            firstday - list with the position of the first entry in a 
                       sequence of 1's longer than minlen
            lenseason - number of conseqcutive entries equal 1.
    '''

    ps=0; lenseason=[]; firstday=[]
    for key,grp in it.groupby(indat):
        tmp=len(list(grp))
        if key == 0:
            pass
        else:
            if(tmp >= minlen):
                lenseason.append(tmp)
                firstday.append(ps)  
        ps=ps+tmp              

    return firstday,lenseason

def FindNrRainy(fd,ln,nyrs,ndays,minlen=25,mxnr=False,tmpmask=False,shr=0.5):
    '''
        Identifies the number of rainy seasons and the expected 
        window that it should occur. The identification is based on
        the frequency (over the years) with each each day is classified
        as having a sustainded precipitation increase (positive gradient
        of the filtered precipitation) considering all years in the dataset.

        INPUT:
            fd     - list with the first day of each period with
                    sustained precipitation increase
            ln     - duration of the period 
            nyrs   - number of years in the dataset
            ndays  - number of days per year (UM model has 360 days)
            minlen - minimum number of days in a period for it to
                    be considered a onset period
            mxnr   - maximum number of rainy seasons allowed

        RETURN:
            frst  - list with the first day of each climatological onset period
            lnt   - durantion of each climatological onset period
            w0    - day of the year considered as the beginning of the water year.
                    Defined as the day with the minimum freqeuncy of days with
                    sustained precipitation increase
            pdfNS - the frequency of years with sustained precipitation increase
                    (per day of the year), after nomalisation (values range from 0
                    to 1) and smoothing (15-day running mean)
    '''

    #mask with periods in onset 
    if(isinstance(tmpmask,bool)):
        tmpmask=np.zeros((nyrs*ndays))
        for ii in range(len(fd)):
            ps=list(range(fd[ii],fd[ii]+ln[ii]))
            tmpmask[ps]=1
        
        tmpmask=np.reshape(tmpmask,(nyrs,ndays))


    #Number of times a day is classified as onset over all years
    if (len(tmpmask.shape) > 1): 
        pdf=np.nansum(tmpmask,0)
    else: pdf=tmpmask

    pdfN=(pdf-pdf.min())/(pdf.max()-pdf.min()) 
    #normalizing such that values range from 0 to 1. Note that 1 does not mean that 
    #the days was in an onset period over all years, but that that was the day classified 
    #as onset the largest number of years

    #locating the minimum and shifting the time series to start at this point
    #Not done for seasonal hindcast because it does not cover one entire year
    w0=np.argmin(pdfN)
    if(w0 >0): pdfN=np.hstack((pdfN[w0:],pdfN[:w0]))

    #smoothing using a running mean
    ws=14
    padded = np.pad(pdfN, (ws//2,ws-ws//2-1), mode='edge')
    pdfNS = np.convolve(padded, np.ones((ws,))/ws, mode='valid')

    #Locating periods with largest frequence of onset (above 0.5)
    #(not really 50% since it is normalized by the maximum)
    pps=[ps for ps in range(ndays-1) if ((pdfNS[ps] < shr) & (pdfNS[ps+1]>=shr))]
    nps=[ps for ps in range(ndays-1) if ((pdfNS[ps] >= shr) & (pdfNS[ps+1]<shr))]

    if(len(pps) == 0): 
        if(pdfNS.max() < 0.5):
            pps=[ps for ps in range(ndays-1) if ((pdfNS[ps] < shr*pdfNS.max()) &\
                                                (pdfNS[ps+1]>=shr*pdfNS.max()))]
            nps=[ps for ps in range(ndays-1) if ((pdfNS[ps] >= shr*pdfNS.max()) & \
                                                (pdfNS[ps+1]<shr*pdfNS.max()))]
        if(pdfNS.min() > 0.5):
            thr=shr*(pdfNS.max()-pdfNS.min())+pdfNS.min()
            pps=[ps for ps in range(ndays-1) if ((pdfNS[ps] < thr) &\
                                                 (pdfNS[ps+1]>=thr))]
            nps=[ps for ps in range(ndays-1) if ((pdfNS[ps] >= thr) & \
                                                 (pdfNS[ps+1]<thr))]
    if(len(pps) == 0):
        return [],[],[]

    #Sanity check point: since the time series was shifted to start at the 
    #minimum value, the curve should cross the y=0.5 an even number of 
    #time (one going up and one going dow)
    if(len(nps) != len(pps)):
        if(len(nps) < len(pps)):
            if(pdfNS[-1] >= shr): nps.append(ndays-1)
            else: 
                if(pdfNS[-1] >= shr*pdfNS.max()): nps.append(ndays-1)
                else: print ('!!!!!CHECK!!!!')
        else:
            if(pdfNS[0] >= shr): pps.insert(0,0)
            else: 
                if(pdfNS[0] >= shr*pdfNS.max()): pps.insert(0,0)
                else: print ('!!!!!CHECK!!!!')
    else:
        if(pps[0] > nps[0]): pps.insert(0,0)
        if(pps[-1] > nps[-1]): nps.append(ndays-1)

            
    #checking for breaks in the periods (value below 0.5 for a few days)
    mrg=[]
    for ii in range(len(pps)-1):
        if(pps[ii+1]-nps[ii] <= 30): #Arbitrary, check if it is a good size
            if (np.min(pdfNS[nps[ii]:pps[ii+1]+1]) > 0.3): mrg.append(ii)
            #in this case, the minimum value between the two periods is 
            #too small, so it is likely a break in the season. Merge both seasons
                
    if(len(mrg) >0):
        mrg.sort()
        for ii in mrg[::-1]: nps.pop(ii); pps.pop(ii+1)

    #01/03/25 - lower thrshold for seasons not present all years (pdfNS.max < 1)
    if((len(pps) > 1) | (len(nps) > 1)):
        ps=(np.asarray([pdfNS[pps[ii]:nps[ii]].max() \
                        for ii in range(len(pps))]) < pdfNS.max()).nonzero()[0]
        if(len(ps) > 0):
            for jj in ps:
                tmp1=[ii for ii in range(ndays-1) \
                            if ((pdfNS[ii] < shr*(pdfNS[pps[jj]:nps[jj]].max())) &\
                                (pdfNS[ii+1] >= shr*(pdfNS[pps[jj]:nps[jj]].max())))]
                tmp2=[ii for ii in range(ndays-1) \
                            if ((pdfNS[ii] >= shr*(pdfNS[pps[jj]:nps[jj]].max())) &\
                                (pdfNS[ii+1] < shr*(pdfNS[pps[jj]:nps[jj]].max())))]
                if(len(tmp2) < len(tmp1)):tmp2.append(ndays-1)
                if(len(tmp1) < len(tmp2)):tmp1.insert(0,0)
                if((len(tmp1) > 0) | (len(tmp2) > 0)):
                    if(tmp1[0] > tmp2[0]): tmp1.insert(0,0)
                    if(tmp1[-1] > tmp2[-1]): tmp2.append(ndays-1)

                    #locating the correct original interval - new first day that is closest
                    #to the original first day
                    kk=np.argmin(abs(pps[jj]-np.asarray(tmp1)))
                    tmp1=tmp1[kk]; tmp2=tmp2[kk]
                    #Check if new period is separate from other candidates
                    keep=True
                    for ii in range(len(pps)):
                        if (ii != jj):
                            #Checking if new envelope overlaps with other periods
                            if(tmp1 in range(pps[ii],nps[ii])): keep=False
                            if(tmp2 in range(pps[ii],nps[ii])): keep=False
                            #Check if new envelope contains other periods
                            if(pps[ii] in range(tmp1,tmp2)): keep=False
                            if(nps[ii] in range(tmp1,tmp2)): keep=False
                    if(keep): pps[jj]=tmp1; nps[jj]=tmp2

    #eliminating periods with less than minlen days 
    ln=np.asarray([nps[ii]-pps[ii]+1 for ii in range(len(pps))])
    if(np.all(ln < minlen)):  #In case all periods are too short, check if 
                              #the longest one is at least 75% of the minlen (19 days)
        if(ln.max() > .75*minlen):
            ps=np.argmax(ln)
            pps=[pps[ps]]; nps=[nps[ps]]
    else:
        rmv=(ln < minlen).nonzero()[0]
        if(len(rmv) > 0):
            rmv.sort()
            for ii in rmv[::-1]: pps.pop(ii); nps.pop(ii)

    #If there is a limit to the number of rainy seasons to be considered
    if(not isinstance(mxnr,bool)):  
        while(len(pps) > mxnr):   #For now curbing the nr of rainy seasons to a max
            #Eliminating the seasons with smaller pdfNS
            mx=np.asarray([pdfNS[pps[ii]:nps[ii]].max() for ii in range(len(pps))])
            ii=np.argmin(mx)
            pps.pop(ii); nps.pop(ii)


    #Returning first day of the onset periods, their length, and 
    #the position of the minimum - this is the day with the least 
    #number of onset candidates
    lnt=np.asarray(nps)-np.asarray(pps)
    frst=np.asarray(pps)+w0
    frst=np.where(frst >=ndays,frst-ndays,frst)
    #sorting according to the first day
    lnt=lnt[np.argsort(frst)]
    frst=frst[np.argsort(frst)]
    pdfNS=np.hstack((pdfNS[-w0:],pdfNS[:(ndays-w0)]))

    return frst,lnt,w0

def WetDrySpells(msk,pp,ns,ons):
    '''
        Function to estimate wet and dry spell statistics during onset periods:
            PrecipRate: Total precipitation divided by the length of the onset period
            PrecipInt: Total precipitation divided by the number of wet days in the onset period
            FrcWetDays: Number of wet days divided by the length of the onset period
            LengthWetSpell: number of wet days (PP >= 1) divided by the number of wet spells (sequence of wet days)
            MaxWetSpell: longest wet spell
            FrcDryDays: Number of dry days divided by the length of the onset period
            LengthDrySpell: number of dry days (PP < 1) divided by the number of dry spells (sequence of dry days)
            MaxDrySpell: longest dry spell
        
        INPUT:
            msk - mask with the onset periods
            pp - precipitation time series
            ns - number of onset periods
            ons - onset periods
        
        RETURNS:
            onsStats - array with the statistics of the onset periods
    '''

    ntot=len(pp)
    nyrs=ons.shape[0]

    onsStats=np.zeros((nyrs,ns,8)); onsStats[:]=np.nan

    yrmsk=np.where(msk == 1,1,np.nan)
    vldpp=pp*yrmsk
    
    vldpp[vldpp < 1]=0.

    wet=np.where(vldpp > 0,1.,0)
    dd0=(np.diff(wet)==1).nonzero()[0]+1; 
    ddf=(np.diff(wet)==-1).nonzero()[0]
    #note that ddf does not have the +1 because diff[i]=x[i+1]-x[i].
    if(wet[0]==1):dd0=np.insert(dd0,0,0)
    if(wet[-1]==1):ddf=np.insert(ddf,len(ddf),ntot-1)
    wetlen=[(ii,jj,vldpp[int(ii):int(jj)+1]) for ii,jj in zip(dd0,ddf)]

    #dry days
    dry=np.where(vldpp == 0,1.,0)
    dd0=(np.diff(dry)==1).nonzero()[0]+1; 
    ddf=(np.diff(dry)==-1).nonzero()[0]
    #note that ddf does not have the +1 because diff[i]=x[i+1]-x[i].
    if(dry[0]==1):dd0=np.insert(dd0,0,0)
    if(dry[-1]==1):ddf=np.insert(ddf,len(ddf),ntot-1)
    drylen=[(ii,jj,vldpp[int(ii):int(jj)+1]) for ii,jj in zip(dd0,ddf)]

    #Matching the characteristics to the periods
    for yr in range(nyrs):
        for ss in range(int(ns)):
            if(np.isfinite(ons[yr,ss,0])):
                #onset period
                di=ons[yr,ss,0]
                df=di+ons[yr,ss,1]                       
                
                #wet spell stats
                vld=[ii for ii in range(len(wetlen)) if ((wetlen[ii][0] >= di) & (wetlen[ii][1] <= df))] 

                if(len(vld) > 0):
                    #total precipitation during the onset period
                    totpp=sum([sum(wetlen[ii][2]) for ii in vld])
                    #note that totpp may differ from the total 
                    #precipitation calculated previously because it 
                    #does not include values below 1mm/day

                    #number of wet days
                    nrwet=sum([len(wetlen[ii][2]) for ii in vld])
                    #Precipitation rate: total precip during wet days 
                    #divided by duration of the onset period
                    onsStats[yr,ss,0]=totpp/ons[yr,ss,1]

                    #Precipitation intensity: total precipitation 
                    #during wet days divided by the number of wet days 
                    #during the onset period
                    onsStats[yr,ss,1]=totpp/float(nrwet)

                    #Fraction of wet days: number of wet days divided 
                    #by the duration of the onset period
                    onsStats[yr,ss,2]=100*(nrwet/ons[yr,ss,1])

                    #Mean duration of the wet spells = number of wet 
                    #days divided by the number of wet spells
                    onsStats[yr,ss,3]=nrwet/float(len(vld))
                    
                    #Longest wet spell in the onset period
                    onsStats[yr,ss,4]=max([len(wetlen[ii][2]) for ii in vld])

                else: onsStats[yr,ss,:5]=0.

                #Dry spell stats
                vld=[ii for ii in range(len(drylen)) if ((drylen[ii][0] >= di) & (drylen[ii][1] <= df))] 

                if(len(vld) > 0):
                    #number of dry days
                    nrdry=sum([len(drylen[ii][2]) for ii in vld])
                    
                    #Fraction of dry days: number of dry days divided 
                    #by the duration of the onset period
                    onsStats[yr,ss,5]=100*(nrdry/ons[yr,ss,1])
                    
                    #Mean duration of the dry spells = number of dry
                    #days divided by the number of dry spells
                    onsStats[yr,ss,6]=nrdry/float(len(vld))
                    #Longest dry spell in the onset period
                    onsStats[yr,ss,7]=max([len(drylen[ii][2]) for ii in vld])
                else: onsStats[yr,ss,5:]=0.
    
    return onsStats


############################
#Functions related to onset date analysis based on Bombardi et al 2019 (B19)
############################

def AnnualCycleB19(pp):

    #Annual cycle
    cycle=pp.groupby('time.dayofyear').mean(dim='time')
    
    yrfft=fft.rfft(cycle.values) 
    #FFT of the annual cycle. 
    yrfft=2*yrfft/cycle.dayofyear.size

    frcharm=cycle.dayofyear.size*(yrfft.real**2+yrfft.imag**2)/\
            (2*(cycle.dayofyear.size-1)*cycle.var(ddof=1).values)
    #variance explained by each harmonic

    #Checking if grid point has only one or two seasonal cycles: if the explained variance 
    #of the first harmonic is larger than the second or third harmonic, the grid point is 
    #analysed, otherwise it is masked out
    if ((frcharm[2] >= frcharm[1]) | (frcharm[3] >= frcharm[1])):
        w0=np.nan
    else:
        Harm=fft.irfft(yrfft[:2],cycle.dayofyear.size)*(cycle.dayofyear.size/2) #reconstructing the annual cycle using only the first harmonic
        w0=(np.argmin(Harm)+1)

    return cycle,w0


def SeasonOnsetDemiseB19(pp,dd0):
    
    ndays=int(365/2)  #serch minimum in the half year after (before if demise date) 
                      #the beginning of the water year (dd0)
    vld=pp.sel(time=pp.time.dt.dayofyear==dd0).time #selecting only days == dd0
    ons=[];dm=[]

    #Onset date
    for t0 in vld:
        if (pp.time[-1].item()-t0.item()).days > 15: #To avoid calculation with short 
                                                     #time series for the last year
    
            tt=ndays-1 if t0+dt(days=ndays-1) < pp.time[-1].item() \
                       else (pp.time[-1].item()-t0.item()).days
            #In the last year, use the number of days available

            if(tt > 40): #minimum number of days in the season
                ppcum=pp.sel(time=slice(t0,t0+dt(days=tt))).cumsum(dim='time')
                #Cumulative precipitation
                ppcum=butterworth(ppcum,cutoff=18) 
                #Smoothing the cumulative precipitation with a low-pass Butterworth filter 
                #with a cutoff of 15 days
                if(ppcum.argmin('time').item() < tt):
                    #Checking if the minimum occurs with 180 days from the beginning of the water year. if not, onset is considered as missing
                    ti=ppcum.isel(time=ppcum.argmin('time').item()).time
                    #Position of the minimum cumulative precipitation, i.e., the onset of the season
                else:
                    ti=xr.DataArray(np.nan)
                ons.append(ti)

    #Demise date
    ndays=int(2*365/3)  #serch minimum in the two thirds of the year before the start of the next water year - accounts for the fact that when serching backwards, it goes through the entire dry season before reaching the wet season
    for t0 in vld[::-1]: #Starting from the end of the year
        if (t0.item()-pp.time[0].item()).days > 15: #To avoid calculation with short time
                                                    #series in the first year (going backwards)
            tt=ndays-1 if t0-dt(days=ndays-1) > pp.time[0].item() \
                       else (t0.item()-pp.time[0].item()).days
            #In the first year, use the number of days available
            
            if(tt > 40): #minimum number of days in the season
                ppcum=pp.sel(time=slice(t0-dt(days=tt),t0))[::-1].cumsum(dim='time')
                #Cumulative precipitation
                ppcum=butterworth(ppcum,cutoff=18) 
                #Smoothing the cumulative precipitation with a low-pass Butterworth filter 
                #with a cutoff of 15 days
                
                if(ppcum.argmin('time').item() < tt):
                    ti=ppcum.isel(time=ppcum.argmin('time').item()).time
                else:
                    ti=xr.DataArray(np.nan)
                #Position of the minimum cumulative precipitation, i.e., the demise of the season

                dm.append(ti)
    dm=dm[::-1] #Reversing the list of demise dates to have the same order as the onset dates

    return ons,dm


########################
#Auxiliary functions
###########################

def DelCandidate(rmv,firstday,lenseason):
    '''
        Auxiliary function to remove the items listed in rmv from the 
        firstday and lenseason lists and returns updated lists 

        INPUT:
            rmv - list with the indices to be removed
            firstday - list with the first day of onset periods
            lenseason - list with the length of onset periods

        RETURNS:
            firstday - updated list with the first day of onset periods
            lenseason - updated list with the length of onset periods
    '''

    fd,ln=firstday.copy(),lenseason.copy()
    if(len(rmv) > 0):
        rmv.sort()
        for ii in rmv[::-1]:
            ln.pop(ii)
            fd.pop(ii)
    return fd,ln


def butterworth(indat,cutoff,ftype='lowpass'):
    '''
        lowpass is the default, but can be also 'highpass' , 'bandpass' or 'bandstop'
    '''
    from scipy.signal import butter, filtfilt, sosfilt
 
    if (isinstance(indat, np.ndarray)):
        padd1 = indat[:int(2*cutoff)][::-1]
        padd2 = indat[-int(2*cutoff):][::-1]
        paddannual = np.hstack((padd1,indat,padd2))
    elif (isinstance(indat, xr.DataArray)):
        padd1 = indat.isel(time=slice(None,int(2*cutoff))).values[::-1]
        padd2 = indat.isel(time=slice(indat.time.size-int(2*cutoff),None)).values[::-1]
        paddannual = np.hstack((padd1,indat.values,padd2))
    else:
        print('Input data type not recognized. Must be numpy array or xarray DataArray')
        return None

    outdat=np.zeros(paddannual.shape); outdat[:]=np.nan

    if(np.any(np.isfinite(paddannual))):
        b, a = butter(2, 1/float(cutoff),btype=ftype)
        y = filtfilt(b, a, paddannual[~np.isnan(paddannual)]) #applies the filter foreward and backward in time to avoid phase shift.
        outdat[~np.isnan(paddannual)]=y

    outdat=outdat[int(2*cutoff):-int(2*cutoff)]
    if isinstance(indat, np.ndarray):
        return outdat
    else:
        return xr.DataArray(outdat, coords=indat.coords, dims=indat.dims)



def CircularQuantile(indat,ndays,qq,dim=['years']):

    def circquantile(indat,ndays,qq):
        '''
        Here indat is a 1D np array
        '''

        indat=np.ravel(indat)  # Flatten the input array

        if np.any(np.isfinite(indat)):

            tmp=np.sort(indat)

            best_span=[np.ptp((tmp-tmp[jj])%ndays) for jj in range(len(tmp))]
            best_span=np.where(np.isnan(best_span),ndays+1,best_span) # to avoid NaN in the argmin
            pst=np.argmin(best_span)
            linear=(tmp-tmp[pst])%ndays
            
            qtllinear = np.nanquantile(linear,qq)
            qtl =(qtllinear+tmp[pst])%ndays

        else:
            qtl = np.full_like(qq, np.nan, dtype=float)

        return qtl
    
    outdat=xr.apply_ufunc(
        circquantile,
        indat,
        input_core_dims=[dim],
        output_core_dims=[["quantile"]],
        kwargs={"ndays": ndays, "qq": qq},
        vectorize=True,
        output_dtypes=[float],
    )
    outdat=outdat.assign_coords(quantile=qq)

    return outdat

def CircularMedian(indat,ndays,dim=['years']):

    out=CircularQuantile(indat,ndays,[0.5],dim=dim)

    out=out.sel(quantile=0.5)  #Selecting the median

    return out


def CircularMean(indat,ndays,dim=['years']):

    angle=(indat/ndays)*2*np.pi
    # Convert to Cartesian coordinates
    angx=np.cos(angle); angy=np.sin(angle)
    aveang=np.arctan2(angy.mean(dim=dim,skipna=True),\
                      angx.mean(dim=dim,skipna=True))
    aved=(aveang/(2*np.pi))*ndays
    aved=np.round(aved % ndays)
    aved=aved.where(aved > 0,1)  # Avoid zero value for first day of onset

    #Mask all missing
    aved=aved.where(indat.notnull().any(dim=dim),np.nan)

    return aved

def CircularVar(indat,ndays,dim=['years']):

    angle=(indat/ndays)*2*np.pi   
    # Convert to Cartesian coordinates
    angx=np.cos(angle); angy=np.sin(angle)
    # Compute mean x and y along the "time" axis
    mean_x = angx.mean(dim=dim,skipna=True)
    mean_y = angy.mean(dim=dim,skipna=True)    
    # Compute resultant vector length R
    R=np.sqrt(mean_x**2 + mean_y**2)
    
    circstd=np.sqrt(-2*np.log(R))*(ndays/(2*np.pi))  # Standard deviation in days
    circvar=circstd**2  # Squared to get variance    # Compute circular variance

    return circvar
 

def poly2mask(indat,poly):

    from matplotlib.path import Path
    
    poly=np.asarray(poly).T  #poly is the list of coordinates from previous function
    if ('lat' in indat.dims):
        msk=indat.where(((indat.lon >= poly[:,0].min()) & (indat.lon <= poly[:,0].max())) &\
                        ((indat.lat >= poly[:,1].min()) & (indat.lat <= poly[:,1].max())))
        tmp=msk.stack(z=('lon','lat'))
    else:
        msk=indat.where(((indat.longitude >= poly[:,0].min()) & (indat.longitude <= poly[:,0].max())) &\
                        ((indat.latitude >= poly[:,1].min()) & (indat.latitude <= poly[:,1].max())),drop=True)
        tmp=msk.stack(z=('longitude','latitude'))

    pts=tmp.where(tmp.notnull(),drop=True).coords['z'].values
    pts=np.asarray([[ii,jj] for ii,jj in pts]) #converting from tuple to array

    path=Path(poly)

    if(len(pts) > 0):
        pvld=path.contains_points(pts)
        rmv=pts[pvld == False]  #points outside the polygon
        if(len(rmv) > 0):
            for xx,yy in rmv:
                #Remove points outside the polygon
                if ('lat' in indat.dims):
                    msk.loc[{'lon':xx,'lat':yy}]=np.nan
                else:
                    msk.loc[{'longitude':xx,'latitude':yy}]=np.nan  

    return msk
