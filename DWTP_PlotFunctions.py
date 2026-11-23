import numpy as np
import sys
import xarray as xr
import DWTP_auxiliary as optb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mplticker
from matplotlib.colors import BoundaryNorm, ListedColormap
from cartopy import feature
import cftime
import seaborn as sns
import cmocean as cm
import colorcet as cc

def MapParam(varn,ndays,ptype='ave'):
    ''' 
        Returns parameters used in Map plot according to variable
    '''
    
    
    match ptype:
        case 'ave':
            match varn:
                case 'NrRainy':
                    shd='RdBu'; shdave=(); unit='';ctl=''
                case v if v in ['FirstDay','Onset','Demise']:
                    #shd=ListedColormap(sns.color_palette("hsv",90,desat=.6).as_hex())
                    shd=plt.get_cmap("cet_cyclic_mybm_20_100_c48")
                    shdave=(1,ndays+1,10,10)
                    unit=''; ctl='navy'
                case 'Length':
                    shd='cividis'; shdave=(20,140,10,20)
                    unit='(days)'; ctl='xkcd:orange'
                case v if v in ['Intensity','OnsIntensity']:
                    shd='BuPu'; shdave=(0,1.6,.2,.4)
                    unit='($x10^{-1} mm.day^{-2}$)'; ctl='navy'
                case v if v in ['TotPrecip','IntensityPP']:
                    shd='BuPu';unit='(mm)'; ctl='navy'
                    shdave=(0,700,50,100) if v == 'TotPrecip' else (0,3000,100,500)
                case v if v in ['PrecipRate','PrecipInt']:
                    shd='YlGnBu'; unit='($mm.day^{-1}$)'; ctl='navy'
                    shdave=(0,12,1,2) if v == 'PrecipRate' else (0,20,2,4)
                case 'FrcWetDays':
                    shd='BuPu'; shdave=(0,100,10,20)
                    unit='(%)'; ctl='xkcd:royal purple'
                case v if v in ['LengthWetSpells','MaxWetSpell']:
                    shd='BuPu'; unit='(days)'; ctl='xkcd:royal purple'
                    shdave=(0,6,1,1) if v == 'LengthWetSpells' else (0,20,2,4)
                case 'FrcDryDays':
                    shd='YlOrBr'; shdave=(0,100,10,20)
                    unit='(%)'; ctl='xkcd:dark brown'
                case v if v in ['LengthDrySpells','MaxDrySpell']:
                    shd='YlOrBr'; unit='(days)'; ctl='xkcd:dark brown'
                    shdave=(0,16,2,4) if v == 'LengthDrySpells' else (0,30,2.5,5)
                case 'Duration':
                    shd='BuPu';shdave=(0,210,10,30);unit='(days)'; ctl='xkcd:royal purple'
                case 'DurationDry': 
                    shd='YlOrBr';shdave=(0,210,10,30);unit='(days)'; ctl='xkcd:dark brown'
                case 'TotalDry':
                    shd='BuPu';shdave=(0,1200,100,300);unit='(mm)'; ctl='navy'
                case 'FrcAmpl':
                    shd='BuPu';shdave=(0,100,10,20);unit='';ctl='navy'
                case 't2mMean':
                    shd='cmo.thermal';shdave=(10,36,2,4);unit='($^oC$)';ctl='green'
                case 't2mRate':
                    shd='cmo.balance';shdave=(-3,3,0.5,1);unit='($x10^{-1}$'+'$^oC.day^{-1}$)';ctl='navy'
                case 't2mFrcXtr':
                    shd='YlOrBr';shdave=(0,20,2,4);unit='(%)';ctl='navy'
                case 'smMean':
                    shd='cmo.rain';shdave=(0,0.6,.05,.1);ctl='navy'
                    unit='($m^3.m^{-3}$)'
                case 'smRate':
                    shd='cmo.balance_r';shdave=(-4,4,0.5,1); ctl='navy'
                    unit='($x10^{-3}m^3.m^{-3}.day^{-1}$)'
                case 'smFrcXtr':
                    shd='YlOrBr';shdave=(0,20,2,4);unit='(%)';ctl='navy'
                case _: 
                    print ('Unkown variable')
                    shd=''; shdave=()
                    unit=''; ctl=''

        case 'diff':
            if varn in ['NrRainy','Demise','Length','Duration','Intensity','OnsIntensity',\
                        'TotPrecip','PrecipRate','PrecipInt','IntensityPP','TotalDry',\
                        'FrcWetDays','LengthWetSpells','MaxWetSpell','FrcAmpl','smMean','smRate']:
                #red to blue
                shd=sns.diverging_palette(5,240,s=80,l=45,sep=4,center='light',as_cmap=True) 
            else: 
                #blue to red
                shd=sns.diverging_palette(240,5,s=80,l=45,sep=4,center='light',as_cmap=True)

            match varn:
                case 'NrRainy':
                    shdave=(-1,2,1,1); lvl=(); unit=''; ctl='navy'
                case v if v in ['FirstDay','Onset','Demise','Duration']: #Blue to red
                    shdave=(-60,60,10,20);unit='(days)';ctl='navy';lvl=(-60,60,10)
                case 'Length':#Red to blue
                    shdave=(-60,60,10,20); lvl=(20,120,20); unit='(days)'; ctl='navy'
                case v if v in ['Intensity','OnsIntensity']:#Red to blue
                    shdave=(-1.5,1.5,.25,.5); lvl=(0,3,.5);unit='($x10^{-1} mm.day^{-2}$)';ctl='navy'
                case 'TotPrecip':#Red to blue
                    shdave=(-300,300,50,100); lvl=(0,500,100);unit='(mm)'; ctl='navy'
                case 'PrecipRate':#Red to blue
                    shdave=(-9,9,1,3); lvl=(0,12,2); unit='($mm.day^{-1}$)';ctl='navy'
                case 'PrecipInt':#Red to blue
                    shdave=(-15,15,2.5,5); lvl=(0,20,4);unit='($mm.day^{-1}$)'; ctl='navy'
                case 'FrcWetDays':#Red to blue
                    shdave=(-45,45,3,9); lvl=(0,100,20);unit='(%)'; ctl='navy'
                case 'LengthWetSpells':#Red to blue
                    shdave=(-6,6,1,2); lvl=(0,6,1);unit='(days)'; ctl='navy'
                case 'MaxWetSpell':#Red to blue
                    shdave=(-15,15,2.5,5); lvl=(0,20,4);unit='(days)'; ctl='navy'
                case 'FrcDryDays':#Blue to red
                    shdave=(-45,45,3,9); lvl=(0,100,20);unit='(%)'; ctl='navy'
                case 'LengthDrySpells':#Blue to red
                    shdave=(-12,12,2,4); lvl=(0,16,4);unit='(days)'; ctl='navy'
                case 'MaxDrySpell':#Blue to red
                    shdave=(-20,20,2.5,5); lvl=(0,30,5);unit='(days)'; ctl='navy'
                case 'IntensityPP': #Red to blue 
                    shdave=(-100,100,10,20);unit='(%)'; ctl='navy';lvl=(-100,100,10)
                case 'DurationDry': #Blue to red
                    shdave=(-60,60,10,20);unit='(days)'; ctl='navy';lvl=(-60,60,10)
                case 'TotalDry': #Red to blue
                    shdave=(0,1200,100,300);unit='(%)'; ctl='navy';lvl=(0,1200,100)
                case 'FrcAmpl': #red to blue
                    shdave=(-100,100,10,20);unit='(%)';ctl='navy';lvl=(-100,100,10)
                case 't2mMean': #Blue to Red
                    shdave=(-10,10,1,2); lvl=(0,40,4);unit='($^oC$)'; ctl='navy'
                case 't2mRate': #Blue to Red
                    shdave=(-3,3,0.5,1);unit='($x10^{-1}$'+'$^oC.day^{-1}$)';ctl='navy'
                    lvl=(-3,3,0.5)
                case 't2mFrcXtr': #Blue to Red
                    shdave=(-40,40,5,10);unit='(%)';ctl='navy'; lvl=(0,20,4,4)
                case 'smMean': #red to blue
                    shdave=(-1,1,0.1,0.2); lvl=(0,0.6,.2,.2);ctl='navy'
                    unit='($m^3.m^{-3}$)'
                case 'smRate': #red to blue
                    shdave=(-4,4,0.5,1); lvl=(-4,4,1,1); ctl='navy'
                    unit='($x10^{-3}m^3.m^{-3}.day^{-1}$)'
                case 'smFrcXtr': #Blue to Red
                    shdave=(-10,10,1,2);lvl=(0,40,10,10);unit='(%)';ctl='navy'
                case _: 
                    print ('Unkown variable')
                    shd=''; shdave=()
                    unit=''; ctl=''; lvl=()       
        
        case 'absdiff':
            match varn:
                case 'BestInit':
                    shd=ListedColormap(sns.color_palette("hsv",90,desat=.6).as_hex()) #Enough to have one colour each 10 days
                    shdave=(1,13,1,1); lvl=(); unit='(month)'; ctl='navy'
                case 'LeadTime':
                    shd='Set2';shdave=(0,7,1,1); lvl=(); unit='(months)'; ctl='navy'
                case 'Overlap':
                    shd='YlGnBu';shdave=(0,100,5,10); lvl=(); unit='(%)'; ctl='navy'
                case 'FirstDay': #Blue to red
                    shd='YlOrRd';shdave=(0,60,5,20); unit=''; ctl='navy'; lvl=()
                case 'Length':#Red to blue
                    shd=sns.diverging_palette(5,240,s=80,l=45,sep=4,center='light',as_cmap=True)
                    shdave=(-30,30,5,10); lvl=(); unit='(days)'; ctl='navy'
                case v if v in ['Intensity','OnsIntensity']:#Red to blue
                    shd=sns.diverging_palette(5,240,s=80,l=45,sep=4,center='light',as_cmap=True)
                    shdave=(-20,20,5,10); lvl=();unit='($x10^{-2} mm.day^{-2}$)';ctl='navy'
                case 'TotPrecip':#Red to blue
                    shd=sns.diverging_palette(5,240,s=80,l=45,sep=4,center='light',as_cmap=True)
                    shdave=(-300,300,50,100); lvl=();unit='(mm)'; ctl='navy'
                case v if v in ['PrecipRate','PrecipInt']:#Red to blue
                    shd=sns.diverging_palette(5,240,s=80,l=45,sep=4,center='light',as_cmap=True)
                    lvl=(); unit='($mm.day^{-1}$)';ctl='navy'
                    shdave=(-12,12,2,4) if v == 'PrecipRate' else (-15,15,2.5,5)
                case 'FrcWetDays':#Red to blue
                    shd=sns.diverging_palette(5,240,s=80,l=45,sep=4,center='light',as_cmap=True)
                    shdave=(-60,60,10,20); lvl=();unit='(%)'; ctl='navy'
                case v if v in ['LengthWetSpells','MaxWetSpell']:#Red to blue
                    shd=sns.diverging_palette(5,240,s=80,l=45,sep=4,center='light',as_cmap=True)
                    lvl=();unit='(days)'; ctl='navy'
                    shdave=(-12,12,2,4) if v == 'LengthWetSpells' else (-18,18,3,4)
                case 'FrcDryDays':#Blue to red
                    shd=sns.diverging_palette(240,5,s=80,l=45,sep=4,center='light',as_cmap=True)
                    shdave=(-60,60,10,20); lvl=();unit='(%)'; ctl='navy'
                case v if v in ['LengthDrySpells','MaxDrySpell']:#Blue to red
                    shd=sns.diverging_palette(240,5,s=80,l=45,sep=4,center='light',as_cmap=True)
                    lvl=();unit='(days)'; ctl='navy'
                    shdave=(-12,12,2,4) if v == 'LengthDrySpells' else (-18,18,3,4)
                case _: 
                    print ('Unkown variable')
                    shd=''; shdave=()
                    unit=''; ctl=''; lvl=()       
        case _:
            print ('Unkown plot type')
            shd=''; shdave=(); unit=''; ctl=''

    if(ptype == 'diff'):  return shd, shdave, unit,ctl,lvl
    else: return shd, shdave, unit,ctl


def Map(sub,indat,mcol,colbar,shd,label,outf,sglvl=0.05,whiteout=False,cbar=True):

    domains={'SAfr': ((-35.0,0.0),(10.0,55.0)),
             'Afr': ((-35.0,15.0),(-18,52.0)),
             'EAfr': ((-10.0,15.0),(25.0,55.0)),
             'SAm':((-60.0,15.0),(275.,330.)),
             'SAAG':((-58.,16.0),(261.,343.)),
             'SAAG-WRF':((-58.,16.0),(261.,343.)),
             'UKMO':((-41.,16.),(274.,331.)),
             }

    cmap=plt.get_cmap(colbar)
    levels=np.arange(shd[0],shd[1]+0.1*shd[2],shd[2]) if (isinstance(shd,tuple)) else shd
    
    norm=BoundaryNorm(levels,ncolors=cmap.N,clip=True)

    plt.ioff()

    minlat,maxlat,minlon,maxlon=domains[sub][0][:]+domains[sub][1][:]
    match sub:
        case 'SAfr': fig = plt.figure(figsize=(3.,2.8)); pst=[0.07,0.17,.9,.82]
        case 'Afr': fig = plt.figure(figsize=(3.,2.8)); pst=[0.07,0.19,.9,.8]
        case 'EAfr': fig = plt.figure(figsize=(5.6,5.1)); pst=[0.03,0.12,.98,.87]
        case 'SAm': fig=plt.figure(figsize=(2.4,2.8)); pst=[0.06,0.06,.8,.92]
        case 'UKMO': fig=plt.figure(figsize=(4,3.8)); pst=[0.03,0.12,.98,.87]
        case 'SAAG' | 'SAAG-WRF': plt.figure(figsize=(4,4.2)); pst=[0.03,0.12,.98,.87]
        case _: print ('Domain not found')

    #basemap
    # PlateCarree (equidistant cylindrical projection)
    ax1 = plt.axes(projection=ccrs.PlateCarree())
    ax1.set_extent([minlon,maxlon,minlat,maxlat], crs=None)
    ax1.set_position(pst)
    ax1.add_feature(feature.BORDERS, linestyle='-', linewidth=0.25,edgecolor=mcol)
    ax1.coastlines(resolution='110m', linewidth=0.5,color=mcol)
    gr=ax1.gridlines(draw_labels=True,linewidth=0.5, color='black', alpha=0.5, linestyle=':')
    gr.top_labels = False
    gr.right_labels = False
    gr.ylabel_style = {'rotation':89,'size':8}
    gr.xlabel_style = {'size':8}
    gr.ylocator = mplticker.FixedLocator(np.arange(-60, maxlat+1, 10))
    gr.yformatter = LATITUDE_FORMATTER
    
    plt.box(on=None)

    #shades
    if('shd' in indat):
        lon,lat=indat['shd']['lon'],indat['shd']['lat']
        ax1.pcolormesh(lon,lat,indat['shd']['ave'],transform=ccrs.PlateCarree(),shading='auto',cmap=cmap,norm=norm)

    #Signifincance 
    if('sct' in indat):
        kk=np.asarray(['sct' in list(indat.keys())[i] for i in range(len(indat))]).nonzero()[0]
        for i in kk:
            dat=indat[list(indat.keys())[i]]
            lon,lat=dat['lon'],dat['lat']

            if(whiteout): #white shade over not significant areas
                if('ccol' in dat): ccol=dat['ccol']
                else: ccol='binary'
                tmp=np.where(dat['ave'] > sglvl,1,np.nan)
                ax1.pcolormesh(lon,lat,tmp,cmap=ccol,transform=ccrs.PlateCarree(),alpha=.75)
            else:   #Black shade in significant areas
                if('ccol' in dat): ccol=dat['ccol']
                else: ccol='binary_r'
                tmp=np.where(dat['ave'] <= sglvl,1,np.nan)
                ax1.pcolormesh(lon,lat,tmp,cmap=ccol,transform=ccrs.PlateCarree(),alpha=.35)

    #Contours
    if('ctl' in indat):
        kk=np.asarray(['ctl' in list(indat.keys())[i] for i in range(len(indat))]).nonzero()[0]
        for i in kk:
            dat=indat[list(indat.keys())[i]]
            lon,lat=dat['lon'],dat['lat']

            if(isinstance(dat['lvl'],tuple)):
                lvls=np.arange(dat['lvl'][0],dat['lvl'][1],dat['lvl'][2])
            else:lvls=dat['lvl']

            if(not np.all(np.isnan(dat['ave']))):
                cto=ax1.contour(lon,lat,dat['ave'],colors=dat['col'],levels=lvls,linewidths=.5,transform=ccrs.PlateCarree())
                if (len(cto.levels) > 1): fm='%i' if ((cto.levels[1]-cto.levels[0]) % 1 == 0) else '%1.2f'
                else: fm='%i' if (cto.levels[0] % 1 == 0) else '%1.2f'
                fts=6
                plt.clabel(cto,fontsize=fts,fmt=fm,inline=1)
                if ('extra' in dat): 
                    ax1.contour(lon,lat,dat['ave'],colors=dat['col'],levels=dat['extra'],linewidths=1,transform=ccrs.PlateCarree())
        
    #Polygons
    if('poly' in indat):
        kk=np.asarray(['poly' in list(indat.keys())[i] for i in range(len(indat))]).nonzero()[0]
        for i in kk:
            dat=indat[list(indat.keys())[i]]
            x,y=dat['poly']
            
            plt.plot(x,y,':',color=dat['ccol'],linewidth=1.5,transform=ccrs.PlateCarree())

    #Extra text
    if('txt' in indat):
        kk=np.asarray(['txt' in list(indat.keys())[i] for i in range(len(indat))]).nonzero()[0]
        for i in kk:
            dat=indat[list(indat.keys())[i]]
            tx,ty=dat['lon'],dat['lat']
            tx1=plt.text(tx,ty,dat['txt'],va='bottom',ha='left',fontsize=6,transform=ccrs.PlateCarree())

    #Extra points
    if('pts' in indat):
        tx,ty=indat['pts']['coords'][:,0],indat['pts']['coords'][:,1]
        plt.plot(tx,ty,' o',color=indat['pts']['pcol'],markersize=1,transform=ccrs.PlateCarree())


    #colorbar
    if(cbar):
        if(sub in ['SAm']): #vertical bar
            cbaxes = fig.add_axes([0.84, 0.06, 0.03, 0.92])
        else: #horizontal bar
            cbaxes = fig.add_axes([0.07, 0.09, 0.9, 0.03])

        lvticks=np.arange(shd[0],shd[1]+shd[2],shd[2]) if(isinstance(shd,tuple)) else levels
        psticks=np.linspace(0,1,len(lvticks))
        dv=shd[2] if isinstance(shd,tuple) else 1./(len(shd)-1)
        if(dv < 1): dv=np.diff(psticks)[0]

        if('cbar' in indat):
            ptxticks=indat['cbar']['pst']
            txticks=indat['cbar']['txt']
        else:
            txticks=np.round(np.arange(shd[0],shd[1]+shd[2],shd[3]),2) if(isinstance(shd,tuple)) else levels
            if((txticks[1]-txticks[0])  % 1 ==0): txticks=txticks.astype(int)
            ptxticks=np.linspace(0,1,len(txticks))

        fts=7
        if(sub in ['SAm']): #vertical bar
            plt.ylim(0,1)
            for ii,lv in enumerate(psticks): 
                plt.fill_between([0.5,1],[lv,lv],[lv+dv,lv+dv],\
                                 color=cmap(norm(lvticks[ii])))
                
            plt.xticks([])
            cbaxes.yaxis.tick_right()
            plt.yticks(ptxticks,txticks,fontsize=fts,rotation=90,va='center')

            cbaxes.tick_params(axis='y',which='both',length=1,pad=2.4)  
            cbaxes.yaxis.set_label_position('right')
            cbaxes.set_ylabel(label,fontsize=fts,labelpad=1.8)


        else:
            plt.xlim(0,1)
            for ii,lv in enumerate(psticks): 
                plt.fill_betweenx([0.5,1],[lv,lv],[lv+dv,lv+dv],color=cmap(norm(lvticks[ii])))
                
            plt.yticks([])
            cbaxes.yaxis.tick_right()
            #Julian day of the first of each month
            plt.xticks(ptxticks,txticks,fontsize=fts,ha='center')

            cbaxes.tick_params(axis='x',which='both',length=1,pad=1.5)  
            cbaxes.set_xlabel(label,fontsize=fts,labelpad=0.5)


    outfile=outf+'.png'
    print(outfile)
    plt.savefig(outfile,dpi=150,format="png")
    plt.close()


def GraphOnsetPeriods(msk,nsea,yr0,ons,ave,lvls,mycol,title,clbl,outf,\
                        ylab='PP ($mm.day^{-1}$)',b19=False,extra=False):
    
    nyrs,ndays=msk.shape

    if nyrs > 20: dyr=10
    elif nyrs > 10: dyr=5
    else:dyr=2

    dd=[0,31,59,90,120,151,181,212,243,273,304,334] if ndays > 360 else list(range(0,ndays,30) )
    dtxt=['J','F','M','A','M','J','J','A','S','O','N','D'] 

    fig=plt.figure(figsize=[3.6,2.5])
    ax1=fig.add_axes([.11,.24,.77,.68])

    levels=list(np.arange(lvls[0],lvls[1]+lvls[2],lvls[2]))
    cmap=plt.get_cmap(mycol)
    norm=BoundaryNorm(levels,ncolors=cmap.N,clip=True)

    plt.pcolormesh(range(ndays+1),range(nyrs+1),msk,cmap=cmap,\
                   norm=norm,shading='auto')

    if(ons.shape[0] == nsea):
        for ss in range(nsea):
            plt.plot([ons[ss,0],ons[ss,0]-1],[0,nyrs],':k') #ons in jday (01/01 == 1)
            ddf=ons[ss,0]+ons[ss,1]
            if(ddf >= ndays): ddf=ddf-ndays
            plt.plot([ddf,ddf],[0,nyrs],':k')
    if(ons.shape[0] >= nyrs):
        plt.contour(range(ndays),range(nyrs+1),ons,levels=np.arange(0,0.55,.1),\
                    cmap='viridis_r',linewidths=.75)
        plt.clim(0,0.5)

    plt.ylim(0,nyrs)
    yrs=np.arange(yr0,yr0+nyrs)
    vyr=(yrs % dyr == 0).nonzero()[0]
    plt.yticks(vyr+.5,yrs[vyr],fontsize=8,rotation=90,va='center')
    plt.ylabel('YEARS',fontsize=8)
    plt.xlim(0,ndays)
    plt.xticks(dd,dtxt,fontsize=8)
    #plt.xlabel('DOY',fontsize=7)
    plt.title(title,fontsize=8)
    
    #If b19==True - plto onset day accoding to Bombardi 2019
    if(not isinstance(b19,bool)):
        for yr in range(nyrs):
            yi=yr; xi=b19[yr]
            if(xi >= ndays): xi=xi-ndays; yi=yr+1
            plt.plot(np.repeat(xi,2),[yi,yi+1],color='orange',alpha=.5,linewidth=1.5)
        plt.plot([ons[0,2],ons[0,2]],[-dyr/10,-dyr/10],'^',color='orange',markersize=5,clip_on=False)


    if isinstance(extra,np.ndarray):
        for yr in range(nyrs):
            for ss in range(nsea):
                lty=':' if ss == 0 else '--'
                if(np.isfinite(extra[yr,ss,0])):
                    yi=yr; xi=extra[yr,ss,0]
                    if(xi > ndays): xi=xi-ndays; yi=yr+1

                    plt.plot(np.repeat(xi,2),[yi,yi+1],lty,color='k',alpha=.5,linewidth=2)

                    xf=xi+extra[yr,ss,1]
                    if(xf > ndays): xf=xf-ndays; yi=yr+1

                    plt.plot(np.repeat(xf,2),[yi,yi+1],lty,color='k',alpha=.5,linewidth=2)

    #Plot climatological cycle as a right y-axis, adjusting to fit in nyrs (max value)
    ax2=ax1.twinx()
    plt.plot(ave,'navy')
    mx=np.ceil(np.nanmax(ave)/10)*10.
    mn=np.floor(np.nanmin(ave)/5)*5.
    if(mn < 0): mn=0
    dx=2 if (mx-mn) < 20 else 5
    mn=0; mx=12; dx=2
    ax2.set_ylim(mn,mx)   
    ax2.set_yticks(list(range(0,int(mn),-dx))[::-1][:-1]+list(range(0,int(mx)+1,dx)))
    ax2.tick_params(axis='y',labelcolor='navy',labelsize=8)
    ax2.set_ylabel(ylab,color='navy',fontsize=8)
    ax2.spines['right'].set_color('navy')

    #Colorbar
    cbaxes = fig.add_axes([0.11, 0.12, 0.77, 0.03])
    plt.xlim(0,1)
    psticks=np.linspace(0,1,len(levels))
    dv=1
    for ii,lv in enumerate(psticks): 
        plt.fill_betweenx([0.5,1],[lv,lv],[lv+dv,lv+dv],\
                          color=cmap(norm(levels[ii])))
        
    plt.yticks([])
    cbaxes.yaxis.tick_right()
    xtcks=list(np.arange(lvls[0],lvls[1]+lvls[2],lvls[3]))
    plt.xticks(np.linspace(0,1,len(xtcks)),np.round(xtcks,1),fontsize=8,ha='center')

    cbaxes.tick_params(axis='x',which='both',length=1,pad=2.4)  
    cbaxes.set_xlabel(clbl,fontsize=8,labelpad=1.8)

    outfile=outf+'.png'
    print(outfile)
    plt.savefig(outfile,format='png',dpi=300)
    plt.close()

def GraphHovmoller(ons,ave,mycol,rng,title,dims,outf,extra=False):

    xvl=ons[dims].values
    yvl=ons.years.values
    ons=ons.values
    
    if yvl.size > 20: dyr=10
    elif yvl.size > 10: dyr=5
    else:dyr=2

    mn,mx=rng

    fig=plt.figure(figsize=[3.6,2.5])
    ax1=fig.add_axes([.11,.23,.79,.7])

    levels=np.arange(mn,mx+1,5)
    cmap=plt.get_cmap(mycol)
    norm=BoundaryNorm(levels,ncolors=cmap.N,clip=True)

    plt.pcolormesh(xvl,yvl,ons,cmap=mycol,norm=norm,shading='auto')
    if(dims == 'lat'): ax1.invert_xaxis()  #latitudes are negative in the south

    plt.ylim(yvl[0]-.5,yvl[-1]+.5)
    pvld=(yvl%dyr == 0).nonzero()[0]
    
    plt.yticks(yvl[pvld],yvl[pvld],fontsize=8,rotation=90,va='center')
    plt.ylabel('YEAR',fontsize=8)

    plt.title(title,fontsize=8)

    if((xvl > 180).any()):  #Only when using lon as main transect
        xlb=np.asarray([(((yy+180)%360)-180) for yy in xvl])  #-180 to 180 dg
    else: xlb=xvl

    if np.ptp(abs(xlb)) > 60: dx=10
    elif np.ptp(abs(xlb)) > 20: dx=5
    elif np.ptp(abs(xlb)) > 5: dx=3
    else: dx=1
    
    xi=int(np.ceil(np.min(xlb)/dx)*dx)
    xf=int(np.floor(np.max(xlb)/dx)*dx)
    xlb=np.arange(xi,xf+dx,dx)
    if(dims == 'lon'):
        lbl=[f"{abs(ii)}{'W' if ii < 0 else 'E'}" for ii in xlb]
    else:
        lbl=[f"{abs(ii)}{'S' if ii < 0 else 'N'}" for ii in xlb]

    xps=xlb+360 if (np.any(xlb < 0) & (dims == 'lon')) else xlb
    plt.xticks(xps,lbl,fontsize=8)

    #line plot as a second y-axis
    ax2=ax1.twinx()
    plt.plot(xvl,ave,color='navy')
    ax2.set_ylim(mn,mx)
    ax2.tick_params(axis='y',labelsize=8,color='navy',labelcolor='navy')
    ax2.set_ylabel('DOY',color='navy',fontsize=8)
    ax2.spines['right'].set_color('navy')

    if(isinstance(extra,bool)):
        pass
    else: #Plotting a second line representing the b19 onset day
        plt.plot(xvl,extra,color='navy',linestyle='--')

    dd0=np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365,396])
    mi=np.argmin(abs(mn-dd0))
    mf=np.argmin(abs(mx-dd0))

    dd=np.asarray([cftime.num2date(ii,calendar='noleap',units='days since 1970-01-01').day\
                    for ii in list(range(dd0[mi],dd0[mf]+1))])
    mm=np.asarray([cftime.num2date(ii,calendar='noleap',units='days since 1970-01-01').month \
                    for ii in list(range(dd0[mi],dd0[mf]+1))])
    doy=np.asarray([cftime.num2date(ii,calendar='noleap',units='days since 1970-01-01').dayofyr\
                    for ii in list(range(dd0[mi],dd0[mf]+1))])
        
    pst=(dd == 1).nonzero()[0]

    dtxt=['%01d/%01d' % (dd[ii],mm[ii]) for ii in pst]
    ax2.set_yticks(doy[pst],dtxt,rotation=90,fontsize=8, va='center')

    #Colorbar
    cbaxes = fig.add_axes([0.11, 0.11, 0.79, 0.03])
    plt.xlim(0,1)
    psticks=np.linspace(0,1,len(levels))
    dv=1
    for ii,lv in enumerate(psticks): 
        plt.fill_betweenx([0.5,1],[lv,lv],[lv+dv,lv+dv],\
                          color=cmap(norm(levels[ii])))
        
    plt.yticks([])
    cbaxes.yaxis.tick_right()
    pst=((dd == 1) | (dd==15)).nonzero()[0]
    xtcks=(doy[pst]-mn)/(mx-mn)
    dtxt=['%01d/%01d' % (dd[ii],mm[ii]) for ii in pst]
    plt.xticks(xtcks,dtxt,fontsize=7,ha='center')


    cbaxes.tick_params(axis='x',which='both',length=1,pad=2.4)  
    cbaxes.set_xlabel('DOY',fontsize=7,labelpad=1.8)

    outfile=outf+'.png'
    print(outfile)
    plt.savefig(outfile,format='png',dpi=300)
    plt.close()


def LineBox(indat,dset,lbl,ylab,nvar,mcol,outf):
    '''
        First array is plotted as a line plot, with interquartile range as envelope and
        minimum and maximum values as dotted lines
        Second array is plotted as a boxplot.
    '''
    obs=dset[0][0]; mdl=dset[1]
    seq=indat[lbl].values

    dx=.9/len(mdl)
    #width of the box

    fig=plt.figure(figsize=[2.5,2])
    if (len(mdl)+1 <= 3): ax1=fig.add_axes([.12,.2,.86,.76])   #caption: 1 row, 3 columns
    elif(len(mdl)+1 <= 6): ax1=fig.add_axes([.12,.25,.86,.71]) #caption: 2 row, 3 columns
    else: ax1=fig.add_axes([.12,.3,.86,.65]) #caption: 3 rows, 3 columns

    #lineplot + envelope only for obs dataset
    seqx=[-.5]+list(range(len(seq)))+[len(seq)-.5]
    dat=indat.sel(mdl=obs).values
    datx=np.vstack((dat[:1,:],dat,dat[-1:,:]))
    #median
    plt.plot(np.asarray(seqx),datx[:,2],linewidth=1,color="#899eaa",ms=4,\
                zorder=10,drawstyle='steps-mid')
    #box
    plt.fill_between(np.asarray(seqx),datx[:,1],datx[:,3],edgecolor='white',linewidth=0,\
                        facecolor='#899eaa',step='mid',alpha=.2)
    #Whiskies
    plt.plot(np.asarray(seqx),datx[:,0],':',linewidth=.5,color='#899eaa',\
                drawstyle='steps-mid')
    plt.plot(np.asarray(seqx),datx[:,-1],':',linewidth=.5,color='#899eaa',\
                drawstyle='steps-mid')

    #models for each area
    mn=indat.min().values; mx=indat.max().values
    if mn > 10: 
        mn=np.floor(mn/10)*10
        mx=np.ceil(mx/10)*10
        
    for im,m in enumerate(seq): 
        plt.plot([im+.5,im+.5],[mn,mx],'--',color='black',linewidth=.1)
        for ib,b in enumerate(mdl):
            xvert=[-.45+im+(ib*dx)+(.1*dx),-.45+im+((ib+1)*dx-0.1*dx)]
            if('FUT' in b):  #Fut scen
                Box(indat.sel(mdl=b,area=m).values,xvert,mcol[ib],contourOnly=True)
            else:
                Box(indat.sel(mdl=b,area=m).values,xvert,mcol[ib],contourOnly=False)

    plt.xlim(-0.5,len(seq)-.5)
    plt.xticks(range(len(seq)),seq,fontsize=6)
    ax1.tick_params(axis='x',which='major',length=1,pad=2.4)  
    plt.xlabel('Area',fontsize=6,labelpad=2)

    plt.ylim(mn,mx); 
    if(nvar in ['FirstDay','Onset']):
        dd0=np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365,396])
        # mi=np.argmin(abs(mn-dd0))
        # mf=np.argmin(abs(mx-dd0))

        dd=np.asarray([cftime.num2date(ii,calendar='noleap',units='days since 1970-01-01').day\
                       for ii in list(range(int(mn),int(mx)+1))])
        mm=np.asarray([cftime.num2date(ii,calendar='noleap',units='days since 1970-01-01').month \
                        for ii in list(range(int(mn),int(mx)+1))])
        
        if(len(dd) < 62): 
            pst=((dd == 1) | (dd == 10) | (dd == 20)).nonzero()[0]
        elif(len(dd) < 90): #labels as days 01 and 15 of each month
            pst=((dd == 1) | (dd == 15)).nonzero()[0]
        else: pst=(dd == 1).nonzero()[0]

        dtxt=['%02d/%02d' % (dd[ii],mm[ii]) for ii in pst]
        plt.yticks(mn+pst,dtxt)

    plt.yticks(fontsize=6,rotation=90,va='center')
    plt.ylabel(ylab,fontsize=6)
    ax1.tick_params(axis='y',which='major',length=2,pad=2.4)  

    #Ploting legend
    if (len(mdl)+1 <= 3): ax=fig.add_axes([0.01,.02,.98,0.05]); nrow=1
    elif(len(mdl)+1 <= 6): ax=fig.add_axes([0.01,.02,.98,0.1]); nrow=2
    else: ax=fig.add_axes([0.01,.02,.98,0.15]); nrow=3

    plt.xlim(0,1); plt.ylim(0,1)  
   
    ncol=int(np.ceil((len(mdl)+1)/float(nrow)))
    dy=1./nrow

    if(ncol <= 2): px=[0.01,.45]
    elif(ncol == 3): px=[0.01,.22,.61] 
    else: px=[0.01,.25,.52,.77]

    ri=0; ci=0
    #Obs
    plt.plot([px[ci],.045+px[ci]],[1-((ri+.5)*dy),1-((ri+.5)*dy)],linewidth=2,\
            color='#899eaa')
    plt.text(0.055+px[ci],1-((ri+.5)*dy),obs,va='center',fontsize=6)
    
    ri+=1
    if(ri == nrow): ri=0; ci+=1
    
    for ib,b in enumerate(mdl):
        if('FUT' in b):  #Fut scen
            fc='white'; ec=mcol[ib]; lw=1
        else:
            fc=mcol[ib]; ec='None'; lw=0

        plt.fill_between([px[ci],.045+px[ci]],[1.-dy*(ri+1)+.02,1.-dy*(ri+1)+.02],\
                         [1.-(dy*ri)-.02,1.-(dy*ri)-.02],facecolor=fc,edgecolor=ec,\
                         linewidth=lw,alpha=.75,zorder=10)
        plt.plot([px[ci],.045+px[ci]],[1-((ri+.5)*dy),1-((ri+.5)*dy)],linewidth=1,\
                    color=mcol[ib],zorder=2)

        plt.text(0.055+px[ci],1-((ri+.5)*dy),b,va='center',fontsize=6)       
        ri+=1
        if(ri == nrow): ri=0; ci+=1
    
    plt.axis('off')

    outfile=outf+'.png'
    print (outfile)
    plt.savefig(outfile,dpi=300,format="png")
    
    plt.close()

def Box(indat,ps,bcol,contourOnly=False):
    '''
        Auxiliary function to plot each box in boxplots
    '''

    mn,lw,md,up,mx=indat

    lst='-'; fc=bcol; ec='None'; lwd=0
    if(contourOnly): fc='white'; ec=bcol;lwd=1

    plt.fill_between(ps,lw,up,edgecolor=ec,facecolor=fc,alpha=.75,zorder=5,linewidth=lwd,linestyle=lst)
        
    plt.plot(ps,np.repeat(md,2),color=bcol,linewidth=1.5,zorder=8,linestyle=lst)
    #Adding 5th and 95 percentiles
    plt.plot(np.repeat(np.nanmean(ps),5),[mn,lw,np.nan,up,mx],color=bcol,linestyle=lst,\
             linewidth=.3)

def BoxH(indat,ps,bcol,contourOnly=False):
    '''
        Auxiliary function to plot each box in boxplots (horizontal version)
    '''

    mn,lw,md,up,mx=indat

    lst='-'; fc=bcol; ec='None'; lwd=0
    if(contourOnly): fc='white'; ec=bcol;lwd=.75

    plt.fill_betweenx(ps,lw,up,edgecolor=ec,facecolor=fc,alpha=.75,zorder=5,linewidth=lwd,linestyle=lst)
        
    plt.plot(np.repeat(md,2),ps,color=bcol,linewidth=1.5,zorder=8,linestyle=lst)
    #Adding 5th and 95 percentiles
    plt.plot([mn,lw,np.nan,up,mx],np.repeat(np.nanmean(ps),5),'|',color=bcol,linestyle=lst,\
             linewidth=.3)

def FDLengthOns(indat,dset,areacol,mdlcol,outf):

    obs=dset[0]; mdl=dset[1]

    fig=plt.figure(figsize=[2,3])
    if (len(mdl) <= 4): ax1=fig.add_axes([.08,.16,.9,.83])   #caption: 1 row, 3 columns
    elif(len(mdl) <= 6): ax1=fig.add_axes([.08,.2,.9,.78]) #caption: 2 row, 3 columns
    else: ax1=fig.add_axes([.08,.24,.9,.74]) #caption: 3 rows, 3 columns
   
    mdn=indat.sel(quantile=0.5)

    mn=np.floor((mdn.sel(nstats='FirstDay').min().item()-10)/10)*10
    mx=np.ceil((mdn.sel(nstats='FirstDay').max().item()+\
                mdn.sel(nstats='Length').max().item()-20)/10)*10

    for ia,aa in enumerate(indat.area.values):
        y0=ia*len(mdl); y1=y0+len(mdl)
        
        #Observations as a background square
        di=mdn.sel(mdl=obs,area=aa,nstats='FirstDay').item()
        df=di+mdn.sel(mdl=obs,area=aa,nstats='Length').item()
        plt.fill_between([di,df],[y0+.1,y0+.1],[y1-.1,y1-.1],edgecolor='k',\
                        facecolor=areacol[ia],alpha=.5,linewidth=.5)
        #Adding interquartile range as dotted lines 
        di=indat.sel(quantile=0.25,mdl=obs,area=aa,nstats='FirstDay').item()
        plt.plot([di,di],[y0+.1,y1-.1],':',color=areacol[ia],linewidth=.5)
        df=di+indat.sel(quantile=0.25,mdl=obs,area=aa,nstats='Length').item()
        plt.plot([df,df],[y0+.1,y1-.1],':',color=areacol[ia],linewidth=.5)
        
        di=indat.sel(quantile=0.75,mdl=obs,area=aa,nstats='FirstDay').item()
        plt.plot([di,di],[y0+.1,y1-.1],':',color=areacol[ia],linewidth=.5,zorder=1)
        df=di+indat.sel(quantile=0.75,mdl=obs,area=aa,nstats='Length').item()
        plt.plot([df,df],[y0+.1,y1-.1],':',color=areacol[ia],linewidth=.5,zorder=1)

        
        for im,mm in enumerate(mdl):
            ddi=mdn.sel(mdl=mm,area=aa,nstats='FirstDay').item()
            if(np.isfinite(ddi)):
                ddf=ddi+mdn.sel(mdl=mm,area=aa,nstats='Length').item()
                ddons=mdn.sel(mdl=mm,area=aa,nstats='Onset').item()

                yyi=y0+im+.1; yyf=yyi+0.8
                if ('FUT' in mm):  #Fut scen 
                    plt.fill_between([ddi,ddf],[yyi,yyi],[yyf,yyf],edgecolor=mdlcol[im],\
                                    linewidth=1,facecolor='white',alpha=.75)
                else:
                    plt.fill_between([ddi,ddf],[yyi,yyi],[yyf,yyf],edgecolor=None,\
                                    linewidth=.1,facecolor=mdlcol[im])
                plt.plot(np.repeat(ddons,2),[yyi,yyf],color='k',\
                         linewidth=1)
                #plotting interquartile range as whiskers
                yym=y0+im+.5
                dil=indat.sel(mdl=mm,area=aa,nstats='FirstDay',quantile=0.25)
                dih=indat.sel(mdl=mm,area=aa,nstats='FirstDay',quantile=0.75)
                plt.plot([dil,dih],[yym,yym],color='k',linewidth=.5)
                dfl=dil+indat.sel(mdl=mm,area=aa,nstats='Length',quantile=0.25)
                dfh=dih+indat.sel(mdl=mm,area=aa,nstats='Length',quantile=0.75)
                plt.plot([dfl,dfh],[yym,yym],color='k',linewidth=.5)

        dons=mdn.sel(mdl=obs,area=aa,nstats='Onset').item()
        plt.plot(np.repeat(dons,2),[y0,y1],':',color=areacol[ia],\
                 linewidth=2,zorder=1)

    plt.xlim(mn,mx)
    dd=np.asarray([cftime.num2date(ii,calendar='noleap',units='days since 1970-01-01').day\
                  for ii in np.arange(mn,mx+1)])
    mm=np.asarray([cftime.num2date(ii,calendar='noleap',units='days since 1970-01-01').month\
                   for ii in np.arange(mn,mx+1)])
    if(len(dd) < 75):
        pst=((dd == 1) | (dd == 15)).nonzero()[0]
    else: pst=(dd == 1).nonzero()[0]
 
    dtxt=[str(dd[ii])+'/'+str(mm[ii]) for ii in pst]
    plt.xticks(mn+pst,dtxt,fontsize=7,ha='center')

    plt.ylim(0,y1)
    plt.yticks(np.linspace(len(mdl)/2.,y1-len(mdl)/2.,indat.area.size),indat.area.values,\
               fontsize=8,va='center',rotation=90)
    
    ax1.tick_params(axis='both',which='major',length=1,pad=2.4)  

    #caption
    if (len(mdl) <= 4): ax=fig.add_axes([0.08,.02,.9,0.08]); nrow=2
    elif(len(mdl) <= 6): ax=fig.add_axes([0.08,.02,.9,0.12]); nrow=3
    else: ax=fig.add_axes([0.08,.02,.9,0.16]); nrow=4
 
    plt.xlim(0,1); plt.ylim(0,1)  
   
    dy=1./nrow
    px=[0.01,.4]

    ri=0; ci=0
    for ib,b in enumerate(mdl):
        if ('FUT' in b):  #Fut scen
            plt.fill_between([px[ci],.1+px[ci]],[1.-dy*(ri+1)+.02,1.-dy*(ri+1)+.02],\
                             [1.-(dy*ri)-.02,1.-(dy*ri)-.02],facecolor='white',\
                              edgecolor=mdlcol[ib],linewidth=1,zorder=10)
        else:
            plt.fill_between([px[ci],.1+px[ci]],[1.-dy*(ri+1)+.02,1.-dy*(ri+1)+.02],\
                             [1.-(dy*ri)-.02,1.-(dy*ri)-.02],facecolor=mdlcol[ib],\
                             edgecolor=None,linewidth=0,alpha=.75,zorder=10)

        plt.text(0.11+px[ci],1-((ri+.5)*dy),b,va='center',fontsize=6)       
        ri+=1
        if(ri == nrow): ri=0; ci+=1

    plt.axis('off')

    outfile=outf+'.png'
    print(outfile)
    plt.savefig(outfile,format='png',dpi=300)
    plt.close()


def Hovmoller(pt1,pt2,indat,dims,outdim,circular=False):
    """
    Function to create a Hovmoller plot along a transect line
    pt1, pt2: coordinates (lon,lat) of the transect line, in the direction you wish to plot
    indat: input data (DataArray) to be averaged along the transect line
    dims: dimensions (in degrees) of the area around each point of the transect line. 
        if len(dims) == 1: use a circle to average around each point
        if len(dims) == 2: use a rectangle perpendicular to the transect line to average. 
                            In this case, the first value is the width of the rectangle along the transect
                            and the second value is the half-distance perpendicular to the transect
    outdim: dimension (lon or lat) of the transect line
    circular: if False, uses regular averaging over the transect area. else, should pass the 
                number of days per year in the dataset to calculate the circular median.
    Returns a DataArray with the Hovmoller data (points along the transect line versus time) and a mask 
    with the ppoints used
    """
    
    npts=np.ceil(np.sqrt((pt2[0]-pt1[0])**2+(pt2[1]-pt1[1])**2)/\
                 (indat.lon[1]-indat.lon[0]).values).astype(int)
    #number of points along the transect line, considering the same resolution as the input data (lon)
    
    lats=np.asarray(np.linspace(pt1[1],pt2[1],npts))
    lons=np.asarray(np.linspace(pt1[0],pt2[0],npts))
    lat_grid, lon_grid = xr.broadcast(indat.lat, indat.lon)
    
    mask_list = []

    match len(dims):
        case 1:  #circular area around each point
            width_deg = dims[0]
            for lat_c, lon_c in zip(lats, lons):
                mask = ((lat_grid - lat_c)**2 + (lon_grid - lon_c)**2)**0.5 < width_deg
                mask_list.append(mask)
                #Each element in the mask_list corresponds to a circle of radius width_deg 
                #centered at the lat_c/lon_c points along the transect line

        case 2:  #rectangular area perpendicular to the transect line
            w,h=dims
            dx = np.diff(lons); dy = np.diff(lats)
            angles = np.arctan2(dy, dx)  # angle of each segment in radians

            # Extend to same length as transect (by repeating last angle)
            angles = np.append(angles, angles[-1])

            for lat_c, lon_c, theta in zip(lats, lons, angles):
                # Translate grid to box center
                dlat = lat_grid - lat_c
                dlon = lon_grid - lon_c

                # Rotate coordinates (rotate the grid around the center)
                x_rot = dlon * np.cos(theta) + dlat * np.sin(theta)
                y_rot = -dlon * np.sin(theta) + dlat * np.cos(theta)

                # Select points inside the rotated box
                mask = (np.abs(x_rot) <= w / 2) & (np.abs(y_rot) <= h / 2)
                mask_list.append(mask)

        case _:
            print('Error: dims must be a list with 1 or 2 elements')
            sys.exit()

    if(circular):
        averaged=[]
        for mask in mask_list:
            tmp=xr.where(mask == 1,indat,np.nan)
            tmpy=np.asarray([yy for yy in mask.lat.values if mask.sel(lat=yy).any()])
            tmpx=np.asarray([xx for xx in mask.lon.values if mask.sel(lon=xx).any()])
            crped=indat.sel(lon=tmpx,lat=tmpy).where(mask.sel(lon=tmpx,lat=tmpy))
            averaged.append(optb.CircularMedian(crped,circular,dim=["lat", "lon"]))
    else:
        averaged = [indat.where(mask).mean(dim=["lat", "lon"]) for mask in mask_list]
    
    outd=lons if outdim == 'lon' else lats
    hovmoller = xr.concat(averaged, dim="distance_idx")
    hovmoller = hovmoller.assign_coords(distance_idx=outd)
    hovmoller=hovmoller.rename({'distance_idx': outdim})


    #This part provides a mask so it is easier to visualize the area
    msk=np.zeros_like(lat_grid)
    for mask in mask_list: msk=msk+mask
    msk=np.where(msk > 0,1,0)

    return hovmoller, msk


