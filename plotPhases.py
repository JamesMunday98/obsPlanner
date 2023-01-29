#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 19:37:43 2022

@author: james
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime
import pandas as pd
import time
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import AltAz, EarthLocation, SkyCoord


def getLightTravelTimes(ra, dec, time_to_correct):
    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = time_to_correct.light_travel_time(target)
    ltt_helio = time_to_correct.light_travel_time(target, 'heliocentric')
    return ltt_bary, ltt_helio

def output_time(input_format, output_format, location, ra, dec, tinp):
    # based on James McCormac's https://github.com/WarwickAstro/time-conversions
    if not input_format == "mjd":
        # set up the astropy time inputs and convert them to JD-UTC-MID
        if input_format == 'jd':
            #print('Input times in JD, applying no initial correction')
            time_inp = Time(tinp, format='jd', scale='utc', location=location)
        elif input_format == 'hjd':
            #print('Input times in HJD, removing heliocentric correction')
            time_inp = Time(tinp, format='jd', scale='utc', location=location)
            _, ltt_helio = getLightTravelTimes(ra, dec, time_inp)
            time_inp = Time(time_inp.utc - ltt_helio, format='jd', scale='utc', location=location)
        elif input_format == 'bjd':
            #print('Input times in BJD, removing barycentric correction')
            time_inp = Time(tinp, format='jd', scale='tdb', location=location)
            ltt_bary, _ = getLightTravelTimes(ra, dec, time_inp)
            time_inp = Time(time_inp.tdb - ltt_bary, format='jd', scale='tdb', location=location).utc
        else:
            print('Unknown input time format, exiting...')
            raise ValueError
    
        # now convert to the output format requested
        if output_format == 'jd':
            #print('Output set to JD_UTC_MID, no further correction required')
            new_time = time_inp.jd
        elif output_format == 'mjd':
            #print('Output set to MJD_UTC_MID, correcting JD --> MJD')
            new_time = time_inp.mjd
        elif output_format == 'hjd':
            #print('Output set to HJD_UTC_MID, adding heliocentric correction')
            _, ltt_helio = getLightTravelTimes(ra, dec, time_inp)
            new_time = (time_inp + ltt_helio).value
        elif output_format == 'bjd':
            #print('Output set to BJD_TDB_MID, adding barycentric correction')
            ltt_bary, _ = getLightTravelTimes(ra, dec, time_inp)
            new_time = (time_inp.tdb + ltt_bary).value
        else:
            print('Unknown output time format, exiting...')

        return new_time
    
    else: return tinp



cwd = os.getcwd()



TNT=EarthLocation(lat=18.57387*u.deg, lon=98.48656*u.deg, height=2457*u.m)

while True:
    fig, ax=plt.subplots()
    
    # define start/end of night, astro dark in UTC
    dateobs_start=[2023,1,31,12,0,0] # make sure that the last 2 values are 0 here (not the same as sunset)
    sunset=[2023,1,31,11,25,0]
    astrodarkstart=[2023,1,31,12,33,0]
    astrodarkend=[2023,1,31,22,44,0]
    sunrise=[2023,1,31,23,52,0]
    midnight = Time('2023-01-31 17:00:00',scale="utc")
    observatory = TNT
    
    
    for dt in [sunset,astrodarkstart,astrodarkend,sunrise,dateobs_start]:
        inp_date = datetime(dt[0],dt[1], dt[2], dt[3], dt[4], dt[5])
        ts = pd.Timestamp(year = int(inp_date.strftime('%Y')),  month = int(inp_date.strftime('%m')), day = int(inp_date.strftime('%d')), 
                          hour = int(inp_date.strftime('%H')), minute = int(inp_date.strftime('%M')), second = int(inp_date.strftime('%S')))
        
        if dt!=dateobs_start:   ax.axvline(ts.to_julian_date()-2400000.5,c='r',ls='--')
        else:                   MJDdateobs_start=ts.to_julian_date()-2400000.5
    
        if dt==sunset:              MJDsunset=ts.to_julian_date()-2400000.5
        elif dt==astrodarkstart:    MJDastrodarkstart=ts.to_julian_date()-2400000.5
        elif dt==astrodarkend:      MJDastrodarkend=ts.to_julian_date()-2400000.5
        elif dt==sunrise:           MJDsunrise=ts.to_julian_date()-2400000.5
    
    
    
    
    
    
    target=[]
    targets_no_ephem=[]
    f = open("targets.lis")
    
    all_names=[]
    
    # handle all the targets
    for cn,lines in enumerate(f):
        if cn%4==0:      name=lines.split("\n")[0]
        elif cn%4==1:    
            ra_dec = lines
            ra_dec_split = ra_dec.split("\n")[0].split(" ")
            try:
                ra=ra_dec_split[0] + ":" + ra_dec_split[1] + ":" + ra_dec_split[2]
                dec = ra_dec_split[3] + ":" + ra_dec_split[4] + ":" + ra_dec_split[5]
            except:
                ra=ra_dec_split[0]
                dec=ra_dec_split[1]
            
            
            source = SkyCoord(ra,dec, frame='icrs', unit=(u.hourangle, u.deg))
        
            
            delta_midnight = np.linspace(-6, 7, 100)*u.hour
        
        
            frame_July13night = AltAz(obstime=midnight+delta_midnight,  location=observatory)
            sourcealtazs_July13night = source.transform_to(frame_July13night)
        
            sourceairmasss_July13night = sourcealtazs_July13night.secz
            #print(sourcealtazs_July13night.alt)
            
        
            mask_below_airmass_3 = (sourceairmasss_July13night < 2) & (sourcealtazs_July13night.alt > 20*u.deg)
            #plt.plot(delta_midnight[mask_below_airmass_3], sourceairmasss_July13night[mask_below_airmass_3])
            
            times_too_low = delta_midnight[mask_below_airmass_3]
            
            tmin = (midnight + np.amin(times_too_low)).jd -2400000.5
            tmax = (midnight + np.amax(times_too_low)).jd -2400000.5
        
        
        elif cn%4==2:
            ephem = lines
            
            
            inp_format = ephem.split(" ")[0]
            
            if ephem.startswith("H") or ephem.startswith("B"):
                if inp_format.startswith("B"): type_time_format = "bjd"
                elif inp_format.startswith("H"): type_time_format = "hjd"
                else: print("not recognised T0 formating"); raise ValueError
                
                
                
                split=ephem.split(" ")
                order=split[1]
                T=split[2]
                HJD_err=split[3]
                T0=float(T)
                if float(split[2])<100000:
                    T0=T0+2400000.5
                else:
                    None
                
                T0 = output_time(type_time_format, "mjd", observatory, ra, dec, T0)
                
                P=split[4]
                P0err=split[5].split("\n")[0]
                P0=float(P)
                
                target.append([name,order,float(T0),float(P0), tmin, tmax, ra, dec])
                if name.startswith("CSS"):
                    print(T0,P0)
                #print(name,T0, P0)
            elif ephem.startswith("null"):
                targets_no_ephem.append([name,tmin,tmax,ra,dec])
        elif cn%4==3:    None
        
        
        if name.startswith("CSS"):
            print(lines)
        
    
    
    
    # handle all targets with a specific interval set
    target_with_interval=[]
    fil = open("targets.prg")
    for cn,lines_fil in enumerate(fil):
        if cn%4==0:      new_entry=True;    name=lines_fil.split("\n")[0]
        elif cn%4==1:    phase_intervals1 = lines_fil.split("\n")[0]
        elif cn%4==2:
            if lines_fil.startswith("0"):    phase_intervals2 = lines_fil
            else:    phase_intervals2 = "none"
        elif cn%4==3:
            target_with_interval.append([name,phase_intervals1,phase_intervals2])
    
    
    
    
    
    
    print("MJDobs start is ", MJDdateobs_start)
    print()
    
    names_with_intervals = np.asarray(target_with_interval).T[0]
    
    
    # print a staralt copy pasta
    staraltnames=np.asarray(target).T[0]
    for cn, n in enumerate(staraltnames):
        if " " in n:
            staraltnames[cn] = n.split(" ")[0] + n.split(" ")[1]
            n = n.split(" ")[0] + n.split(" ")[1]
        if "+" in n:
            staraltnames[cn] = n.split("+")[0] + n.split("+")[1]
            n = n.split("+")[0] + n.split("+")[1]
        if "-" in n:
            staraltnames[cn] = n.split("-")[0] + n.split("-")[1]
            n = n.split("-")[0] + n.split("-")[1]
        if "." in n:
            staraltnames[cn] = n.split(".")[0]
            n = n.split(".")[0]
        
        
        
    staraltra=np.asarray(target).T[6]
    staraltdec=np.asarray(target).T[7]
    
    np.savetxt("staralt.txt", np.array([staraltnames,staraltra,staraltdec]).T,fmt="%s")
    file=open("staralt.txt")
    for lin in file.readlines():
        print(lin.split("\n")[0])
    
    
    
    for k in target:
        name_tar = k[0]
        tmin = k[4]
        tmax = k[5]
        
        
        #print(tmin, tmax, MJDsunset, MJDsunrise)
        
        if name_tar in names_with_intervals:
            arg = np.argwhere(name_tar==np.asarray(target_with_interval).T[0])[0][0]
            this_interval_set = target_with_interval[arg]
        else:
            this_interval_set=["none","none","none"]
            
        interval1=this_interval_set[1]
        if interval1!="none":
            int1_min = float(interval1.split(" ")[0])
            int1_max = float(interval1.split(" ")[1])
        else:
            int1_min="off"; int1_max="off"
    
        interval2=this_interval_set[2]
        if interval2!="none":
            int2_min = float(interval2.split(" ")[0])
            int2_max = float(interval2.split(" ")[1])
        else:
            int2_min="off"; int2_max="off"
        
        #print(k)
        perTarget = k[3]
        
        
        trial_MJDs = np.linspace(tmin,tmax,43200)
        
        ax.plot([tmin,tmax], [name_tar,name_tar], c='grey', alpha=0.5)
        phase = ((trial_MJDs - k[2])%perTarget)/perTarget
        
        if isinstance(int1_min, float):
            # handle primary eclipse phases
            if int1_min>0:
                print("warning, this code expects the minimum of the interval to be negative, system", name_tar, "may be incorrect")
            to_plot_mask_prim = ((phase-1>int1_min)    |     ((phase<int1_max) & (phase>int1_min)) | ((phase+1<int1_max) & (phase+1>int1_min)))
                
            ax.scatter(trial_MJDs[to_plot_mask_prim],np.full((len(trial_MJDs[to_plot_mask_prim]),),name_tar),c='k',marker='s',s=8)
            if isinstance(int2_min, float):
                # handle primary eclipse phases
                to_plot_mask_sec = (phase<int2_max) & (phase>int2_min)
                ax.scatter(trial_MJDs[to_plot_mask_sec],np.full((len(trial_MJDs[to_plot_mask_sec]),),name_tar),c='b',marker='s',s=8)
            
        #else:
        #    plt.plot([trial_MJDs[0],trial_MJDs[1]],[name_tar,name_tar],c='k',zorder=100)
    
    for k2 in targets_no_ephem:
        a_name,a_tmin,a_tmax,a_ra,a_dec=k2
        trial_MJDs = np.linspace(a_tmin,a_tmax,43200)
        plt.plot(trial_MJDs,np.full((len(trial_MJDs),),a_name), c='b')
        
        
        
    # plot the MJD we currently are on
    ts = pd.Timestamp(year = int(datetime.utcnow().strftime('%Y')),  month = int(datetime.utcnow().strftime('%m')), day = int(datetime.utcnow().strftime('%d')), 
                      hour = int(datetime.utcnow().strftime('%H')), minute = int(datetime.utcnow().strftime('%M')), second = int(datetime.utcnow().strftime('%S'))) 
    
    MJDnow = ts.to_julian_date()-2400000.5
    
    
    if MJDnow>=MJDsunset-0.2:
        ax.axvline(MJDnow)
    
    #interval=2*u.hr
    midnight_hour_for_utc = midnight.value[11:13]
    #print(midnight.jd - 2400000.5)
    
    midnight_mjd = midnight.jd - 2400000.5
    dhour = 1
    interval = dhour/24
    
    xticks_vals = []
    xticks_labels = []
    for i in range(-10,10):
        toappend=midnight_mjd+(i*interval)
        if toappend>MJDsunset-0.5*interval and toappend<MJDsunrise+0.5*interval:
            xticks_vals.append(toappend)
            xticks_labels.append(int(midnight_hour_for_utc)+(i*dhour))
    
    #midnight_mjd-(4*interval),midnight_mjd-(2*interval),midnight_mjd,midnight_mjd+(2*interval), midnight_mjd+(4*interval)
    
    ax.xaxis.set_ticks(xticks_vals)
    ax.set_xticklabels(xticks_labels)
    
    plt.title(str(midnight))
    plt.xlabel('Time (UTC)')
    #plt.xlim(0,1)
    plt.gca().invert_yaxis()
    os.chdir(cwd)
    plt.tight_layout()
    plt.grid(axis='x')
    #plt.savefig('VisitedTargets.png',dpi=400)
    plt.show(block=False)
    plt.pause(10*60)
        
    plt.clf()
    plt.close()


