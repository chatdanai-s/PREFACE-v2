#Part 7 of overall pipeline.
#Also known as Phase Two, transit predictor.
#Take in .csv for each planet and find viable events in a given observing window.
#Now telescope, filter, location and time-dependent.

from __future__ import print_function

#Over-riding astropy_mpl_style makes things pretty!
#import os
import sys
#import timeit

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import ephem

#Spyder imports matplotlib automatically, using an interactive backend (my set default).
#We must switch backends after this to a non-interactive environment - we only need to write plots to directory!
plt.switch_backend('Agg')
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
plt.rcParams.update({'font.size': 18}) #Font size over-ride
from astropy.time import Time
import astropy.units as u
import astropy.coordinates as ac
import scipy.integrate as integrate
from scipy.interpolate import interp1d

from dateutil.relativedelta import relativedelta

csvinpath='../CSV_Bank/P2_CSV_InputParts' #'/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/P2_CSV_InputParts' 	
csvoutpath='../CSV_Bank/P2_CSV_OutputParts' #'/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/P2_CSV_OutputParts' 	

PT8=[30]*1000 #Min obs. altitude for air mass = 2.
LIV=[88]*1000 #Max obs. altitude for the Liverpool instrument
Cross=[-18]*1000 #Marker for night start/end

#Blank arrays for use with time cut/observation ranking
TargetTimes=[]
ErrTimes=[]
Sep=[]
x_hold_fall=[]
x_hold_rise=[]
z_hold=[]
Cross_hold=[]
ErrTot = None

#Print behaviour stored as a variable.
p = print

def get_print():
    return print

#Master function for Phase Two - calculation of event times/metrics.
#def P2Predictor(csv_in, Name, ObsStart, ObsEnd, ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, JobChunk, CoreNo, PlanetNo):
def P2Predictor(csv_in, Name, ObsStart, ObsEnd, ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, PlanetNo):
    #Graphs will live here.
    PlotPath='/mirror/scratch/jmorgan/Codes/Pipeline/P2_Outputs/%s_%s-band_%s_to_%s' % (Inst, Filter, ObsStart.strftime('%b %d %Y'), ObsEnd.strftime('%b %d %Y'))
    
    Cols=['Planet','R*','Rp','RA:HH','RA:MM','RA:SS','Dec:DD','Dec:MM','Dec:SS','T14','Depth','T0 (HJD or BJD)','T0_err',
          'P (day)','P_err', 'Impact Parameter','a_Calc', 'Previous Study Flag', 'Rank', 'Habitable_Rank', 
          'Multi_Transit_Rank','Multi_Transit_Habitable_Rank']
    #csv file containing one planet.
    di=pd.read_csv(csv_in, usecols=Cols, skipinitialspace=True)
    Loc=ac.EarthLocation(lat=(ds['Lat'][S])*u.deg, lon=(ds['Long'][S])*u.deg, height=(ds['Alt'][S])*u.m) #Location info
    #p(Loc)
    
    def AltChecker(row): 
        #Specify co-ordinates of target.
        if row['Dec:DD'] >= 0:
            LockOn = ac.SkyCoord((row['RA:HH']+(row['RA:MM']/60.)+(row['RA:SS']/3600.)), (row['Dec:DD']+(row['Dec:MM']/60.)+(row['Dec:SS']/3600.)),
                     unit=(u.hourangle, u.deg), frame='icrs') 
        else:
            LockOn = ac.SkyCoord((row['RA:HH']+(row['RA:MM']/60.)+(row['RA:SS']/3600.)), (row['Dec:DD']-(row['Dec:MM']/60.)-(row['Dec:SS']/3600.)),
                     unit=(u.hourangle, u.deg), frame='icrs') 
        
        #p(LockOn.dec.degree) 
        
        return LockOn.ra.hour, LockOn.dec.degree
                        
    def T1Calc(row):      
        #Internal conversion to recover midpoint of first observation in BJD_TDB.
        #Scale is not JD, but the format (245XXXX.XXX...) is!
        eph_obj=Time(row['T0 (HJD or BJD)'], format='jd').datetime
        #p(eph_obj)
        
        #Conversion of ephemeris error.
        ESecErr=(row['T0_err']*u.d).to(u.s)
        #p(ESecErr)
        
        #Conversion of error in planet's determined period.
        PSecErr=(row['P_err']*u.d).to(u.s)
        #p(PSecErr)
        
        #Recovers transit start time from published midpoint, returning a timedelta object.
        StartTime=eph_obj-dt.timedelta(0.5*row['T14'])   
        #p(StartTime)
        
        #Builds range of n relevant for target/window combination.
        #nstart finds final transit before our window begins, so you gotta add +1!
        nStart=int((ObsStart-StartTime).total_seconds()/((row['P (day)']*u.d).to(u.s).value))+1 
        nEnd=int((ObsEnd-StartTime).total_seconds()/((row['P (day)']*u.d).to(u.s).value))       
        Ran=np.linspace(nStart,nEnd,((nEnd-nStart)+1))
        #p(Ran)
        
        #Increase the number of orbits n to find future transit start times in our window.
        for n in Ran:
            global ErrTot
            NextTransitTime=(StartTime+(int(n)*dt.timedelta(row['P (day)'])))
            #p(('Transit begins at %s BJD_TDB') % NextTransitTime)
            TargetTimes.append(NextTransitTime)
            #p('Event appended!')
            
            #Calculation of total error.
            #Uncertainty should scale overall as no. of orbits n, so this needs to be inside the squared term!
            ErrTot=np.sqrt(((ESecErr)**2)+(n*PSecErr)**2).value              
            ErrTimes.append(ErrTot)           
            #p(ErrTot)
                
        return ErrTot
    
    #Calculation of key times.
    def TimesCalc(row):
        
        try:
            with np.errstate(invalid='raise'):        
                T23=(row['P (day)']/np.pi)*np.arcsin(np.sqrt((((row['R*']*u.solRad)-((row['Rp']*u.R_jup).to(u.solRad)))**2.)
                -((row['Impact Parameter']*(row['R*']*u.solRad))**2.))/((row['a_Calc']*u.AU).to(u.solRad)))
                #p('T23 = %s' % T23)
                
                #Recovers length of ingress/egress. abs is a necessary guess if T23 > T14 (unphysical soln.)     
                T12=abs((row['T14']-(T23.value))/2.) 
                #p(T12)
                
                #Time at end of ingress
                T2=row['T1']+dt.timedelta(T12)
                #p('T2 at %s BJD_TDB' % T2)
                  
                #Calculates midpoint, end and length of transit event.
                TransitMid=row['T1']+dt.timedelta(0.5*row['T14'])   
                #p('T0 at %s BJD_TDB' % TransitMid)
                
                TransitEnd=row['T1']+dt.timedelta(row['T14'])   
                #p('T4 at %s BJD_TDB' % TransitEnd)
                
                #Start of egress
                T3=TransitEnd-dt.timedelta(T12)
                #p('T3 at %s BJD_TDB' % T3)
                
                #TESS targets have dodgy parameters - this calculation will catch them.
                #Forces a FloatingPointError so the transit is treated as grazing, passed to except block.
                np.sqrt((T3-T2).total_seconds())
                
                #A solid baseline is important for a good transit curve!
                BaseOne=row['T1']-dt.timedelta(0,0,0,0,60)
                #p('For effective baseline, begin observations at %s BJD_TDB' % BaseOne)
                
                BaseTwo=TransitEnd+dt.timedelta(0,0,0,0,60)
                #p('For effective baseline, end observations at %s BJD_TDB' % BaseTwo)
        
        #Grazing transits have no T23, so a different approach is needed.
        except FloatingPointError:
            T23=0.
            #p('Grazing exception occurred!')
                          
            #Calculates midpoint, end and length of transit event.
            TransitMid=row['T1']+dt.timedelta(0.5*row['T14'])   
            #p('T0 at %s BJD_TDB' % TransitMid)
            
            #Time at end of ingress
            T2=TransitMid
            #p('T2 at %s BJD_TDB' % T2)
            
            #Start of egress
            T3=TransitMid
            #p('T3 at %s BJD_TDB' % T3)
            
            TransitEnd=row['T1']+dt.timedelta(row['T14'])   
            #p('T4 at %s BJD_TDB' % TransitEnd)
                        
            #A solid baseline is important for a good transit curve!
            BaseOne=row['T1']-dt.timedelta(0,0,0,0,60)
            #p('For effective baseline, begin observations at %s BJD_TDB' % BaseOne)
            
            BaseTwo=TransitEnd+dt.timedelta(0,0,0,0,60)
            #p('For effective baseline, end observations at %s BJD_TDB' % BaseTwo)
        except ValueError:
            sys.exit('Time error for %s (%s) transiting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
            
        return T2, TransitMid, T3, TransitEnd, BaseOne, BaseTwo
    
    ###WE GET INTO CALCULATION AND PLOTTING HERE###       
    
    def EventMetric(row): 
        #Specify co-ordinates of target, as previously converted.
        LockOn = ac.SkyCoord(row['RA:HMS'], row['Dec:Deg'], unit=(u.hourangle, u.deg), frame='icrs') 
        #p(LockOn)
        
        #Chooses appropriate midnight for each useful time (for graphing purposes).
        #Use relativedelta to avoid end of months/years breaking the code. 
        if row['T1'].hour >= 12:
            MidnightSet=row['T1'].replace(hour=0,minute=0,second=0)+relativedelta(days=1)
        else:
            MidnightSet=row['T1'].replace(hour=0,minute=0,second=0)
        
        midnight=Time(MidnightSet) #Astropy object conversion for compatibility.
        delta_midnight = np.linspace(-12, 18, 1000)*u.hour
        Moon_step = np.linspace(-12,18,21)*u.hour
        times = midnight + delta_midnight #Times for Moon/target calc
        Moon_times = midnight + Moon_step #Times for Moon/target separation calc
        frame = ac.AltAz(obstime=times, location=Loc) #Creates obs frames
        Moon_frame = ac.AltAz(obstime=Moon_times, location=Loc)
        
        #Where is the sun from midday to evening?
        sunaltazs = ac.get_sun(times).transform_to(frame) #AltAz pairs created
        
        #Where is the Moon from midday to evening?
        #VERSION NOTE: ONLY COMPATIBLE WITH ASTROPY VERSION 1.2.1
        moonaltazs = ac.get_moon(times).transform_to(frame) #AltAz pairs created
        moonsepaltazs = ac.get_moon(Moon_times).transform_to(Moon_frame)
        
        #Where is the target from midday to evening?
        Targetaltazs = LockOn.transform_to(frame) #AltAz pairs created
        Targetsepaltazs = LockOn.transform_to(Moon_frame)
        
        #Calculates separation between target and Moon at 21 intervals.
        for e in range(len(Moon_step)):
            S1=ac.SkyCoord(Targetsepaltazs.az,Targetsepaltazs.alt, unit = "deg")
            #p(S1)
            S2=ac.SkyCoord(moonsepaltazs.az,moonsepaltazs.alt, unit = "deg")
            S3=S1.separation(S2)
        
        #Clears separations for previous transit and formats current one.
        if len(Sep) == 21:
            del Sep[:]
        else:
            pass
        
        Sep.extend(S3.value)
        FormSep = [int(elem) for elem in Sep]   
        #p(FormSep)
        
        ###LUNAR ILLUMINATION CALC
        
        m=ephem.Moon()
        m.compute(row['T1'])
        Phase=(m.moon_phase)*100.
        #p('Lunar disk illumination is %.3f %%.' % Phase) #%% prints a precentage sign.
        
        #CLEARANCE PROTOCOLS FOR MULTIPLE-EVENT RUNNING
        
        if len(x_hold_fall) > 0:
            del x_hold_fall[:]
        else:
            pass
        
        if len(x_hold_rise) > 0:
            del x_hold_rise[:]
        else:
            pass
        
        if len(Cross_hold) > 0:
            del Cross_hold[:]
        else:
            pass
        
        #When will the transit, baseline and night begin and end relative to midnight? 
        #Finds differences either side of 00:00:00.
        #NOTE: Take only the value when plotting! You can't plot a unit.
        
        #When is transit midpoint relative to midnight?    
        dif0=((midnight-Time(row['T0'])).to('hr').value)*-1.   
        #p('T14 0 = %.2f hrs relative to midnight. ' % dif0)
            
        difh1=((midnight-Time(row['Baseline_Start'])).to('hr'))*-1.
        dif1=difh1.value
        #p('Baseline start = %.2f hrs relative to midnight.' % dif1)
        
        difh2=((midnight-Time(row['Baseline_End'])).to('hr'))*-1.
        dif2=difh2.value
        #p('Baseline end = %.2f hrs relative to midnight.' % dif2)
        
        difh3=((midnight-Time(row['T1'])).to('hr'))*-1.
        dif3=difh3.value
        #p('T1 = %.2f hrs relative to midnight.' % dif3)
        
        difh4=((midnight-Time(row['T4'])).to('hr'))*-1.
        dif4=difh4.value
        #p('T4 = %.2f hrs relative to midnight.' % dif4)
        
        difh5=((midnight-Time(row['T2'])).to('hr'))*-1.
        dif5=difh5.value
        #p('T2 = %.2f hrs relative to midnight.' % dif5)
        
        difh6=((midnight-Time(row['T3'])).to('hr'))*-1.
        dif6=difh6.value
        #p('T3 = %.2f hrs relative to midnight.' % dif6)
        
        #When does each event night start and end? 
        #Checks when sun path crosses twilight marker at y=-18, appends to x_hold_fall or _rise (depending on gradient).
        CrossCheck = sunaltazs.alt.value - Cross
        #p(CrossCheck)
        for i in range(len(CrossCheck) - 1):
            if CrossCheck[i] == 0. or CrossCheck[i] * CrossCheck[i + 1] < 0.:
                # crossover at i
                x_ = delta_midnight.to('min').value[i]
                #p("Crossing point at %.2f mins" % x_)
                if sunaltazs.alt.value[i+1] - sunaltazs.alt.value[i] < 0.:
                    x_hold_fall.append(x_)
                elif sunaltazs.alt.value[i+1] - sunaltazs.alt.value[i] > 0.:
                    x_hold_rise.append(x_)
        
        #p(x_hold_fall)
        #p(x_hold_rise)
        
        #For "two-night events", with a bit at the start and end!
        #It is the second condition in particular which catches these.
        NGen1=(len(x_hold_fall) == len(x_hold_rise), x_hold_fall[0] > x_hold_rise[0], dif0 < 0.)
        NGen2=(len(x_hold_fall) == len(x_hold_rise), x_hold_fall[0] > x_hold_rise[0], dif0 > 0.)
        #Condition for two falls, one rise, transit in first night.
        NGen3=(len(x_hold_fall) > len(x_hold_rise), dif0 < 0.)
        
        #p(NGen1)
        #p(NGen2)
        #p(NGen3)
        
        if all(NGen1):
            x_hold_fall[0]=-780.
        elif all(NGen2):
            x_hold_rise[0]=1160.
        elif all(NGen3):
            x_hold_rise.append(1160.)
        else:
            pass
        
        #p(x_hold_fall)
        #p(x_hold_rise)
        
        if len(x_hold_rise) > len (x_hold_fall):
            TNight=abs(x_hold_rise[1]-x_hold_fall[0])
        elif len(x_hold_rise) < len (x_hold_fall):
            TNight=abs(x_hold_rise[0]-x_hold_fall[1])  
        else:
            TNight=abs(x_hold_rise[0]-x_hold_fall[0]) #Night length
        
        #p('TNight = %.2f mins' % TNight)
                 
        NGen4=(len(x_hold_fall) > len(x_hold_rise), dif0 > 0.)
        NGen5=(len(x_hold_fall) < len(x_hold_rise))
        NGen6=(len(x_hold_fall) < len(x_hold_rise), dif0 > 0.)
        
        if all(NGen4): #Two falls, one rise, transit during second night.
            x_hold_fall[0]=x_hold_fall[1]
            x_hold_rise[0]=x_hold_fall[1]+TNight
        elif all(NGen6): #Exception for one fall, two rises, transit in second night.
            x_hold_rise[0]=x_hold_rise[1]    
        elif NGen5 == True: #Exception for one fall, two rises, transit in first night.
            FallOne = x_hold_rise[0]-TNight
            x_hold_fall.insert(0, FallOne)
        else:
            pass
        
        #p(x_hold_fall)
        #p(x_hold_rise)

        TNight0=x_hold_rise[0]-TNight/2. #Night midpoint
        #p('TNight0 = %.2f mins relative to midnight.' % TNight0)
                               
        #HEIGHT AND NIGHT SPOT CHECKS   
        
        #Generates time-space for overlap between transit and night.
        Ett=[difh3.value,(x_hold_fall[0]/60.)]
        Tvo=[difh4.value,(x_hold_rise[0]/60.)]
        Tre=midnight + np.linspace(max(Ett),min(Tvo),1000)*u.hr
        
        TDur = np.linspace(0, ((row['T14']*1440.)+120.), 1000)*u.min #Generates transit time-space
        #p(TDur)
        Ttimes = Time(row['Baseline_Start']) + TDur #Iterates across time-space
        #p(Ttimes)
        Tframe = ac.AltAz(obstime=Ttimes, location=Loc) #Creates new obs frames  
        SC = LockOn.transform_to(Tframe) #New AltAz pairs created
        #p(SC.alt.value[999])
                
        Nt=np.linspace(x_hold_fall[0],x_hold_rise[0],1000)*u.min #Generates night time-space
        #p(Nt)
        Ntimes = midnight + Nt #Iterates across time-space
        #p(Ntimes)
        Nframe = ac.AltAz(obstime=Ntimes, location=Loc) #Creates new obs frames  
        NC = LockOn.transform_to(Nframe) #New AltAz pairs created
        #p(NC.alt.value[999])
        
        Kframe = ac.AltAz(obstime=Tre, location=Loc) #Creates new obs frames  
        KC = LockOn.transform_to(Kframe)
        #p(KC.alt.value)
        
        HSCGen=((a >= PT8[0] for a in SC.alt.value), (a >= PT8[0] for a in NC.alt.value))
        HSC=any(HSCGen) #Is the target at observable height anywhere in the run?
        HFC=all(a >= PT8[0] for a in SC.alt.value) #Is the target at observable height throughout the run?
                                      
        NSC=any(a >= PT8[0] for a in NC.alt.value) #Is the target at observable height at some point in the night?
        
        KSC=any(a >= PT8[0] for a in KC.alt.value) #Is the target observable at height, at night? (K2-99 exception.)
        #p(KC.alt.value)
        
        #IMPLEMENT CONDITIONS 1-3        
        #Bank of generator expressions (S denotes "strict", for which target must always be sufficiently high)
        C1Gen=((abs((dif0*60.)-TNight0) <= (((row['T14']*1440.)+TNight)/2.)),HSC, KSC)
        C1GenS=((abs((dif0*60.)-TNight0) <= (((row['T14']*1440.)+TNight)/2.)),HFC,NSC)
        C2Gen=((((row['T14']*1440.) + (2.*abs((dif0*60.)-TNight0)))<= TNight),HSC, KSC)
        C2GenS=((((row['T14']*1440.) + (2.*abs((dif0*60.)-TNight0)))<= TNight),HFC,NSC)
        C3Gen=(((((TNight-row['T14']*1440.)/2.)-abs((dif0*60.)-TNight0))>= 60.),HSC, KSC)
        C3GenS=(((((TNight-row['T14']*1440.)/2.)-abs((dif0*60.)-TNight0))>= 60.),HFC,NSC)

        #Assigns internal rank. If one condition is untrue, an all statement will return 'False'.
        #Events which fail an all statement are passed down to the next statement until 'True' is returned and they stop.
        if all(C3GenS):
            Internal_Rank='03_F'
        elif all(C3Gen):
            Internal_Rank='03_P'
        elif all(C2GenS):
            Internal_Rank='02_F'
        elif all(C2Gen):
            Internal_Rank='02_P'
        elif all(C1GenS):
            Internal_Rank='01_F'
        elif all(C1Gen):
            Internal_Rank='01_P'
        else:
            Internal_Rank='X'   
        #p(Internal_Rank)
            
        #Choose integral limits and return key times.
        #We have one integration to do, for transit + baselines for air mass marker.
        
        #If the target is always above 30 degrees, there are no crossing points to worry about!
        
        if 'F' in Internal_Rank:
            global L1, L2
            try:
                #FOR TRANSIT + BASELINE
                if difh1.value < (x_hold_fall[0])/60.:
                    L1 = (x_hold_fall[0])/60.
                else:
                    L1 = difh1.value
                    
                #p("Lower air mass integral limit is %.2f hrs relative to midnight." % L1)
               
                if difh2.value < (x_hold_rise[0])/60.:
                    L2 = difh2.value               
                else:
                    L2 = (x_hold_rise[0])/60.
                         
                #Get-out for super-marginal events (eg HAT-P-66).                    
                if L2 < L1:
                    Internal_Rank='X'
                else:
                    pass
                       
                #p("Upper air mass integral limit is %.2f hrs relative to midnight." % L2) 
                
            except NameError:
                sys.exit('%s (%s) has thrown an exception for a full event at %s - look at me!' % (Name, dj['Planet'][0], row['T1']))
                
        #If there are crossing points, passed to bank of generator expressions.
        #any() will stop iterating as soon as it returns a "True" case.
        #For lower integral limit:
        elif 'P' in Internal_Rank:    
            try:
                PGen=(SC.alt.value[0] < PT8[0], SC.alt.value[999] < PT8[0])   
                if any(PGen):
                    PT8CC = SC.alt.value - PT8[0]
                    for i in range(len(PT8CC) - 1):
                        if PT8CC[i] == 0. or PT8CC[i] * PT8CC[i + 1] < 0.:
                            LT1 = ((TDur.to('hr').value[i])-1.) #relative to start of transit
                            #p("Crossing point is %.2f hrs relative to start of transit." % LT1)
                            LTC1=difh3.value+LT1 #convert to be relative to midnight.
                            #p("Crossing point is %.6f hrs relative to midnight." % LTC1)
                            Cross_hold.append(LTC1)
                            #p(Cross_hold)
                            #p(len(Cross_hold))
                            
                #For targets that are barely visible, we get two crossing points!
                #These observing limits are easier to calculate.
                if len(Cross_hold) == 2 :
                    
                    if Cross_hold[0] < (x_hold_fall[0])/60.:
                        L1=(x_hold_fall[0])/60
                    else:
                        L1=Cross_hold[0]
                    
                    #p("Lower air mass integral limit (double-cross) is %.6f hrs relative to midnight." % L1)   
                    
                    if (x_hold_rise[0])/60. < Cross_hold[1]:
                        L2=(x_hold_rise[0])/60
                    else:
                        L2=Cross_hold[1]
                                   
                    #p("Upper air mass integral limit (double-cross) is %.6f hrs relative to midnight." % L2)
                
                #One crossing point, calculate observing limits:
                elif len(Cross_hold) == 1 :
                    PGen2=(difh1.value < (x_hold_fall[0])/60., SC.alt.value[999] < PT8[0])
                    #LTC1 < difh2.value < (x_hold_rise[0])/60., 
                    #p(PGen2)
                    PGen12=(difh3.value < (x_hold_fall[0])/60. < LTC1 < difh4.value < (x_hold_rise[0])/60., SC.alt.value[999] < PT8[0])         
                    #p(PGen12) #For transits that are too long on fall side for one night.
                    PGen15=(difh3.value < (x_hold_fall[0])/60. < LTC1 < (x_hold_rise[0])/60. < difh4.value, SC.alt.value[999] < PT8[0])
                    #p(PGen15) #This condition is for transits that are too long for 1 night!
                    PGen3=(all(PGen2),all(PGen12), all(PGen15), (x_hold_fall[0])/60. > LTC1)
                    PGen4=((x_hold_fall[0])/60. < LTC1 < difh2.value, (x_hold_fall[0])/60. > difh1.value)
                    PGen5=(difh1.value < (x_hold_rise[0])/60. < LTC1, difh1.value < LTC1 < (x_hold_rise[0])/60.)
                    PGen6=((x_hold_fall[0])/60. < LTC1 < difh2.value, (x_hold_rise[0])/60. > difh1.value, SC.alt.value[999] > PT8[0])
                    PGen7=(difh1.value < (x_hold_fall[0])/60. < LTC1 < difh2.value, SC.alt.value[999] > PT8[0])
                    PGen14=(all(PGen4),all(PGen6), all(PGen7))
                    #any () & all() only accept one argument, so it is necessary to nest generator expressions.
                    #Pass denotes that no action should be taken. != means 'is not equal to'.
                    if any(PGen3):
                        L1 = (x_hold_fall[0])/60.
                    elif any(PGen14):
                        L1 = LTC1
                    elif any(PGen5):
                        L1 = difh1.value   
                    else:
                        sys.exit('Bad lower integration limit for %s (%s) transting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
                
                    #p("Lower air mass integral limit is %.6f hrs relative to midnight." % L1)   
                    
                    #############################################################################
    
                    PGen8=(round(difh1.value,6) <= round(LTC1, 6) < difh2.value < x_hold_rise[0]/60., SC.alt.value[999] > PT8[0])      
                    PGen9=(LTC1 < (x_hold_rise[0])/60. < difh2.value, LTC1 < difh2.value)   
                    PGen10=(round(difh1.value,6) <= round(LTC1, 6) < (x_hold_rise[0])/60. < difh2.value, SC.alt.value[999] > PT8[0])
                    PGen13=((x_hold_fall[0])/60. < difh3.value < LTC1 < (x_hold_rise[0])/60. < difh4.value, SC.alt.value[999] > PT8[0])                               
                    #p(PGen13) #For transits that are too long on rise side for one night.
                    PGen11=(all(PGen10), all(PGen13), LTC1 > (x_hold_rise[0])/60.)
    
                    if all(PGen8): 
                        L2 = difh2.value 
                    elif any(PGen11):
                        L2 = (x_hold_rise[0])/60.
                    elif any(PGen9):
                        L2 = LTC1
                    else:
                        sys.exit('Bad upper integration limit for %s (%s) transting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
               
                   #p("Upper air mass integral limit is %.6f hrs relative to midnight." % L2)
            
                #Get-out for super-marginal events (eg HAT-P-66).                
                if L2 < L1:
                    Internal_Rank='X'
                else:
                    pass
                
            except NameError:
                sys.exit('%s (%s) has thrown an exception for a partial event at %s - look at me!' % (Name, dj['Planet'][0], row['T1']))
                
        else:
            pass #This will throw an UnboundLocalError if an event has been incorrectly flagged.     
                
        #Choose all segment limits for event visibility metric.
        if Internal_Rank != 'X':
            
            #First baseline visibility limits       
            PGen1 = (L2 > difh1.value >= L1, L1 < difh3.value < L2) #Both start and end visible
            #p(PGen1)       
            PGen2 = (L2 > difh1.value > L1, difh3.value >= L2) #Only start visible
            #p(PGen2)
            PGen3 = (difh1.value <= L1, L1 < difh3.value < L2) #Only end visible
            #p(PGen3)
            PGen7 = (any(PGen2),any(PGen3))
            PGen4 = (difh1.value < L1, difh3.value < L1) #None visible
            #p(PGen4)
            PGen5 = (difh1.value > L2, difh3.value > L2)
            #p(PGen5)
            PGen6 = (all(PGen4),all(PGen5))
            #p(PGen6)
            
            if all(PGen1):
                L3=difh1.value
                L4=difh3.value
            elif all(PGen2):
                L3=difh1.value
                L4=L2
            elif all(PGen3):
                L3=L1
                L4=difh3.value
            elif all(PGen7):
                L3=L1
                L4=L2
            elif any(PGen6):
                L3=0.0
                L4=0.0
            else: #Exception to fire if something goes wrong...
                sys.exit('Bad first baseline limits for %s (%s) transiting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
                
            #p("Lower first baseline limit is %.2f hrs relative to midnight." % L3)
            #p("Upper first baseline limit is %.2f hrs relative to midnight." % L4)
            
            #Ingress visibility limits
            
            PGen1 = (L2 >= difh3.value > L1, L1 < difh5.value < L2)
            #p(PGen1)
            PGen2 = (L2 >= difh3.value > L1, difh5.value > L2)
            #p(PGen2)
            PGen3 = (difh3.value < L1, L1 < difh5.value < L2)
            #p(PGen3)
            PGen7 = (any(PGen2),any(PGen3))
            PGen4 = (difh3.value < L1, difh5.value < L1)
            #p(PGen4)
            PGen5 = (difh3.value >= L2, difh5.value > L2)
            #p(PGen5)
            PGen6 = (all(PGen4),all(PGen5))
            
            if all(PGen1):
                L5=difh3.value
                L6=difh5.value
            elif all(PGen2):
                L5=difh3.value
                L6=L2
            elif all(PGen3):
                L5=L1
                L6=difh5.value
            elif all(PGen7):
                L5=L1
                L6=L2
            elif any(PGen6):
                L5=0.0
                L6=0.0
            else:
                sys.exit('Bad ingress limits for %s (%s) transiting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
                
            #p("Lower T12 limit is %.2f hrs relative to midnight." % L5)
            #p("Upper T12 limit is %.2f hrs relative to midnight." % L6)
            
            #Full transit limits
            
            PGen1 = (L2 > difh5.value > L1, L1 < difh6.value < L2)
            #p(PGen1)
            PGen2 = (L2 > difh5.value > L1, difh6.value > L2)
            #p(PGen2)
            PGen3 = (difh5.value < L1, L1 < difh6.value < L2)
            #p(PGen3)        
            PGen4 = (difh5.value < L1 < L2 < difh6.value) #Crossing point occurs
            #p(PGen4)
            PGen5 = (difh5.value < L1, difh6.value < L1)
            #p(PGen5)
            PGen6 = (difh5.value > L2, difh6.value > L2)
            #p(PGen6)
            PGen7 = (all(PGen5),all(PGen6))
            #p(PGen7)
            
            if all(PGen1):
                L7=difh5.value
                L8=difh6.value
            elif all(PGen2):
                L7=difh5.value
                L8=L2
            elif all(PGen3):
                L7=L1
                L8=difh6.value
            elif PGen4 == True:
                L7=L1
                L8=L2
            elif any(PGen7):
                L7=0.0
                L8=0.0
            else:
                sys.exit('Bad full transit limits for %s (%s) transiting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
                
            #p("Lower T23 limit is %.2f hrs relative to midnight." % L7)
            #p("Upper T23 limit is %.2f hrs relative to midnight." % L8)
            
            #Egress limits
            
            PGen1 = (L2 > difh6.value > L1, L1 < difh4.value < L2)
            #p(PGen1)
            PGen2 = (L2 > difh6.value > L1, difh4.value > L2)
            #p(PGen2)
            PGen3 = (difh6.value < L1, L1 < difh4.value < L2)
            #p(PGen3)
            PGen7 = (any(PGen2),any(PGen3))
            PGen4 = (difh6.value < L1, difh4.value < L1)
            #p(PGen4)
            PGen5 = (difh6.value > L2, difh4.value > L2)
            #p(PGen5)
            PGen6 = (all(PGen4),all(PGen5))
            #p(PGen6)
            
            if all(PGen1):
                L9=difh6.value
                L10=difh4.value
            elif all(PGen2):
                L9=difh6.value
                L10=L2
            elif all(PGen3):
                L9=L1
                L10=difh4.value
            elif all(PGen7):
                L9=L1
                L10=L2
            elif any(PGen6):
                L9=0.0
                L10=0.0
            else:
                sys.exit('Bad egress limits for %s (%s) transiting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
                
            #p("Lower T34 limit is %.2f hrs relative to midnight." % L9)
            #p("Upper T34 limit is %.2f hrs relative to midnight." % L10)
            
            #Second baseline limits
            
            PGen1 = (L2 > difh4.value > L1, L1 < difh2.value <= L2)
            #p(PGen1)
            PGen2 = (L2 > difh4.value > L1, difh2.value > L2)
            #p(PGen2)
            PGen3 = (difh4.value < L1, L1 < difh2.value < L2)
            #p(PGen3)
            PGen4 = (difh4.value < L1, difh2.value < L1)
            #p(PGen4)
            PGen5 = (difh4.value > L2, difh2.value > L2)
            #p(PGen5)
            PGen6 = (all(PGen4),all(PGen5))
            #p(PGen6)
            
            if all(PGen1):
                L11=difh4.value
                L12=difh2.value
            elif all(PGen2):
                L11=difh4.value
                L12=L2
            elif all(PGen3):
                L11=L1
                L12=difh2.value
            elif any(PGen6):
                L11=0.0
                L12=0.0
            else:
                sys.exit('Bad second baseline limits for %s (%s) transiting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))
                
            #p("Lower second baseline limit is %.2f hrs relative to midnight." % L11)
            #p("Upper second baseline limit is %.2f hrs relative to midnight." % L12)
            
            #What times should you put on your observing proposal? For this, UTC is convenient.
            #All previous calculations use BJD_TDB - these are best for publications!
            #This step was previously handled by P2_Timesplitter, but has been shrunk to go here.
            PropStart = midnight.datetime+dt.timedelta(0,0,0,0,0,L1)
            ObsStart_BJD = Time(PropStart, scale='tdb', location=Loc) #BJD_TDB (implicit)
            ltt_bary_Start = ObsStart_BJD.light_travel_time(LockOn) #JD to BJD conversion factor.
            ObsStart = Time(PropStart, scale='tdb')-ltt_bary_Start #JD_TDB
            ObsStart_JD = ObsStart.utc
            
            PropEnd = midnight.datetime+dt.timedelta(0,0,0,0,0,L2)
            ObsEnd_BJD = Time(PropEnd, scale='tdb', location=Loc) #BJD_TDB (implicit)
            ltt_bary_End = ObsEnd_BJD.light_travel_time(LockOn) #JD to BJD conversion factor.
            ObsEnd = Time(PropEnd, scale='tdb')-ltt_bary_End #JD_TDB
            ObsEnd_JD = ObsEnd.utc
            
            #Recovery of base for air mass integration.
            #If it's a two-night event, a slightly extended integration base is needed. Try/except catches this.
            try:
                with np.errstate(invalid='raise'):
                    f3 = interp1d(delta_midnight, abs(np.sin(((Targetaltazs.alt).degree)*(np.pi/180)))**(-0.6), kind='cubic')
                    Sval, Serr = integrate.quad(f3,L1,L2)
                    Air_Mass_Marker = max(0.000,3.*(((L2-L1)/Sval)-(2./3.))) #Returns an air mass weighting between 0 and 1, where 30 degrees = 0.
                    #p(("Air mass metric is %.3f \n") % Air_Mass_Marker)    
            except ValueError:
                ext_delta_midnight=np.linspace(-13,19,1000)*u.hour
                f3 = interp1d(ext_delta_midnight, abs(np.sin(((Targetaltazs.alt).degree)*(np.pi/180)))**(-0.6), kind='cubic')
            		
            try:
                with np.errstate(invalid='raise'): 
                    Sval, Serr = integrate.quad(f3,L1,L2)
                    Air_Mass_Marker = max(0.000,3.*(((L2-L1)/Sval)-(2./3.))) #Returns an air mass weighting between 0 and 1, where 30 degrees = 0.
                    #p(("Air mass metric is %.3f \n") % Air_Mass_Marker) 
            except ValueError: #Second exception for very long transiters (eg. Kepler-432)
                Sval, Serr = integrate.quad(f3,L1-24,L2-24)
                Air_Mass_Marker = max(0.000,3.*((((L2-24)-(L1-24))/Sval)-(2./3.)))
                #p(("Long air mass metric is %.3f \n") % Air_Mass_Marker)
                
            #Weight by percentage of event + baseline captured.    
            try:
                with np.errstate(invalid='raise'):
                    WBase=np.around((((L4-L3)+(L12-L11))/2.),4) #Baseline
                    #p('Baseline weight = %s' % WBase)
                    WTrans=np.around(((L8-L7)/(difh6.value-difh5.value)),4) #Full depth
                    #p('Full transit weight = %s' % WTrans)
                    WInOut=np.around((((L10-L9)+(L6-L5))/(2*(difh5.value-difh3.value))),4) #Ingress-egress
                    #p('Ingress-egress weight = %s' % WInOut)
                           
                    EventPercent=np.around((WBase*WTrans*WInOut),4)
                    #p(("Event weight is %.4f \n") % EventPercent)                    
            #except ZeroDivisionError:
            except FloatingPointError:
                WBase=np.around((((L4-L3)+(L12-L11))/2.),4) #Baseline
                WTrans=np.nan #Full depth (never achieved for grazing events!)
                WInOut=np.around((((L10-L9)+(L6-L5))/(2*(difh5.value-difh3.value))),4) #Ingress-egress
                    
                EventPercent=np.around((WBase*WInOut),4)
                #p(("Grazing event weight is %.4f \n") % EventPercent)
            #except FloatingPointError:
            except ZeroDivisionError:   
                sys.exit('ZeroDivisionError for %s (%s) transiting at %s - flag me!' % (Name, dj['Planet'][0], row['T1']))

        else:
            ObsStart_JD=np.nan
            ObsEnd_JD=np.nan
            Air_Mass_Marker = 0.000 #Unobservable 
            WBase=0.0
            WTrans=0.0
            WInOut=0.0
            EventPercent = 0.00
        
        # #Plot target path and transit viability for all events.
        # #Scatter c argument denotes colour, lw=linewidth, default s=20, zorder stops fill from mucking up the colour of your other lines.
        # #Use plt.close(fig) to stop all the figures from eating your RAM.
        
        plt.figure(figsize=(18, 10))
        #plt.plot(delta_midnight, sunaltazs.alt, color='r', label='Sun') #Solar path
        plt.plot(delta_midnight, moonaltazs.alt, 'g--', label='Lunar Path') #Moon path
        # plt.plot(Moon_step, moonsepaltazs.alt, 'go', label='Moon Step') #Moon path
        plt.plot(delta_midnight, PT8 , label='Min. Altitude')
        plt.plot(delta_midnight, LIV , label='Max. Altitude')    
        plt.scatter(delta_midnight, Targetaltazs.alt, c = Targetaltazs.az.value, label = 'Target Path', lw=0, s=8, cmap='winter') #Target path
        plt.plot(Moon_step, Targetsepaltazs.alt, 'bo', label='Lunar Separation (deg)')
        # #Attaches separations to markers along target path.
        for w, txt in enumerate(FormSep):
            plt.annotate(FormSep[w], (Moon_step.value[w], Targetsepaltazs.alt.value[w]), color = 'orange', xytext=(-7.5, -22), textcoords='offset pixels')   
        # #plt.fill_between(delta_midnight.to('hr').value, 0, 90, sunaltazs.alt <= 0*u.deg, color='navy', zorder=0, label='Twilight') #Bounded by astronomical twilight
        # #Dark fill starts once sun reaches -18 deg alt, superimposed on twilight zone.
        plt.fill_between(delta_midnight.to('hr').value, 0, 90, sunaltazs.alt <= -18*u.deg, color='black', zorder=0, label='Night') #Bounded by astronomical twilight
        # #Red fill denotes transit in progress, superimposed on other zones. alpha denotes transparency.
        plt.axvspan(difh1.value,difh2.value, 0,90, alpha=0.3, color='r', zorder=0, label='Observation')
        plt.axvspan(difh3.value,difh4.value, 0,90, alpha=0.3, color='r', hatch='X', zorder=0, label='Transit Active')
        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper right')
        plt.text(23,85,'Moon Illumination:\n %.0f %%' % Phase)
        plt.text(23,75,'Closest Lunar\n Separation: %.0f deg' % min(FormSep))
        plt.xlim(-12, 18)
        plt.xticks(np.arange(15)*2 -12)
        plt.ylim(0, 90)
        plt.xlabel('Hours from BJD_TDB Midnight')
        plt.ylabel('Altitude [deg]')
        plt.title(r'%s Transit Observability for %s BJD_TDB, %s'% (row['Planet'], row['T1'].strftime("%b %d %Y, %H:%M:%S"), Inst))
        plt.savefig(r'%s/%.3f_%s_%s_%s_%s.jpg' % (PlotPath, EventPercent, Internal_Rank, row['Planet'], Inst, row['T1'].strftime("%b_%d_%Y_%H_%M_%S")))         
        plt.close() #This is still needed to release system resources after plot is saved. Must be explicit!

        #Air mass plot, only for viable events.
        #If you wish to plot interpolation, move this block into the f3 for loop.
#        plt.figure()
#        plt.scatter(delta_midnight, abs(np.sin(((Targetaltazs.alt).degree)*(np.pi/180)))**(-0.6), c = Targetaltazs.az, label = 'Target Path', lw=0, s=8,cmap='winter') #Target path
#        #plt.plot(delta_midnight, f3(delta_midnight), 'r--', label= 'Interpolation')
#        #Dark fill starts once sun reaches -18 deg alt, superimposed on twilight zone.
#        plt.fill_between(delta_midnight.to('hr').value, 0, 90,sunaltazs.alt < -18*u.deg, color='black', zorder=0, label='Astronomical Twilight') 
#        #Red fill denotes transit in progress, superimposed on other zones. alpha denotes transparency.
#        plt.axvspan(difh1.value,difh2.value, 0,90, alpha=0.3, color='r', zorder=0, label='Observation')
#        plt.axvspan(difh3.value,difh4.value, 0,90, alpha=0.3, color='r', hatch='X', zorder=0, label='Transit Active')
#        plt.colorbar().set_label('Azimuth [deg]')
#        plt.legend(loc='upper right')
#        plt.xlim(-12, 18)
#        plt.xticks(np.arange(15)*2 -12)
#        plt.ylim(0, 60)
#        plt.xlabel('Hours from BJD_TDB Midnight')
#        plt.ylabel('Air Mass')
#        plt.title(r'%s - Air Mass Plot for %s BJD_TDB' % (row['Planet'], row['T1'].strftime("%b %d %Y, %H:%M:%S")))
#        plt.close()
            
        return Phase, min(Sep), Internal_Rank, ObsStart_JD, ObsEnd_JD, Air_Mass_Marker, WBase, WTrans, WInOut, EventPercent
        #Note the conversion to UTC for our observing start/end times!
    
    def FinalMetric(row):
        
        FMet = row['%s' % Metric_Mode]*row['Air_Mass_Metric']*row['Event_Weight']
        
        return FMet

    #Altitude spot-check - will the target ever be high enough to be observed?
    di[['RA:HMS','Dec:Deg']]=di.apply(lambda row: pd.Series(AltChecker(row)), axis=1)
    
    if PT8[0] >= (90.-abs((np.float(ds['Lat'][S[0]:int(S[0]+1.)]))-(di['Dec:Deg'][0]))):
        p((r'Phase Two fired for %s - %s - not observable from chosen location.') % (di['Planet'][0], PlanetNo+1) ) # CHANGED TO REMOVE MP
        #p((r'Phase Two fired for %s - %s of %s targets for Core %s - not observable from chosen location.') % 
          #(di['Planet'][0], PlanetNo+1, len(JobChunk[CoreNo]), CoreNo+1))  
    else:    
        #Generates bank of T1s and errors, writes as two columns.
        #Then fills in existing in planet information with fillna - one row becomes many!
        di['DummyCol']=di.apply(lambda row: T1Calc(row),axis=1)
        di.drop('DummyCol',axis=1,inplace=True)
            
        if len(TargetTimes) == 0.:
            p((r'Phase Two fired for %s - %s - no transits in chosen window.') % (di['Planet'][0], PlanetNo+1) ) # CHANGED TO REMOVE MP
            #p((r'Phase Two fired for %s - %s of %s targets for Core %s - no transits in chosen window.') % 
              #(di['Planet'][0], PlanetNo+1, len(JobChunk[CoreNo]), CoreNo+1))
        else:    
            di=pd.concat([di,pd.DataFrame(TargetTimes)], axis=1) #Set axis=1 to combine column-wise. 
            dj=pd.concat([di,pd.DataFrame(ErrTimes)], axis=1)
            dj.columns=['Planet','R*','Rp','RA:HH','RA:MM','RA:SS','Dec:DD','Dec:MM','Dec:SS','T14','Depth','T0 (HJD or BJD)',
            'T0_err','P (day)','P_err', 'Impact Parameter','a_Calc','Previous Study Flag','Rank','Habitable_Rank',
            'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank', 'RA:HMS','Dec:Deg','T1','T1_err']
            dj.fillna(method='ffill',inplace=True)
            
            #Writes additional key times.
            dj[['T2','T0','T3','T4','Baseline_Start','Baseline_End']]=dj.apply(lambda row: pd.Series(TimesCalc(row)), axis=1)
            dj[['Lunar_Illumination','Closest_Lunar_Approach','Internal Rank','Observation_Start_(UTC)','Observation_End_(UTC)',
                'Air_Mass_Metric', 'Baseline_Weight','Transit_Curve_Weight','Ingress-Egress_Weight', 'Event_Weight']] = dj.apply(lambda row: pd.Series(EventMetric(row)), axis=1)
           
            #Block to catch colour-bar issues, if needed in future.
            # try:
            #     dj[['Lunar_Illumination','Closest_Lunar_Approach','Internal Rank','Observation_Start_(UTC)','Observation_End_(UTC)',
            #     'Air_Mass_Metric', 'Baseline_Weight','Transit_Curve_Weight','Ingress-Egress_Weight', 'Event_Weight']] = dj.apply(lambda row: pd.Series(EventMetric(row)), axis=1)
            # except u.core.UnitConversionError:
            #     sys.exit('Planet %s has done a bad thing - Core No. %s dropped.' % (dj['Planet'][0],CoreNo+1))    
                
            dj['Final_%s' % Metric_Mode]=dj.apply(lambda row: FinalMetric(row), axis=1)
            dk=dj.sort_values(by ='Final_%s' % Metric_Mode, ascending=False)
            dk.to_csv((r'%s/%s') % (csvoutpath, Name),index=False)
            
            VGen=(dj['Internal Rank'] == 'X').all()
            if VGen == True:
                p((r'Phase Two fired for %s - %s - no viable events found.') % (di['Planet'][0], PlanetNo+1) ) # CHANGED TO REMOVE MP
                #p((r'Phase Two fired for %s - %s of %s targets for Core %s - no viable events found.') % 
                  #(di['Planet'][0], PlanetNo+1, len(JobChunk[CoreNo]), CoreNo+1))
            else:
                p((r'Phase Two fired for %s - %s') % (dj['Planet'][0], PlanetNo+1) ) # CHANGED TO REMOVE MP
                #p((r'Phase Two fired for %s - %s of %s targets for Core %s.') % 
                  #(dj['Planet'][0], PlanetNo+1, len(JobChunk[CoreNo]), CoreNo+1))
    
    #CLEARANCE PROTOCOLS FOR MULTIPLE-TARGET RUNNING   
    del TargetTimes[:]
    del ErrTimes[:]
    del z_hold[:] 
    
    return get_print() #Return all printed statements as part of main line.