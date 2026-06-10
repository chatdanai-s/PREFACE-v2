# Part 3 of overall pipeline.
# Takes in assembled TEPCat catalogue and calculates all required telescope/filter/location/time-independent parameters.
# Updates data set in conjunction with ModCheck, outputs as WorkingTEPSet.

import pandas as pd
import numpy as np
import astropy.units as u

from bisect import bisect_left

# Functions to calculate area ratio, row-by-row.
def AreaRatio(row):
    Rp = (row['Rp'] * u.jupiterRad).to(u.solRad)
    Rs = (row['R*'] * u.solRad)
    AreaRatio = (Rp/Rs) ** 2
    return AreaRatio
    
# Proxy for scale height, units of km.
# Assumes atmosphere of solar metallicity ideal gas.
# One entry has gravity of 0; exception is needed to avoid dividing by 0!
def H_Proxy(row):
    kB = 1.380649e-23
    mu_m = 2.33 * 1.6605e-27
    gp = row['gp']
    Teq = row['Teq']

    if gp == 0.:
        return np.nan
    
    if Teq == np.nan:
        Teq = row['Teq_Calc']
        
    scaleHeight = ((kB * Teq) / (mu_m * gp)) * u.m
    scaleHeight = (scaleHeight.to(u.km)).value
    return scaleHeight

# Calculates equilibrium temperature for all planets using stellar temps.
def Teq_Calc(row):
    # Recover semi-major axis from period and stellar mass. (Assume Mp << M*)
    P = (row['P (day)'] * u.d).to(u.yr)
    Ms = row['M*']
    a = (np.cbrt(P**2 * Ms)).value
    
    # For Bond albedo (across all EM wavelengths) = 0 and full heat redistribution.
    Teff = row['Teff']
    Rs = (row['R*'] * u.solRad).to(u.AU)
    Teq = Teff * np.sqrt(Rs / (2 * a*u.AU))

    return a, Teq

# Returns a flag based on the three supplied references for each object.
# If all different, well-studied. If all identical, not followed up at all!
def RefNo_Flag(row):
    all_identical = row['Discovery_Reference'] == row['Recent_Reference'] == row['Ephemeris_Reference']
    
    all_different = all([row['Discovery_Reference'] != row['Recent_Reference'],
                         row['Recent_Reference'] != row['Ephemeris_Reference'],
                         row['Discovery_Reference'] != row['Ephemeris_Reference']
                        ])
        
    if all_different == True:
        RefNo = 3.  # Well-studied!
    elif all_identical == True:
        RefNo = 1.  # Poorly studied :(
    else:
        RefNo = 2.  # Indecisive... 笑いましょう (TL Note: Let's laugh)
        
    return RefNo

# Bisection of temp list - repeatedly halves list and works out which half the temperature has to be in.
# Converges on answer over repeated bisections - very efficient! (Time on order log(n))
def takeClosest(Tsort, i):
    pos = bisect_left(Tsort, i)
    if pos == 0:
        return Tsort[0]
    
    if pos == len(Tsort):
        return Tsort[-1]
    
    before = Tsort[pos-1]
    after = Tsort[pos]

    if after - i < i - before:
        return after
    else:
        return before
     
# Central function to build WorkingTEPSet.
def WorkBuilder(CSV_core_folder):
    # Locate files (formerly imported as di, dj, dc)
    tepset_csv_path = CSV_core_folder / 'FullTEPSetWithAllImpacts.csv'
    speccolor_csv_path = CSV_core_folder / 'SpecConverter.csv'
    vmagcorr_csv_path = CSV_core_folder / 'V_Spectral_Corrections_Bank.csv'
    workingtep_csv_path = CSV_core_folder / 'WorkingTEPSet.csv'


    # TEPCat only keeps V and K-magnitudes - any others must be calculated for each star.
    # See www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt for conversion table.          
    # Reads in bank of spectral conversions, creates sorted list of temperature cut-offs.
    # Also reads in a table of corrections such that V_TEPCat = V_calculated.

    df_tep = pd.read_csv(tepset_csv_path, skipinitialspace=True)

    Fields = ['Teff_SpT', 'B-V', 'U-B', 'V-Rc', 'V-Ic', 'V-Ks', 'J-H', 'H-K']
    df_spec = pd.read_csv(speccolor_csv_path, usecols=Fields, skipinitialspace=True, na_values=0)
    df_spec['U-B'] = df_spec['U-B'].replace('...', np.nan)

    df_corr = pd.read_csv(vmagcorr_csv_path)


    Tmps = df_spec['Teff_SpT'].tolist()     # Descending temperature of spectral type
    Cors = df_corr['Vmag_List'].tolist()    # List of Vmags with correction
    Tsort = Tmps[::-1]                      # Make Tmps ascendingly sorted
    

    def MagCalc(row):            
        # Returns appropriate conversion factors for each entry and calculates new magnitudes.
        # Pandas does not like indexing with floats, so an integer must be used.
        # It must then go back to being a float for calculation purposes!
        # Correction C_f only applied for stars fainter than Vmag 13.
        # Calculated K-values from V become too faint - minus sign gives a brightening correction.
        # Rather than go to a third step, J-mag recovery builds on the H-mag transformation.
        # See SpecConverterExplorer.py for more.
        
        # Retrieve relevant rows of two dataframes
        R = np.where(df_spec['Teff_SpT'] == takeClosest(Tsort, row['Teff']))[0][0]  # Row of df_spec with closest Teff
        C = np.where(df_corr['Vmag_List'] == takeClosest(Cors, row['Vmag']))[0][0]  # Correction index of df_corr
        
        # Gather conversion factor
        B_f = df_spec['B-V'].iloc[R]
        R_f = df_spec['V-Rc'].iloc[R]
        I_f = df_spec['V-Ic'].iloc[R]
        K_f = df_spec['V-Ks'].iloc[R]
        
        U_f = df_spec['U-B'].iloc[R]
        K_f = df_spec['V-Ks'].iloc[R]
        H_f = df_spec['H-K'].iloc[R]
        J_f = df_spec['J-H'].iloc[R]

        C_f = df_corr['Correction'].iloc[C]

        # Calculate magnitudes
        Vmag_uncalc = row['Vmag']

        Bmag = Vmag_uncalc + B_f - C_f
        Rmag = Vmag_uncalc - R_f - C_f
        Imag = Vmag_uncalc - I_f - C_f
        Kmag = Vmag_uncalc - K_f - C_f

        # If entries have missing K-mags (eg. two OGLE objects), use the calculated ones instead!
        Kmag_uncalc = row['Kmag']
        if row['Kmag'] == '':
            Kmag_uncalc = Kmag

        # Calculate more magnitudes
        Umag = Bmag + U_f        # Correction already applied in V->B transformation.
        Vmag = Kmag_uncalc + K_f + C_f
        Hmag = Kmag_uncalc + H_f + C_f
        Jmag = Hmag + J_f

        return Bmag, Rmag, Imag, Kmag, Umag, Vmag, Hmag, Jmag
    

    # Applies functions to filtered data set and creates new columns.
    # For functions returning multiple variables, count over them as a Series to write to individual columns.
    magCols = ['Bmag','Rmag','Imag','Kmag_Calc','Umag','Vmag_Calc','Hmag','Jmag']

    df_tep['Depth_Calc'] = df_tep.apply(lambda row: AreaRatio(row), axis=1)
    df_tep[magCols] = df_tep.apply(lambda row: MagCalc(row),
                                   axis=1, result_type='expand')
    df_tep[['a_Calc','Teq_Calc']] = df_tep.apply(lambda row: Teq_Calc(row),
                                                 axis=1, result_type='expand') 
    df_tep['ScaleHeight_Prox'] = df_tep.apply(lambda row: H_Proxy(row), axis=1)    
    df_tep['Previous Study Flag'] = df_tep.apply(lambda row: RefNo_Flag(row), axis=1)

    # If entries have missing K-mags (eg. two OGLE objects), use the calculated ones instead!
    Kmag_missing = (df_tep['Kmag'] == '')
    df_tep.loc[Kmag_missing, 'Kmag'] = df_tep.loc[Kmag_missing, 'Kmag_Calc']


    # Only recovered for systems with no b-value yet in literature.
    # Seager/ORM model (2003) is simple - if it breaks, assign b of 0.5.
    # Can subsequently recover inclination in degrees.
    def ImpactCalc(row):
        a_diff = ((row['a (AU)'] - row['a_Calc']) / row['a_Calc']) * 100.
        SizeRatio = (row['Rp'] * u.jupiterRad).to(u.solRad) / (row['R*'] * u.solRad)

        # Use known a-values if available
        if pd.notnull(row['a (AU)']):
            a = row['a (AU)']
        else:
            a = row['a_Calc']
        
        SepRatio = (a * u.AU).to(u.solRad) / (row['R*'] * u.solRad)

        # To stop -ve square roots from ruining your day
        try:
            with np.errstate(invalid='raise'):
                piT14_P = np.pi * row['T14'] / row['P (day)']
                b = ((1+SizeRatio) / np.cos(piT14_P)) ** 2 - (SepRatio * np.tan(piT14_P)) ** 2
                b = np.sqrt(b)
        
        except FloatingPointError:
            b = np.nan
        
        return a_diff, SizeRatio, SepRatio, b
    
    # Check for validity of spectral conversions.        
    def SpecCheck(row):
        delta_Vmag = row['Vmag'] - row['Vmag_Calc']
        delta_Kmag = row['Kmag'] - row['Kmag_Calc']
        return delta_Vmag, delta_Kmag
    

    # More columns
    df_tep[['a-a_Calc % Error','Rp/R*','a_Calc/R*', 'Impact Parameter (ORM)']] = df_tep.apply(lambda row: ImpactCalc(row),
                                                                                              axis=1, result_type='expand')
    df_tep[['V-Vcalc', 'K-Kcalc']] = df_tep.apply(lambda row: SpecCheck(row),
                                                  axis=1, result_type='expand')
    
    # Merge impact parameter columns and tidy up.   
    df_tep['Impact Parameter (ORM)'] = df_tep['Impact Parameter (ORM)'].fillna(value=0.5)
    df_tep['Impact Parameter'] = df_tep['Impact Parameter'].fillna(df_tep['Impact Parameter (ORM)'])
    df_tep = df_tep.drop(['Impact Parameter (ORM)'], axis=1)
    

    # Recover ugriz Sloan magnitudes for TNT, LT IO:O and others.
    # Equations used are from Table 3 of Jordi et al. (2006)
    def SloanCalc(row):
        if row['Vmag'] - row['Rmag'] <= 0.93:
            rMinusR = 0.267 * (row['Vmag']-row['Rmag']) + 0.088
        else:
            rMinusR = 0.77 * (row['Vmag']-row['Rmag']) - 0.37
        rmag = row['Rmag'] + rMinusR

        iMinusI = 0.247 * (row['Rmag']-row['Imag']) + 0.329
        imag = row['Imag'] + iMinusI
                  
        gminusr = 1.646 * (row['Vmag']-row['Rmag']) - 0.139
        gmag = rmag + gminusr
        
        zmag = ((row['Imag'] - imag + 0.397) / 0.386) + imag
        umag = ((row['Bmag'] - gmag - 0.15) / 0.175) + gmag
                      
        return umag, gmag, rmag, imag, zmag
    
    df_tep[['umag','gmag','rmag','imag','zmag']] = df_tep.apply(lambda row: SloanCalc(row),
                                                                axis=1, result_type='expand')
                
    # RISE uses a 720nm longpass filter - I and Z have strong overlap, so can't just be averaged.
    # RISE does not take routine photometric standards, so take these numbers with a pinch of salt!
    
    # SPRAT zero-points are given in approx. V-band (blue arm) and R-band (red arm), so we need these for the calculations.
    def RISECalc(row):
        F_I = 10. ** (-0.4 * (row['Imag']-24.814))
        F_z = 10. ** (-0.4 * (row['zmag']-24.840)) #z and Z are roughly equivalent (Bessell, 2005)
    
        F_S = F_z + F_I        
        m_RISE = row['Imag'] + (2.5 * np.log10(F_I/F_S))
        
        m_SPRAT_Red = row['Rmag']
        m_SPRAT_Blue = row['Vmag']
        
        return m_RISE, m_SPRAT_Red, m_SPRAT_Blue
    
    # Recover magnitudes in 600RI+19 grism for VLT FORS2 spectroscopic observations and JWST.
    # From ESO's exposure time calculator, mzp_(R+I) per dispersed element is approximately 19.2.
    def VLTGrismCalc(row):
        F_R = 10. ** (-0.4 * (row['Rmag']-28.13))
        F_I = 10. ** (-0.4 * (row['Imag']-27.23))
        
        F_S = F_R + F_I
        
        m_600RI = row['Rmag'] + (2.5 * np.log10(F_R/F_S))
        # m_600RI = (row['Rmag'] + row['Imag']) / 2.
        
        m_HK = (row['Hmag'] + row['Kmag']) / 2.
        m_G750M = row['Vmag']
        
        return m_600RI, m_HK, m_G750M
    
    #Mag recovery for NTT EFOSC2, using broad GR01 grism.
    def NTTGrismCalc(row):    
        F_U = 10. ** (-0.4 * (row['Umag']-24.39))
        F_B = 10. ** (-0.4 * (row['Bmag']-26.16))
        F_V = 10. ** (-0.4 * (row['Vmag']-26.28))
        F_R = 10. ** (-0.4 * (row['Rmag']-26.36))
        F_I = 10. ** (-0.4 * (row['Imag']-25.48))
        
        if np.isnan(F_U):
            F_S = F_B + F_V + F_R + F_I
        else:
            F_S = F_U + F_B + F_V + F_R + F_I
        
        m_GR01 = row['Vmag'] + (2.5 * np.log10(F_V/F_S))
        
        return m_GR01
    
    # Finalize working TEP set
    df_tep[['RISEmag', 'SPRAT_Redmag', 'SPRAT_Bluemag']] = df_tep.apply(lambda row: RISECalc(row),
                                                                        axis=1, result_type='expand')
    df_tep[['600RI+19mag', 'HKmag', 'G750Mmag']] = df_tep.apply(lambda row: VLTGrismCalc(row),
                                                                axis=1, result_type='expand')
    df_tep['GR01mag'] = df_tep.apply(lambda row: NTTGrismCalc(row), axis=1)
    
    # Write dataframe to csvs.
    df_tep.to_csv(workingtep_csv_path, index=False)

    print('[WorkingTEPSetBuilder] New Working TEPSet constructed.')
