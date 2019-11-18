import numpy as np 
import matplotlib.pyplot as plt 
from scipy.signal import savgol_filter
import sys


"""
Stuff that we can get from the outputfiles:
cycles = np.genfromtxt(inputfile, usecols=0, skip_header=1)
counter = np.genfromtxt(inputfile, usecols=1, skip_header=1)
temp = np.genfromtxt(inputfile, usecols=2, skip_header=1)
energy_avg = np.genfromtxt(inputfile, usecols=3, skip_header=1)
heat_cap = np.genfromtxt(inputfile, usecols=4, skip_header=1)    
magn_mom = np.genfromtxt(inputfile, usecols=5, skip_header=1)
magn_susp = np.genfromtxt(inputfile, usecols=6, skip_header=1)
magn_mom_abs = np.genfromtxt(inputfile, usecols=7, skip_header=1)
variance_E = np.genfromtxt(inputfile, usecols=8, skip_header=1)
variance_M = np.genfromtxt(inputfile, usecols=9, skip_header=1)
energy = np.genfromtxt(inputfile, usecols=10, skip_header=1)
"""



def plotting_4c(inputfile, inputfile_random, temp, stable_index):
    cycles = np.genfromtxt(inputfile, usecols=0, skip_header=1)[:stable_index]
    counter = np.genfromtxt(inputfile, usecols=1, skip_header=1)[:stable_index]
    energy_avg = np.genfromtxt(inputfile, usecols=3, skip_header=1)[:stable_index]
    magn_mom_abs = np.genfromtxt(inputfile, usecols=7, skip_header=1)[:stable_index]

    cycles_random = np.genfromtxt(inputfile_random, usecols=0, skip_header=1)[:stable_index]
    counter_random = np.genfromtxt(inputfile_random, usecols=1, skip_header=1)[:stable_index]
    energy_avg_random = np.genfromtxt(inputfile_random, usecols=3, skip_header=1)[:stable_index]
    magn_mom_abs_random = np.genfromtxt(inputfile_random, usecols=7, skip_header=1)[:stable_index]
   
    ax1 = plt.subplot(311)
    plt.plot(cycles, (energy_avg)/cycles, label=r'$E_{avg}$', color='olive')
    plt.plot(cycles_random, (energy_avg_random)/cycles_random, label=r'$E_{avg}$, random', color='purple')
    plt.legend(loc='best')
    plt.title(r'Most likely state reached for T = {}$[kT/J]$'.format(temp))
    plt.ylabel('Energy')
    
    ax2 = plt.subplot(312, sharex = ax1)
    plt.plot(cycles, (magn_mom_abs)/cycles, label=r'Magnetic moment', color='olive')
    plt.plot(cycles_random, (magn_mom_abs_random)/cycles_random, label=r'Magnetic moment, random', color='purple')
    plt.legend(loc='best')
    plt.ylabel('Magnetic moment')

    ax3 = plt.subplot(313, sharex = ax1)
    plt.plot(cycles, (counter)/cycles, label='Accepted configurations', color='olive')
    plt.plot(cycles_random, (counter_random)/cycles_random, label='Accepted configurations, random', color='purple')
    plt.legend(loc='best')
    plt.title(r'Total number of accepted configurations for T = {}$[kT/J]$'.format(temp))
    plt.ylabel('No. of configurations')
    plt.xlabel('Monte Carlo cycles')
    plt.show()

def plotting_4d(inputfile, temp, stable_index):
    energy = np.genfromtxt(inputfile, usecols=10, skip_header=1)
    variance_E = np.genfromtxt(inputfile, usecols=8, skip_header=1)
    print(len)
    print('Mean variance in energy: ', np.mean(variance_E))
    plt.figure()
    plt.hist(energy[stable_index:], bins=40)
    plt.title(r'Probability distribution of the energy T = {}$[kT/J]$'.format(temp))
    plt.ylabel('No. of events')
    plt.xlabel('Energy')
    plt.show()


def plotting_4e(L40, L60, L80, L100):
    temp_40 = np.genfromtxt(L40, usecols=2, skip_header=1)
    energy_avg_L40 = np.genfromtxt(L40, usecols=3, skip_header=1)
    heat_cap_L40 = np.genfromtxt(L40, usecols=4, skip_header=1)   
    magn_susp_L40 = np.genfromtxt(L40, usecols=6, skip_header=1) 
    magn_mom_abs_L40 = np.genfromtxt(L40, usecols=7, skip_header=1)
    
    energy_avg_L60 = np.genfromtxt(L60, usecols=3, skip_header=1)
    heat_cap_L60 = np.genfromtxt(L60, usecols=4, skip_header=1)    
    magn_susp_L60 = np.genfromtxt(L60, usecols=6, skip_header=1) 
    magn_mom_abs_L60 = np.genfromtxt(L60, usecols=7, skip_header=1)

    energy_avg_L80 = np.genfromtxt(L80, usecols=3, skip_header=1)
    heat_cap_L80 = np.genfromtxt(L80, usecols=4, skip_header=1)    
    magn_susp_L80 = np.genfromtxt(L80, usecols=6, skip_header=1) 
    magn_mom_abs_L80 = np.genfromtxt(L80, usecols=7, skip_header=1)

    energy_avg_L100 = np.genfromtxt(L100, usecols=3, skip_header=1)
    heat_cap_L100 = np.genfromtxt(L100, usecols=4, skip_header=1)    
    magn_susp_L100 = np.genfromtxt(L100, usecols=6, skip_header=1) 
    magn_mom_abs_L100 = np.genfromtxt(L100, usecols=7, skip_header=1)

    smooth_fac = 15
    order = 3
    heat_cap_L40_smoothed = savgol_filter(heat_cap_L40, window_length=smooth_fac, polyorder=order)
    heat_cap_L60_smoothed = savgol_filter(heat_cap_L60, window_length=smooth_fac, polyorder=order)
    heat_cap_L80_smoothed = savgol_filter(heat_cap_L80, window_length=smooth_fac, polyorder=order)
    heat_cap_L100_smoothed = savgol_filter(heat_cap_L100, window_length=smooth_fac, polyorder=order)

    magn_susp_L40_smoothed = savgol_filter(magn_susp_L40, window_length=smooth_fac, polyorder=order)
    magn_susp_L60_smoothed = savgol_filter(magn_susp_L60, window_length=smooth_fac, polyorder=order)
    magn_susp_L80_smoothed = savgol_filter(magn_susp_L80, window_length=smooth_fac, polyorder=order)
    magn_susp_L100_smoothed = savgol_filter(magn_susp_L100, window_length=smooth_fac, polyorder=order)


    indx_40 = np.argmax(heat_cap_L40_smoothed)
    Cv_40 = heat_cap_L40_smoothed[indx_40]
    Tc_40 = temp_40[indx_40]

    indx_60 = np.argmax(heat_cap_L60_smoothed)
    Cv_60 = heat_cap_L60_smoothed[indx_60]
    Tc_60 = temp_40[indx_60]

    indx_80 = np.argmax(heat_cap_L80_smoothed)
    Cv_80 = heat_cap_L80_smoothed[indx_80]
    Tc_80 = temp_40[indx_80]

    indx_100 = np.argmax(heat_cap_L100_smoothed)
    Cv_100 = heat_cap_L100_smoothed[indx_100]
    Tc_100 = temp_40[indx_100]    

    plt.subplot(221)
    plt.title(r'Expectaion values for the energy $\langle E \rangle$')
    plt.plot(temp_40, energy_avg_L40, label=r'$L = 40$', color='purple')
    plt.plot(temp_40, energy_avg_L60, label=r'$L = 60$', color='olive')
    plt.plot(temp_40, energy_avg_L80, label=r'$L = 80$', color='coral')
    plt.plot(temp_40, energy_avg_L100, label=r'$L = 100$', color='teal')
    plt.legend(loc='best')
    plt.ylabel(r'$\langle E \rangle$')
    plt.xlabel(r'$T[K]$')

    plt.subplot(222)
    plt.title(r'Expectaion values for the magnetisation $\langle |M| \rangle$')
    plt.plot(temp_40, magn_mom_abs_L40, label=r'$L = 40$', color='purple')
    plt.plot(temp_40, magn_mom_abs_L60, label=r'$L = 60$', color='olive')
    plt.plot(temp_40, magn_mom_abs_L80, label=r'$L = 80$', color='coral')
    plt.plot(temp_40, magn_mom_abs_L100, label=r'$L = 100$', color='teal')
    plt.legend(loc='best')
    plt.ylabel(r'$\langle |M| \rangle$')
    plt.xlabel(r'$T[K]$')

    plt.subplot(223)
    plt.title(r'Spesific heat $C_V$')
    plt.plot(temp_40, heat_cap_L40_smoothed, label=r'$L = 40$', color='purple')
    plt.plot(temp_40, heat_cap_L60_smoothed, label=r'$L = 60$', color='olive')
    plt.plot(temp_40, heat_cap_L80_smoothed, label=r'$L = 80$', color='coral')
    plt.plot(temp_40, heat_cap_L100_smoothed, label=r'$L = 100$', color='teal')
    plt.plot(Tc_40, Cv_40,'o', label=r'$T_C \approx ${} for $L = 40$'.format(Tc_40), color='purple')
    plt.plot(Tc_60, Cv_60,'o', label=r'$T_C \approx ${} for $L = 60$'.format(Tc_60), color='olive')
    plt.plot(Tc_80, Cv_80,'o', label=r'$T_C \approx ${} for $L = 80$'.format(Tc_80), color='coral')
    plt.plot(Tc_100, Cv_100,'o', label=r'$T_C \approx ${} for $L = 100$'.format(Tc_100), color='teal')
    plt.legend(loc='best')
    plt.ylabel(r'$C_V$')
    plt.xlabel(r'$T[K]$')

    plt.subplot(224)
    plt.title(r'Magnetic susceptibility $\chi$')
    plt.plot(temp_40, magn_susp_L40_smoothed, label=r'$L = 40$', color='purple')
    plt.plot(temp_40, magn_susp_L60_smoothed, label=r'$L = 60$', color='olive')
    plt.plot(temp_40, magn_susp_L80_smoothed, label=r'$L = 80$', color='coral')
    plt.plot(temp_40, magn_susp_L100_smoothed, label=r'$L = 100$', color='teal')
    plt.legend(loc='best')
    plt.ylabel(r'$\chi$')
    plt.xlabel(r'$T[K]$')
    plt.show()

    L_lst = [40, 60, 80, 100]
    Tc_lst = [Tc_40, Tc_60, Tc_80, Tc_100]

    L_array = np.flip(L_lst)
    Tc_array = np.flip(Tc_lst)
    plt.plot(L_array, Tc_array, 'o')
    plt.show()

    Tc_infty = np.polyfit(L_lst, Tc_lst, deg=1)[1]
    print('Critical temperature, numerical ',Tc_infty)


T10 = 'output_T10.txt'
T10_random = 'output_T10_random.txt'

T24 = 'output_T24.txt'
T24_random = 'output_T24_random.txt'


L40 = 'output_L40_10e5.txt'
L60 = 'output_L60_10e5.txt'
L80 = 'output_L80_10e5.txt'
L100 = 'output_L100_10e5.txt'

L40_random = 'output_L40_10e5_random.txt'
L60_random = 'output_L60_10e5_random.txt'
L80_random = 'output_L80_10e5_random.txt'
L100_random = 'output_L100_10e5_random.txt'


plotting_4c(T10, T10_random, 1.0, 3000)
plotting_4c(T24, T24_random, 2.4, 7500)
#plotting_4d(T10, 1.0, 7500)
#plotting_4d(T24, 2.4, 7500)

#plotting_4e(L40, L60, L80, L100)
#plotting_4e(L40_random, L60_random, L80_random, L100_random)
