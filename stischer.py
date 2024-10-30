import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tkinter import filedialog
import os
import time
# This code assumes that your spectra are saved in text files with columns labeled 'wave', 'flux', and 'error'.
# The data files should be sorted by wavelength, and the columns separated by tabs.

class Spectrum:
    #A linked list of spectra
    def __init__(self, df_in):
        self.wave = df_in['wave'].to_numpy()
        self.flux = df_in['flux'].to_numpy()
        self.error = df_in['error'].to_numpy()
        self.next = None

    def concatenate(self):
        self.wave = np.concatenate([self.wave, self.next.wave])
        self.flux = np.concatenate([self.flux, self.next.flux])
        self.error = np.concatenate([self.error, self.next.error])
        self.next = self.next.next
        print("concatenation complete!")
        return self

    def overlap_mean(self):
        wave1 = self.wave
        flux1 = self.flux
        error1 = self.error

        wave2 = self.next.wave
        flux2 = self.next.flux
        error2 = self.next.error
        flux2_interp = np.interp(wave1,wave2,flux2)
        error2_interp = np.interp(wave1,wave2,error2)

        weight1 = np.divide(1,np.power(error1,2))
        weight2 = np.divide(1,np.power(error2_interp,2))

        mean_flux = np.divide(np.add(np.multiply(flux1,weight1),np.multiply(flux2_interp,weight2)),np.add(weight1,weight2))
        mean_error = np.sqrt(np.divide(1,np.add(weight1,weight2)))
        next_next = self.next.next

        mean_spectrum = Spectrum(pd.DataFrame({'wave':wave1,'flux':mean_flux,'error':mean_error}))
        mean_spectrum.next = next_next
        return mean_spectrum

    def stitch(self,is_overlap=False):
        fig, ax = plt.subplots()
        ax.plot(self.wave, self.flux)

        plt.show(block=False)
        if self.next is None:
            print("This linked list is done!")
            return self
        #print(self.next.wave)
        fig,ax = plt.subplots()
        ax.plot(self.wave, self.flux)
        ax.plot(self.next.wave, self.next.flux)
        plt.xlim(min(min(self.wave),min(self.next.wave)), max(max(self.wave),max(self.next.wave)))
        plt.ylim(0,max(self.flux))
        plt.show(block=False)
        #check if the first two spectra are ordered, or if the spectra were called in as overlaps
        if min(self.wave) > min(self.next.wave) and (is_overlap is False):
            curr_wave = self.wave
            curr_flux = self.flux
            curr_error = self.error
            next_wave = self.next.wave
            next_flux = self.next.flux
            next_error = self.next.error
            self.wave = next_wave
            self.flux = next_flux
            self.error = next_error
            self.next.wave = curr_wave
            self.next.flux = curr_flux
            self.next.error = curr_error
            return self.stitch()
        elif max(self.wave) < min(self.next.wave):
            concatenated =  self.concatenate()
            return concatenated.stitch()
        #OVERLAP TIME!
        else:
            #Check for gaps in the first spectrum

            for i in range(len(self.wave)-1):
                if (self.wave[i+1]-self.wave[i]) > 10: #set the maximum tolerable width of a gap to be 10 angstroms
                    if self.wave[i+1] > min(self.next.wave):
                        wave1 = self.wave[0:i+1]
                        flux1 = self.flux[0:i+1]
                        error1 = self.error[0:i+1]
                        wave2 = self.next.wave
                        flux2 = self.next.flux
                        error2 = self.next.error
                        wave3 = self.wave[i+1:]
                        flux3 = self.flux[i+1:]
                        error3 = self.error[i+1:]
                        next_next = self.next.next
                        spectrum2 = Spectrum(pd.DataFrame({'wave':wave2, 'flux':flux2, 'error':error2}))
                        spectrum3 = Spectrum(pd.DataFrame({'wave':wave3, 'flux':flux3, 'error':error3}))
                        spectrum3.next = next_next
                        spectrum2.next = spectrum3
                        self.next = spectrum2
                        self.wave= wave1
                        self.flux = flux1
                        self.error = error1
                        fig,ax = plt.subplots()
                        ax.plot(wave1, flux1)
                        ax.plot(wave2, flux2)
                        ax.plot(wave3, flux3)
                        plt.show(block=False)
                        return self.stitch(is_overlap=is_overlap)
            #Check for gaps in the second spectrum
            for j in range(len(self.next.wave)-1):
                if (self.next.wave[j+1]-self.next.wave[j]) > 10:
                    wave2 = self.next.wave[0:j+1]
                    flux2 = self.next.flux[0:j+1]
                    error2 = self.next.error[0:j+1]
                    wave3 = self.next.wave[j+1:]
                    flux3 = self.next.flux[j+1:]
                    error3 = self.next.error[j+1:]
                    next_next = self.next.next
                    spectrum2 = Spectrum(pd.DataFrame({'wave':wave2, 'flux':flux2, 'error':error2}))
                    spectrum3 = Spectrum(pd.DataFrame({'wave':wave3, 'flux':flux3, 'error':error3}))
                    spectrum3.next = next_next
                    spectrum2.next = spectrum3
                    self.next = spectrum2
                    fig, ax = plt.subplots()
                    ax.plot(self.wave, self.flux)
                    ax.plot(wave2, flux2)
                    ax.plot(wave3, flux3)
                    plt.show(block=False)
                    return self.stitch(is_overlap=is_overlap)
            #No gaps now
            if self.next is not None and self.wave[0] == self.next.wave[0]:
                is_overlap = True
            if is_overlap is False:
                index1 = max(np.nonzero(self.wave<min(self.next.wave))[0])
                wave1 = self.wave[0:index1+1]
                flux1 = self.flux[0:index1+1]
                error1 = self.error[0:index1+1]
                wave2 = self.wave[index1+1:]
                flux2 = self.flux[index1+1:]
                error2 = self.error[index1+1:]
                spectrum2 = Spectrum(pd.DataFrame({'wave':wave2, 'flux':flux2, 'error':error2}))
                next_next = self.next.next
                spectrum3 = self.next
                spectrum2.next = spectrum3
                fig, ax = plt.subplots()
                ax.plot(wave1, flux1)
                ax.plot(wave2, flux2)
                ax.plot(spectrum3.wave, spectrum3.flux)
                plt.show(block=False)
                stitched_spectrum2 = spectrum2.stitch(is_overlap=True)
                spectrum1 = Spectrum(pd.DataFrame({'wave':wave1, 'flux':flux1, 'error':error1}))
                spectrum1.next = stitched_spectrum2
                fig, ax = plt.subplots()
                ax.plot(spectrum1.wave, spectrum1.flux)
                ax.plot(stitched_spectrum2.wave, stitched_spectrum2.flux)
                stitched_spectrum1 = spectrum1.stitch(is_overlap=False)
                self.wave = stitched_spectrum1.wave
                self.flux = stitched_spectrum1.flux
                self.error = stitched_spectrum1.error
                self.next = next_next
                return self.stitch(is_overlap=False)

            if is_overlap:
                gap_checked = False
                if max(self.wave) > max(self.next.wave):
                    index3 = min(np.nonzero(self.wave > max(self.next.wave))[0])
                    wave1 = self.wave[0:index3-1]
                    flux1 = self.flux[0:index3-1]
                    error1 = self.error[0:index3-1]
                    wave2 = self.next.wave
                    flux2 = self.next.flux
                    error2 = self.next.error
                    wave3 = self.wave[index3:]
                    flux3 = self.flux[index3:]
                    error3 = self.error[index3:]
                    gap_checked = True
                elif max(self.next.wave) > max(self.wave):
                    index3 = min(np.nonzero(self.next.wave > max(self.wave))[0])
                    wave1 = self.wave
                    flux1 = self.flux
                    error1 = self.error
                    wave2 = self.next.wave[0:index3-1]
                    flux2 = self.next.flux[0:index3-1]
                    error2 = self.next.error[0:index3-1]
                    wave3 = self.next.wave[index3:]
                    flux3 = self.next.flux[index3:]
                    error3 = self.next.error[index3:]
                    gap_checked = True
                if gap_checked:
                    fig,ax = plt.subplots()
                    ax.plot(wave1, flux1)
                    ax.plot(wave2, flux2)
                    ax.plot(wave3, flux3)
                    next_next = self.next.next
                    print(wave3)
                    spectrum3 = Spectrum(pd.DataFrame({'wave':wave3, 'flux':flux3, 'error':error3}))
                    spectrum3.next = next_next
                    spectrum2 = Spectrum(pd.DataFrame({'wave':wave2, 'flux':flux2, 'error':error2}))
                    spectrum2.next = spectrum3
                    self.wave = wave1
                    self.flux = flux1
                    self.error = error1
                    self.next = spectrum2
                overlap_mean = self.overlap_mean()
                return overlap_mean.stitch()
        return self.stitch()


print("Welcome to STIScher. Please make sure that your spectrum files are collected in the same folder.")
print("Please enter the folder that contains the files.")
print("In the file selection dialog, make sure you're *in* the folder that contains the data files before clicking 'OK'.")
print("Simply selecting the folder may not work.")
directory_name = filedialog.askdirectory()
print(directory_name)
same_grating = input("Are the spectra all of the same grating? (y/n)")
if same_grating=='y' or same_grating=='Y':
    overlap = True
else:
    overlap = False

spectra = None
spectrum_current = None
for root, dirs, files in os.walk(directory_name):
    for file in files:
        print(os.path.join(root,file))
        if ".txt" in file:
            df = pd.read_csv(os.path.join(root, file), sep='\t')
            if spectra is None:
                spectra = Spectrum(df)
                spectrum_current = spectra
            else:
                spectrum_current.next = Spectrum(df)
                spectrum_current = spectrum_current.next

stitched_spectrum = spectra.stitch(overlap)
final_df = pd.DataFrame({'wave':stitched_spectrum.wave,'flux':stitched_spectrum.flux,'error':stitched_spectrum.error})
final_df.to_csv(directory_name+'/stitched_spectrum.txt', sep='\t',index=False)