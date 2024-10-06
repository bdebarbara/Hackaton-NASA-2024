# Import libraries
import numpy as np
import pandas as pd
from obspy import read
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import scipy

from scipy import signal
from matplotlib import cm    

def get_seismic_events(cat_directory, test_filename):
    data_directory = cat_directory
    mseed_file = f'{data_directory}{test_filename}.mseed'
    st = read(mseed_file)

    tr = st.traces[0].copy()
    tr_times = tr.times()
    tr_data = tr.data

    # Set the minimum frequency
    minfreq = 0.5
    maxfreq = 1

    # Going to create a separate trace for the filter data
    st_filt = st.copy()
    # st_filt.filter('bandpass',freqmin=minfreq,freqmax=maxfreq)
    tr_filt = st_filt.traces[0].copy()
    tr_times_filt = tr_filt.times()
    tr_data_filt = tr_filt.data

    f, t, sxx = signal.spectrogram(tr_data_filt, tr_filt.stats.sampling_rate)

    newsxx = np.transpose(sxx)

    smth_rng = 5
    for i in range(len(newsxx)):
        for j in range(len(newsxx[i])):
            for k in range(smth_rng):
                if i+k < len(newsxx):
                    newsxx[i][j] += newsxx[i+k][j]
            newsxx[i][j] = newsxx[i][j]/smth_rng

    ampl = []
    for i in range(len(newsxx)):
        ampl.append(np.sum(newsxx[i]))

    # smothing
    ampl = signal.savgol_filter(ampl, 5, 2)
    for i in range(len(ampl)):
        amp = 0
        for j in range(5):
            if i+j < len(ampl):
                amp += ampl[i+j]
        ampl[i] = amp/5

    avr_window = 1000

    # Find peaks
    peaks, _ = scipy.signal.find_peaks(ampl, distance=33)  # You can adjust parameters as needed

    def calc_local_avreage(avr_window, pos: int, ampl):
        print(type(pos))
        print(type(avr_window))
        _avr = 0
        avr_window_adj = avr_window
        for j in range(int(pos - int(avr_window/2)), int(pos + int(avr_window/2))):
            if j < 0 or j >= len(ampl):
                avr_window_adj -= 1
                continue
            _avr += ampl[j]
        _avr /= avr_window_adj
        return _avr

    _peaks = []
    for i in range(len(peaks)):
        # claculate the local average of the peak
        avr = 0
        avr = calc_local_avreage(avr_window, peaks[i], ampl)

        if ampl[peaks[i]] > avr:
            _peaks.append(peaks[i])
    peaks = _peaks

    avr = 0

    for i in peaks:
        avr += ampl[i]
    avr = avr / len(peaks)

    _peaks = []
    for i in range(len(peaks)):
        if ampl[peaks[i]] > avr:
            _peaks.append(peaks[i])
    peaks = _peaks

    _peaks = []
    newPeaks=[]

    for peak in peaks:
        if ampl[peak] > ampl[peak-10] and ampl[peak] > ampl[peak+10]:
            _peaks.append(peak)

    print(_peaks)

    for peak in _peaks:

        #print(peak, _peaks.index(peak), _peaks)
        print(ampl[peak])
        minAvg = peak
        maxAvg = peak
        avr = calc_local_avreage(avr_window, peak, ampl)
        while ampl[minAvg] > avr and minAvg > 0:
            minAvg -= 1

        while ampl[maxAvg] > avr and maxAvg < 2554:
            maxAvg += 1

        if maxAvg - minAvg > 30:
            print(maxAvg - minAvg)
            newPeaks.append(peak)

    
    return tr_times, tr_data, t, newPeaks, avr, sxx, tr_times_filt, f, ampl, newsxx


_cat_directory = './data/lunar/test/data/S12_GradeB/'

import os
import glob
import pandas as pd  # Assuming you want to work with pandas for CSV files

def process_csv_file(file_path):
    # This is where you can process the CSV file
    # For demonstration, let's just read and print the first few rows
    df = pd.read_csv(file_path)
    print(f"Processing {file_path}")
    print(df.head())  # Print the first few rows of the dataframe

def main(folder_path):
    # Use glob to find all .csv files in the folder
    csv_files = glob.glob(os.path.join(folder_path, '*.csv'))

    for n, csv_file in enumerate(csv_files):
        _test_filename = csv_file[:-4].split('\\')[-1]
        tr_times, tr_data, t, newPeaks, avr, sxx, tr_times_filt, f, ampl, newsxx = get_seismic_events(_cat_directory, _test_filename)

        fig = plt.figure(figsize=(10, 10))
        ax = plt.subplot(2, 1, 1)

        # Plot trace
        ax.plot(tr_times,tr_data)

        # Mark detection
        for pos in t[newPeaks]:
            plt.axvline(x=pos, color='red', linestyle='--', linewidth=2)

        # Make the plot pretty
        ax.set_xlim([min(tr_times),max(tr_times)])
        ax.set_ylabel('Velocity (m/s)')
        ax.axhline(y=avr, color='black', linestyle='--', linewidth=1)
        ax.set_xlabel('Time (s)')
        ax.set_title(f'{_test_filename}', fontweight='bold')

        ax2 = plt.subplot(2, 1, 2)
        vals = ax2.pcolormesh(t, f, sxx, cmap=cm.jet, vmax=1e-18)
        ax2.set_xlim([min(tr_times_filt),max(tr_times_filt)])
        ax2.set_ylim([0, max(f)])
        ax2.set_xlabel(f'Time (Day Hour:Minute)', fontweight='bold')
        ax2.set_ylabel('Frequency (Hz)', fontweight='bold')

        os.makedirs("images", exist_ok=True)
        output_path = os.path.join("images", f'{n}.png')
        plt.savefig(output_path)
        plt.close()



if __name__ == "__main__":
    folder_path = input("Enter the folder path containing CSV files: ")
    main(_cat_directory)


