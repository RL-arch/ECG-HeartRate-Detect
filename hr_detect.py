#coding： Ronald B Liu   15.11.2019
## IDE PyCharm Ubuntu 18.04
# ECG filter & Heartrate detection (single channel)


import numpy as np
import matplotlib.pylab as plt


class FIR_filter:
    def __init__(self, fir_input):
        self.offset = 0
        self.p = 0
        self.coeff = 0
        self.buffer = np.zeros(number_of_taps)
        self.input = fir_input

    def dofilter(self, v):
        #lms update, get tap input power, buffer 
        output = 0
        self.buf_val = self.p + self.offset
        self.buffer[self.buf_val] = v
        while self.buf_val >= self.p:
            output += (self.buffer[self.buf_val] * self.input[self.coeff])
            self.buf_val = self.buf_val - 1
            self.coeff = self.coeff + 1

        self.buf_val = self.p + number_of_taps - 1

        while self.coeff < number_of_taps:
            output += (self.buffer[self.buf_val] * self.input[self.coeff])
            self.buf_val = self.buf_val - 1
            self.coeff = self.coeff + 1

        self.offset = self.offset + 1

        if self.offset >= number_of_taps:
            self.offset = 0

        self.coeff = 0
        return output

def filterShift(f0,f1,f2):
    #filter inlcudes band stop (to remove 50hz) and high pass (to remove DC)
    fir_f_resp = np.ones(number_of_taps)
    #---------------band stop-------------------------------
    fir_f_resp[int((f1 / fs) * number_of_taps):int((f2 / fs) * number_of_taps) + 1] = 0
    fir_f_resp[number_of_taps - int((f2 / fs) * number_of_taps):number_of_taps - int((f1 / fs) * number_of_taps) + 1] = 0 #do mirror
    #---------------low stop (high pass)-------------------------
    fir_f_resp[0:int((f0 / fs) * number_of_taps) + 1] = 0
    fir_f_resp[number_of_taps - int((f0 / fs) * number_of_taps):number_of_taps] = 0 #do mirror

    fir_t = np.fft.ifft(fir_f_resp) #frequency response to time
    h_real = np.real(fir_t) #time real ~ casual
    h_shifted = np.zeros(number_of_taps)
    h_shifted[0:int(number_of_taps / 2)] = h_real[int(number_of_taps / 2):number_of_taps]
    h_shifted[int(number_of_taps / 2):number_of_taps] = h_real[0:int(number_of_taps / 2)]
    return h_shifted

class matched_filter(FIR_filter):
    def __init__(self, ecg_input):
        self.input = ecg_input
    def Rpeak_detection(self, template, oringin):
        fir_coeff = template[::-1]
        detected_array = np.zeros(len(oringin))
        fir_template = FIR_filter(fir_coeff)
        for i in range(len(oringin)):
            detected_array[i] = fir_template.dofilter(self.input[i])
        detected_output = detected_array * detected_array  # The signal is squared to improve the output
        return detected_output

class generateTemplate:
    def __init__(self):
        self

    def mexicanhat(self):
        t = np.linspace(-250, 250, 500) #in the range of taps
        temp = (2 / np.sqrt(3 * 35) * np.pi ** (1 / 4)) * \
                     (1 - (t ** 2 / 35 ** 2)) * np.exp((-t ** 2) / (2 * 35 ** 2))
        return temp

    def gaussian1OD(self):
        t = np.linspace(-250, 250, 500)
        temp = -t * np.exp((-t ** 2) / 50) / (125 * np.sqrt(2 * np.pi))
        return temp

    def gaussian(self):
        t = np.linspace(-250, 250, 500)
        temp = np.exp((-t ** 2) / 50) / (5 * np.sqrt(2 * np.pi))
        return temp

    def shannon(self):
        t = np.linspace(-250, 250, 500)
        temp = np.sqrt(100) * np.sinc(100 * t) * np.exp(2 * 1j * t * np.pi * 4)
        return temp

class detectMomentaryHeartRate:
    def __init__(self, inputlist):
        self.inputlist = inputlist

    def detectMomentaryHeartRate(self,Fs):
        list = self.inputlist  # Output from Matched filter
        BPM = []  # It will be the array of Peaks
        counter = 0  # auxiliary counter
        threshold = max(self.inputlist) * 0.05
        for i in range(len(list)):
            if list[i] > threshold:
                differenceTime = (i - counter)  # difference of time T in second, f = 1/T
                counter = i
                bpm = 1 / differenceTime * (60 * Fs) #1min
                if 200 > bpm > 40:  # Limits for the detection of momentary heart rate
                    BPM.append(bpm)  # Add this peak to the BPM array
        BPM = np.delete(BPM, 0) #remove
        return BPM

def do_detection(input, Fs, oringin):
    template = generateTemplate()  # Create the class TemplateMaker

    # generate the different templates
    gaussian = template.gaussian()
    devgaussian = template.gaussian1OD()
    shannon = template.shannon()
    mexicanHat = template.mexicanhat()

    # Matching Filtering the R peak of the signal
    calculHeartBeat = matched_filter(input)
    detgaussian = calculHeartBeat.Rpeak_detection(gaussian, oringin)
    det1ODgaussian = calculHeartBeat.Rpeak_detection(devgaussian, oringin)
    detshannon = calculHeartBeat.Rpeak_detection(shannon, oringin)
    detmexicanHat = calculHeartBeat.Rpeak_detection(mexicanHat, oringin)

    #  calculate its coefficients of the matched filter analytically from a mathematical formula
    momentary_heart_rate_gaussian = detectMomentaryHeartRate(detgaussian)
    momentary_heart_rate_gaussian10d = detectMomentaryHeartRate(det1ODgaussian)
    momentary_heart_rate_shannon = detectMomentaryHeartRate(detshannon)
    momentary_heart_rate_mexicanhat = detectMomentaryHeartRate(detmexicanHat)

    # get the result value
    MHRGaussian = momentary_heart_rate_gaussian.detectMomentaryHeartRate(Fs)
    MHRGaussian1OD = momentary_heart_rate_gaussian10d.detectMomentaryHeartRate(Fs)
    MHRShannon = momentary_heart_rate_shannon.detectMomentaryHeartRate(Fs)
    MHRMexicanHat = momentary_heart_rate_mexicanhat.detectMomentaryHeartRate(Fs)
    return MHRShannon

if __name__ == '__main__':
    print("Warten Sie mal, bitte die Warnung ignorieren....")
    #----------------------------initial---------------------------------------------------
    # Define the frequencies coefficients for the FIR filter
    f0 = 1  #remove DC
    f1 = 45 #Fstop
    f2 = 55 #Fpass
    # Initialise the script
    fs = 1000  # 1000 Hz sampling rate
    number_of_taps = 500 #taps

    #--------------------------------exclute---------------------------------------------------
    # Chosen the interval 10,000-40,000 that are 20 seconds due to fs = 1000.
    data = np.loadtxt("ecg2.dat")
    #channel 0 ~ time; channel 1 ~ left hand(?); channel 2 ~ right hand; channel 3 ~ ankle
    time = data[10000:40000, 0]    #[from 10000 to 40000, time chose channel 0]
    unfilteredSignal = data[10000:40000, 2]   #[from 10000 to 40000, ecg chose channel 3]

    Shift = filterShift(f0, f1, f2)
    a = FIR_filter(Shift)
    filteredSignal = np.zeros(len(unfilteredSignal))
    #do filter
    for i in range(len(unfilteredSignal)):
        filteredSignal[i] = a.dofilter(unfilteredSignal[i])

    moment_heartrate = do_detection(filteredSignal, fs, unfilteredSignal)
    #---------------plot-------------------------------------------------------------------

    plt.figure(1)
    plt.plot(time, unfilteredSignal)
    plt.title("Pre Filtered")
    plt.xlabel("Time (s)")
    plt.ylabel('Amplitude')

    plt.figure(2)
    plt.plot(filteredSignal)
    plt.title("Filtered")
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")

    plt.figure(3)
    plt.plot(moment_heartrate)
    plt.title("Shannon Derivative")
    plt.xlabel("Time (s)")
    plt.ylabel("Rate r(t)")

    print("Fertig!")
    plt.show()

