#coding： Ronald B Liu   15.11.2019
## IDE PyCharm Ubuntu 18.04
# ECG filter (single channel)

import numpy as np
import matplotlib.pylab as plt


class FIR_filter:
    def __init__(self, FIR_input):
        self.offset = 0
        self.p = 0
        self.coeff = 0
        self.buffer = np.zeros(number_of_taps)
        self.input = FIR_input

    def dofilter(self, v):
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
    FIR_f_resp = np.ones(number_of_taps)
    FIR_f_resp[int((f1 / fs) * number_of_taps):int((f2 / fs) * number_of_taps) + 1] = 0
    FIR_f_resp[number_of_taps - int((f2 / fs) * number_of_taps):number_of_taps - int((f1 / fs) * number_of_taps) + 1] = 0
    FIR_f_resp[0:int((f0 / fs) * number_of_taps) + 1] = 0
    FIR_f_resp[number_of_taps - int((f0 / fs) * number_of_taps):number_of_taps] = 0

    FIR_t = np.fft.ifft(FIR_f_resp)
    h_real = np.real(FIR_t)
    h_shifted = np.zeros(number_of_taps)
    h_shifted[0:int(number_of_taps / 2)] = h_real[int(number_of_taps / 2):number_of_taps]
    h_shifted[int(number_of_taps / 2):number_of_taps] = h_real[0:int(number_of_taps / 2)]
    return h_shifted

if __name__ == '__main__':
    print("Warten Sie mal...")
    # Define the frequencies for the FIR filter
    f0 = 1
    f1 = 45
    f2 = 55
    # Initialise
    fs = 1000  # 1000 Hz sampling rate
    number_of_taps = 500
    # Chosen the interval 10,000-40,000 that are 20 seconds due to fs = 1000.
    data = np.loadtxt("ecg1.dat")
    time = data[10000:40000, 0]
    oringin = data[10000:40000, 1]

    Shifted = filterShift(f0,f1,f2)
    a = FIR_filter(Shifted)
    y = np.zeros(len(oringin))
    for i in range(len(oringin)):
        y[i] = a.dofilter(oringin[i])

    plt.figure(1)
    plt.plot(time, oringin)
    plt.title("Pre Filtered")
    plt.xlabel("Time (s)")
    plt.ylabel('Amplitude')

    plt.figure(2)
    plt.title("Filtered")
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.plot(y)

    print("Fertig!")
    plt.show()

