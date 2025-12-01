// MIT License
//
// Copyright (c) 2025 Jussi Lind
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef FFTBALANCE_HPP
#define FFTBALANCE_HPP

#include <fftw3.h>
#include <sndfile.h>
#include <string>
#include <vector>

struct Band
{
    double low;
    double high;
    double fCenter;
};

class FFTBalance
{
public:
    FFTBalance(int numBands, double maxBoostDb, double maxCutDb);
    ~FFTBalance();

    void loadTracks(const std::string & infile, const std::string & reffile);
    void processAndWrite(const std::string & outfile);

private:
    // Helper functions
    using SampleVector = std::vector<double>;
    struct InputData
    {
        SampleVector monoData;
        SampleVector fullBuffer;
        SF_INFO info;
    };
    InputData readFile(const std::string & filepath) const;
    std::vector<Band> generateBands(double samplerate, size_t bandCount) const;
    double getBandCenterFrequency(double low, double high) const;
    size_t freqToBin(double frequency, double samplerate, size_t binCount) const;
    SampleVector calculateBandAmp(const SampleVector & monoData, const SF_INFO & sfinfoIn, const std::vector<Band> & bands, size_t firstNonZeroFrame, size_t lastNonZeroFrame) const;
    void calculateGains();

    // Member variables
    int m_bandCount;
    double m_maxBoostDb;
    double m_maxCutDb;

    InputData m_inputData;
    InputData m_referenceData;

    std::vector<Band> m_bands;
    SampleVector m_inputBandAmp;
    SampleVector m_refBandAmp;
    SampleVector m_gains;
};

#endif // FFTBALANCE_HPP
