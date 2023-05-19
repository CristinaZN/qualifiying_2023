//! Rohde & Schwarz Engineering Competition 2023
//!
//! This is the code to speed up. Enjoy!

#pragma once

#include "ec2023/ec2023.h"
#include <iostream>      
#include <iomanip>
#include <vector>
#include <complex>
static constexpr float OVERLAP_RATIO = 0.75;
static constexpr size_t WINDOW_SIZE = 1024;

void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag);
unsigned int bitReverse(unsigned int x, int log2n) {
  int n = 0;
  int mask = 0x1;
  for (int i=0; i < log2n; i++) {
    n <<= 1;
    n |= (x & 1);
    x >>= 1;
  }
  return n;
}

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{  
  const size_t numSamples = inputSignal.size();
  const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
  const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
  const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;
  const ec::Float PI = 3.14159265358979323846f;

  std::vector<ec::Float> signalWindow(WINDOW_SIZE);
  std::vector<ec::Float> signalFreqReal(WINDOW_SIZE);
  std::vector<ec::Float> signalFreqImag(WINDOW_SIZE);
  std::vector<ec::Float> spectrumWindow(sizeSpectrum);
  std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

  size_t idxStartWin = 0;

  for (size_t J = 0; J < numWins; J++)
  {
    for (size_t I = 0; I < WINDOW_SIZE; I++)
    {
      ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));                                         
      blackmanWinCoef = blackmanWinCoef + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));

      signalWindow[I] = inputSignal[I + idxStartWin] * blackmanWinCoef;
    }

    compute_fourier_transform(signalWindow, signalFreqReal, signalFreqImag);

    for (size_t I = 0; I < sizeSpectrum; I++)
    {
      ec::Float freqVal = signalFreqReal[I] * signalFreqReal[I] + signalFreqImag[I] * signalFreqImag[I];
      freqVal = ec_sqrt(freqVal);
      freqVal = freqVal / ec::Float(WINDOW_SIZE);

      if (I > 0 && I < sizeSpectrum - 1) freqVal = freqVal * 2.0f;

      freqVal = freqVal * freqVal;

      freqVal = 10.0f * ec_log10(1000.0f * freqVal);

      outputSpectrum[I] = ec_max(outputSpectrum[I], freqVal);
    }

    idxStartWin += stepBetweenWins;

  }

  return outputSpectrum;
}

void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag)
{
  const ec::Float PI = 3.14159265358979323846f;

  size_t inputSize = input.size();

  outputReal.clear();
  outputReal.resize(inputSize, 0.0f);
  outputImag.clear();
  outputImag.resize(inputSize, 0.0f);

  int log2n=log(inputSize)/log(2);
  int n = 1 << log2n;
  //change input into complex version as x

  std::vector<std::complex<ec::Float>>x;
  for (int i=0;i<8;i++){
      x.push_back({input[i],0});
  }
  //change the order of sequences
  std::complex<ec::Float>b[inputSize];
  for (unsigned int i=0; i < n; ++i) {
    b[bitReverse(i, log2n)] = x[i];
  }

  for (int s = 1; s <= log2n; ++s) {
    int m = 1 << s;
    int m2 = m >> 1;
    std::complex<ec::Float> w(1, 0);
    std::complex<ec::Float> wm (ec_cos(-1.0f*PI / m2),ec_sin(-1.0f*PI / m2));
    for (int j=0; j < m2; ++j) {
      for (int k=j; k < n; k += m) {
        std::complex<ec::Float>  t = w * b[k + m2];
        std::complex<ec::Float>  u = b[k];
        b[k] = u + t;
        b[k + m2] = u - t;
      }
     w *= wm;
    }
  }
  
  for(int i=0;i<inputSize;i++){
    outputReal[i]=b[i].real();
    outputImag[i]=b[i].imag();
  }
  return;
}

