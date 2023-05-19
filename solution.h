//! Rohde & Schwarz Engineering Competition 2023
//!
//! This is the code to speed up. Enjoy!

#pragma once

#include "ec2023/ec2023.h"
#include <iostream>
#include <iomanip>
#include <vector>

static constexpr float OVERLAP_RATIO = 0.75;
static constexpr size_t WINDOW_SIZE = 1024;
static constexpr size_t WINDOW_SIZE_2 = 1024 * 1024;
static constexpr float log_10_window_size_2 = 6.02059991327962f;
static constexpr float log_10_4 = 0.6020599913279623f;


void compute_fourier_transform(const std::vector<ec::Float> &input, std::vector<ec::Float> &outputReal, std::vector<ec::Float> &outputImag);

std::vector<ec::Float> process_signal(const std::vector<ec::Float> &inputSignal)
{
  const size_t numSamples = inputSignal.size();
  const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
  const size_t vecHW_block_size_32 = sizeSpectrum / 32 + 1;

  const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
  const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;
  const ec::Float PI = 3.14159265358979323846f;
    const float log_10_bias = 3.0f - log_10_window_size_2 + log_10_4;
    const float log_10_bias_times_10 = 10.0f * log_10_bias;
    const float log_10_bias_for_HT = 3.0f - log_10_window_size_2;
    const float log_10_bias_for_HT_times_10 = log_10_bias_for_HT * 10.0f;

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

    ec::VecHw &vecHW2 = *ec::VecHw::getSingletonVecHw();
    vecHW2.resetMemTo0();

    int used_index = 0;

    int signalFreqReal_index[2];
    vecHW2.copyToHw(signalFreqReal, 0, WINDOW_SIZE, used_index);
    signalFreqReal_index[0] = used_index;
    used_index += WINDOW_SIZE;
    signalFreqReal_index[1] = used_index;

    int signalFreqImag_index[2];
    vecHW2.copyToHw(signalFreqImag, 0, WINDOW_SIZE, used_index);
    signalFreqImag_index[0] = used_index;
    used_index += WINDOW_SIZE;
    signalFreqImag_index[1] = used_index;

    int signalFreqReal_square_index[2];
    for (int mul_index = 0; mul_index < sizeSpectrum; mul_index += 32)
    {
      // find signalFreqReal[i] * signalFreqReal[i]
      vecHW2.mul32(signalFreqReal_index[0] + mul_index, signalFreqReal_index[0] + mul_index, used_index + mul_index);
    }
      signalFreqReal_square_index[0] = used_index;
    used_index += vecHW_block_size_32 * 32;
      signalFreqReal_square_index[1] = used_index;

    int signalFreqImag_square_index[2];
    for (int mul_index = 0; mul_index < sizeSpectrum; mul_index += 32)
    {
      // find signalFreqImag[i] * signalFreqImag[i]
      vecHW2.mul32(signalFreqImag_index[0] + mul_index, signalFreqImag_index[0] + mul_index, used_index + mul_index);
    }
      signalFreqImag_square_index[0] = used_index;
    used_index += vecHW_block_size_32 * 32;
      signalFreqImag_square_index[1] = used_index;

      int freqVal_vec_index[2];
    for (int add_index = 0; add_index < sizeSpectrum; add_index += 32)
    {
      vecHW2.add32(signalFreqReal_square_index[0] + add_index, signalFreqImag_square_index[0] + add_index, used_index+add_index);
    }
      freqVal_vec_index[0] = used_index;
      used_index += vecHW_block_size_32 * 32;
      freqVal_vec_index[1] = used_index;



    std::vector<ec::Float> freqVal_vec(sizeSpectrum, ec::Float(0));
    vecHW2.copyFromHw(freqVal_vec, freqVal_vec_index[0], sizeSpectrum, 0);


      size_t i = 0; // i = 0
      freqVal_vec[0] = 10.0f * ec::ec_log10(freqVal_vec[0]) + log_10_bias_for_HT_times_10;
      outputSpectrum[0] = ec::ec_max(outputSpectrum[0], freqVal_vec[0]);

      i = sizeSpectrum - 1; // i = 512
      freqVal_vec[512] = 10.0f * ec::ec_log10(freqVal_vec[512]) + log_10_bias_for_HT_times_10;
      outputSpectrum[512] = ec::ec_max(outputSpectrum[512], freqVal_vec[512]);



    for (i = 1; i < sizeSpectrum-1; i++)
    {
        freqVal_vec[i] = log_10_bias_times_10 + 10 * ec::ec_log10(freqVal_vec[i]);
        outputSpectrum[i] = ec::ec_max(outputSpectrum[i], freqVal_vec[i]);
    }


    idxStartWin += stepBetweenWins;
  }

  return outputSpectrum;
}

void compute_fourier_transform(const std::vector<ec::Float> &input, std::vector<ec::Float> &outputReal, std::vector<ec::Float> &outputImag)
{
  const ec::Float PI = 3.14159265358979323846f;

  size_t inputSize = input.size();

  outputReal.clear();
  outputReal.resize(inputSize, 0.0f);
  outputImag.clear();
  outputImag.resize(inputSize, 0.0f);

  for (size_t I = 0; I < inputSize; ++I)
  {
    for (size_t J = 0; J < inputSize; ++J)
    {
      const ec::Float angleTerm = (-2.0f * PI) * ec::Float(I) * J * (1.0f / ec::Float(inputSize));

      outputReal[I] += input[J] * ec_cos(angleTerm);
      outputImag[I] += input[J] * ec_sin(angleTerm);
    }
  }

  return;
}
