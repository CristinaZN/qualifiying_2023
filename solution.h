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
size_t opt_num = 32;
size_t cos_opt_num = 4;
void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag);

std::vector<ec::Float> valueVector(ec::Float number, size_t size){
  std::vector<ec::Float> *rlt = new std::vector<ec::Float>(size,number);
  return *rlt;
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

  // Copy input signal to the HW
  ec::VecHw& hwInputSingal= *ec::VecHw::getSingletonVecHw();
  hwInputSingal.copyToHw(inputSignal,0,inputSignal.size(),0);

  //Initialize 
  // ec::VecHw& signalWindow = *ec::VecHw::getSingletonVecHw();
  // ec::VecHw& signalFreqReal = *ec::VecHw::getSingletonVecHw();
  // ec::VecHw& signalFreqImag = *ec::VecHw::getSingletonVecHw();
  // ec::VecHw& spectrumWindow = *ec::VecHw::getSingletonVecHw();
  // ec::VecHw& outputSpectrum = *ec::VecHw::getSingletonVecHw();

  size_t idxStartWin = 0;
  size_t step = WINDOW_SIZE;

  // Initilize a vector I = [0 , 1 , 2 , 3, ...., WINDOW_SIZE - 1]
  std::vector<ec::Float> vecI(WINDOW_SIZE*4,0);
  for(size_t i = 0 ; i < WINDOW_SIZE; i++){
    vecI[i] = i;}
  // Initialize a vector A = [a, a, ..., a]
  ec::Float constant = 2.0f * PI / (WINDOW_SIZE -1);//这一步应该没必要用vechw？
  std::vector<ec::Float> vecA = valueVector(constant,WINDOW_SIZE);

  ec::VecHw& hwI = *ec::VecHw::getSingletonVecHw();
  // hwI = [vecI, vecA, 0, 0]
  hwI.copyToHw(vecI,0,WINDOW_SIZE,0);
  hwI.copyToHw(vecA,0,WINDOW_SIZE,WINDOW_SIZE);
  // compute (ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1)), copy to hwI[WinSize+1:2*WinSize]
  // hwI = [vecI, result, 0, 0]
  for (size_t i = 0; i < WINDOW_SIZE / opt_num; i++){
    hwI.mul32(i*opt_num,WINDOW_SIZE+i*opt_num,WINDOW_SIZE+i*opt_num,opt_num);
  }
  // compute cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1)), copy to hwI[WinSize+1:2*WinSize]
  // hwI = [vecI, result, 0, 0]
  for (size_t i = 0; i < WINDOW_SIZE / cos_opt_num; i++){
    hwI.cos4(WINDOW_SIZE+i*cos_opt_num,WINDOW_SIZE+i*cos_opt_num,cos_opt_num);
  }
  // compute - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1)), copy to hwI[WinSize+1:2*WinSize]
  // hwI = [vector(-0.5), result, 0, 0]
  constant = -0.5f;
  vecA = valueVector(constant,WINDOW_SIZE);
  hwI.copyToHw(vecA,0,WINDOW_SIZE,0);
  for (size_t i = 0; i < WINDOW_SIZE / opt_num; i++){
    hwI.mul32(i*opt_num,WINDOW_SIZE+i*opt_num,WINDOW_SIZE+i*opt_num,opt_num);
  }
  // hwI.mul32(0,WINDOW_SIZE,WINDOW_SIZE,WINDOW_SIZE);
  // compute 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1)), copy to hwI[WinSize+1:2*WinSize]
  // hwI = [vector(0.42), result, 0, 0]
  constant = 0.42f;
  vecA = valueVector(constant,WINDOW_SIZE);
  hwI.copyToHw(vecA,0,WINDOW_SIZE,0);
  for (size_t i = 0; i < WINDOW_SIZE / opt_num; i++){
    hwI.add32(i*opt_num,WINDOW_SIZE+i*opt_num,WINDOW_SIZE+i*opt_num,opt_num);
  }
  // hwI.add32(0,WINDOW_SIZE,WINDOW_SIZE,WINDOW_SIZE);
  // hwI = [vector(0.42), result1, 0, result1]
  // result1 = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1))
  for (size_t i = 0 ; i< WINDOW_SIZE / opt_num; i++){
    hwI.assign32(WINDOW_SIZE+i*opt_num,WINDOW_SIZE*3 + i * opt_num, opt_num);
  }
  // hwI.assign32(WINDOW_SIZE,WINDOW_SIZE*3,WINDOW_SIZE);


  constant = 4.0f * PI / (WINDOW_SIZE -1);
  vecA = valueVector(constant,WINDOW_SIZE);
  // hwI = [vecI,vecA,0,result1]
  hwI.copyToHw(vecI,0,WINDOW_SIZE,0);
  hwI.copyToHw(vecA,0,WINDOW_SIZE,WINDOW_SIZE);
  //compute (ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1)), copy to hwI[WinSize+1:2*WinSize]
  // hwI = [vecI,result(vecI.*vecA),0,result1]
  for (size_t i = 0; i < WINDOW_SIZE / opt_num; i++){
    hwI.mul32(i*opt_num,WINDOW_SIZE+i*opt_num,WINDOW_SIZE+i*opt_num,opt_num);
  }
  // hwI.mul32(0,WINDOW_SIZE,WINDOW_SIZE,WINDOW_SIZE);
  //compute cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1)), copy to hwI[WinSize+1:2*WinSize]
  // hwI = [vecI,result(cos(vecI*vecA)),0,result1]
  for (size_t i = 0; i < WINDOW_SIZE / cos_opt_num; i++){
    hwI.cos4(WINDOW_SIZE+i*cos_opt_num,WINDOW_SIZE+i*cos_opt_num,cos_opt_num);
  }
  // hwI.cos4(WINDOW_SIZE,WINDOW_SIZE,WINDOW_SIZE);
  //compute 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1)), copy to hwI[WinSize+1:2*WinSize]
  constant = 0.08f;
  vecA = valueVector(constant,WINDOW_SIZE);
  hwI.copyToHw(vecA,0,WINDOW_SIZE,0);
  // hwI = [vector(0.08),result(0.08*cos(vecI*vecA)),0,result1]
  for (size_t i = 0; i < WINDOW_SIZE / opt_num; i++){
    hwI.mul32(i*opt_num,WINDOW_SIZE+i*opt_num,WINDOW_SIZE+i*opt_num,opt_num);
  }
  // hwI.mul32(0,WINDOW_SIZE,WINDOW_SIZE,WINDOW_SIZE);
  // compute blackmanWinCoef + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));
  // hwI = [final_result,result(0.08*cos(vecI*vecA)),0,result1]
  // final_result = result1 + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1))
  for (size_t i = 0; i < WINDOW_SIZE / opt_num; i++){
    hwI.add32(WINDOW_SIZE + i * opt_num,WINDOW_SIZE * 3 + i * opt_num,i*opt_num,opt_num);
  }
  // hwI.add32(WINDOW_SIZE,WINDOW_SIZE*3,0,WINDOW_SIZE);
  // ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));                                         
  // blackmanWinCoef = blackmanWinCoef + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));


  for (size_t J = 0; J < numWins; J++)
  {
    
    // for (size_t I = 0; I < WINDOW_SIZE; I++)
    // {
      // ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));                                         
      // blackmanWinCoef = blackmanWinCoef + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));
      // signalWindow[I] = inputSignal[I + idxStartWin] * blackmanWinCoef;
    // }
    
    // hwI = [final_result,signalWindow[idxStartWin:idxStartWin+WinSize-1],0,result1]
    hwI.copyToHw(inputSignal,idxStartWin,WINDOW_SIZE,WINDOW_SIZE);
    // hwI = [final_result,signalWindow[idxStartWin:idxStartWin+WinSize-1],result,result1]
    for(size_t i = 0 ; i < WINDOW_SIZE/opt_num ; i++){
      hwI.mul32(i * opt_num, WINDOW_SIZE + i * opt_num, WINDOW_SIZE*2+i*opt_num,opt_num);
    }
    // hwI.mul32(0,WINDOW_SIZE,WINDOW_SIZE*2,WINDOW_SIZE);
    hwI.copyFromHw(signalWindow,WINDOW_SIZE*2,WINDOW_SIZE,0);
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

