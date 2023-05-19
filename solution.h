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
static constexpr size_t WINDOW_SIZE_2 = 1024 * 1024;
static constexpr float log_10_window_size_2 = 6.02059991327962f;
static constexpr float log_10_4 = 0.6020599913279623f;
void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag);


std::vector<ec::Float> valueVector(ec::Float number, size_t size){
  std::vector<ec::Float> *rlt = new std::vector<ec::Float>(size,number);
  return *rlt;
}

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{  

  const size_t numSamples = inputSignal.size();
  const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
  const size_t vecHW_block_size_32 = sizeSpectrum / 32 + 1; // vecHW calculation, 32 elements per block

  const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
  const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;
  const ec::Float PI = 3.14159265358979323846f;

  // after reformulation of the expression, these params are used directly
    const float log_10_bias = 3.0f - log_10_window_size_2 + log_10_4;
    const float log_10_bias_times_10 = 10.0f * log_10_bias;
    const float log_10_bias_for_HT = 3.0f - log_10_window_size_2;
    const float log_10_bias_for_HT_times_10 = log_10_bias_for_HT * 10.0f;

  std::vector<ec::Float> signalWindow(WINDOW_SIZE);
  std::vector<ec::Float> signalFreqReal(WINDOW_SIZE);
  std::vector<ec::Float> signalFreqImag(WINDOW_SIZE);
  std::vector<ec::Float> spectrumWindow(sizeSpectrum);
  std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());
    ec::VecHw& hwInputSingal= *ec::VecHw::getSingletonVecHw();
    hwInputSingal.copyToHw(inputSignal,0,inputSignal.size(),0);
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
    hwI.resetMemTo0(0,WINDOW_SIZE*4);
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
    std::vector<ec::Float> vec1(WINDOW_SIZE,0);
    hwI.copyFromHw(vec1,0,WINDOW_SIZE,0);
  
  for (size_t J = 0; J < numWins; J++)
  {

    
     for (size_t I = 0; I < WINDOW_SIZE; I++)
     {
       ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));
       blackmanWinCoef = blackmanWinCoef + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));
       signalWindow[I] = inputSignal[I + idxStartWin] * blackmanWinCoef;
     }


    // hwI = [final_result,signalWindow[idxStartWin:idxStartWin+WinSize-1],0,result1]
    hwI.copyToHw(inputSignal,idxStartWin,WINDOW_SIZE,WINDOW_SIZE);
    hwI.copyToHw(vec1,0,WINDOW_SIZE,0);
    // hwI = [final_result,signalWindow[idxStartWin:idxStartWin+WinSize-1],result,result1]
    for(size_t i = 0 ; i < WINDOW_SIZE/opt_num ; i++) {
        hwI.mul32(i * opt_num, WINDOW_SIZE + i * opt_num, WINDOW_SIZE * 2 + i * opt_num, opt_num);
    }
    
  for (size_t I = 0; I < WINDOW_SIZE; I++)
     {
       ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));
       blackmanWinCoef = blackmanWinCoef + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));
        // std::cout <<vec[I].toFloat()<<' '<<blackmanWinCoef.toFloat()<<std::endl;
     }
    // hwI.mul32(0,WINDOW_SIZE,WINDOW_SIZE*2,WINDOW_SIZE);
      std::vector<ec::Float> signalWindow_dan(WINDOW_SIZE);
    hwI.copyFromHw(signalWindow_dan,WINDOW_SIZE*2,WINDOW_SIZE,0);

   for(size_t i = 0; i < WINDOW_SIZE; i++){
       if(signalWindow_dan[i] - signalWindow[i] > ec::Float(0.0000001f)
       || signalWindow_dan[i] - signalWindow[i] < ec::Float(0.0000001f)) {
           std::cout << "signalWindow_dan[i] != signalWindow[i] at: " << i <<std::endl;
           printf("%lf %lf\n",signalWindow_dan[i].toFloat(),signalWindow[i].toFloat());
       }
   }


    compute_fourier_transform(signalWindow_dan, signalFreqReal, signalFreqImag);

    int memory_used = WINDOW_SIZE;


    ec::VecHw &vecHW2 = *ec::VecHw::getSingletonVecHw();
    vecHW2.resetMemTo0();

    int used_index = 0; // used for counting index

    int signalFreqReal_index[2]; // store the index of Sig_Re in HW_mem
    vecHW2.copyToHw(signalFreqReal, 0, sizeSpectrum, used_index);
    signalFreqReal_index[0] = used_index;
    used_index += sizeSpectrum;
    signalFreqReal_index[1] = used_index;

    int signalFreqImag_index[2];  // store the index of Sig_Im in HW_mem
    vecHW2.copyToHw(signalFreqImag, 0, sizeSpectrum, used_index);
    signalFreqImag_index[0] = used_index;
    used_index += sizeSpectrum;
    signalFreqImag_index[1] = used_index;


    int signalFreqReal_square_index[2]; // store the index of Sig_Re^2 in HW_mem
    for (int mul_index = 0; mul_index < sizeSpectrum; mul_index += 32)
    {
      // find signalFreqReal[i] * signalFreqReal[i]
      vecHW2.mul32(signalFreqReal_index[0] + mul_index, signalFreqReal_index[0] + mul_index, used_index + mul_index);
    }
    signalFreqReal_square_index[0] = used_index;
    used_index += vecHW_block_size_32 * 32;
    signalFreqReal_square_index[1] = used_index;

    int signalFreqImag_square_index[2]; // store the index of Sig_Im^2 in HW_mem
    for (int mul_index = 0; mul_index < sizeSpectrum; mul_index += 32)
    {
      // find signalFreqImag[i] * signalFreqImag[i]
      vecHW2.mul32(signalFreqImag_index[0] + mul_index, signalFreqImag_index[0] + mul_index, used_index + mul_index);
    }


    signalFreqImag_square_index[0] = used_index;
    used_index += vecHW_block_size_32 * 32;
    signalFreqImag_square_index[1] = used_index;

    int freqVal_vec_index[2]; // store the index of freqVal in HW_mem
    for (int add_index = 0; add_index < sizeSpectrum; add_index += 32)
    {
      vecHW2.add32(signalFreqReal_square_index[0] + add_index, signalFreqImag_square_index[0] + add_index, used_index+add_index);
    }
      freqVal_vec_index[0] = used_index;
      used_index += vecHW_block_size_32 * 32;
      freqVal_vec_index[1] = used_index;



    std::vector<ec::Float> freqVal_vec(sizeSpectrum, ec::Float(0));
    vecHW2.copyFromHw(freqVal_vec, freqVal_vec_index[0], sizeSpectrum, 0);

    // Head and Tail is not mul_by_2
      size_t i = 0; // i = 0
      freqVal_vec[0] = 10.0f * ec::ec_log10(freqVal_vec[0]) + log_10_bias_for_HT_times_10;
      outputSpectrum[0] = ec::ec_max(outputSpectrum[0], freqVal_vec[0]);

      i = sizeSpectrum - 1; // i = 512
      freqVal_vec[512] = 10.0f * ec::ec_log10(freqVal_vec[512]) + log_10_bias_for_HT_times_10;
      outputSpectrum[512] = ec::ec_max(outputSpectrum[512], freqVal_vec[512]);


    // cancel square_root and square
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
