//! Rohde & Schwarz Engineering Competition 2023
//!
//! (c) Rohde & Schwarz GmbH & Co. KG, Munich
//!
//! This code is needed to run your solution locally. Any changes will be lost once you upload your solution!

#pragma once
#include <vector>

static constexpr size_t VEC_HW_MEM_SIZE = 4 * 1024; // [4B-values]

namespace ec
{

  class VecHw final
  {
  private:
    std::vector<float> m_mem; // cannot be vector<ec::Float>, because then all ops, even initialisation, are measured individually.

    static inline VecHw* pSingletonVecHw = nullptr;

    inline VecHw()
      : m_mem(VEC_HW_MEM_SIZE, 0.0f)
    {
    }

    // VecHw object不可以复制
    VecHw(const VecHw& other) = delete;

  public:
    // 创建一个VecHw object
    static inline VecHw* getSingletonVecHw()
    {
      if (pSingletonVecHw == nullptr)
      {
        pSingletonVecHw = new VecHw;
      }

      return  pSingletonVecHw;
    }

    /* 把m_mem这个向量的某一部分初始化为零向量
     * index_start: 从m_mem向量的第index_start个位置开始
     * num: 将之后的num个元素设置为零
     */
    inline void resetMemTo0(size_t idx_start, size_t num)
    {
      if (idx_start + num > m_mem.size())
        throw std::runtime_error("Size of memory of VecHw to reset does not fit the size of the memory.");

      std::fill(m_mem.begin() + static_cast<long int>(idx_start), m_mem.begin() + static_cast<long int>(idx_start + num), 0.0f);
    }


    inline void resetMemTo0()
    {
      resetMemTo0(0, m_mem.size());
    }

    /* 将操作数拷贝到VecHw object的m_mem向量里面
     * data: 外部的操作数，即要被拷贝的向量
     * idx_from: 从data的哪个位置开始拷贝
     * num: 拷贝多少个data的元素
     * idx_to: 拷贝到m_mem的哪里
     */
    inline void copyToHw(const std::vector<ec::Float>& data, size_t idx_from, size_t num, size_t idx_to)
    {
      Measurement::inc(Measurement::MeasType::initCopy_VecHw);

      if (idx_from + num > data.size())
        throw std::runtime_error("Size of vector to copy does not fit the number of elements to be copied.");

      if (idx_to + num > m_mem.size())
        throw std::runtime_error("Size of memory of VecHw to copy to does not fit the number of elements to be copied.");

      for (size_t I = 0; I < num; I++)
      {
        Measurement::inc(Measurement::MeasType::copy4B_VecHw);
        m_mem[idx_to + I] = data[idx_from + I].m_value;
      }

    }

    /* 把结果拷贝到外部
     * data: 外面的向量，用来保存结果
     * idx_from: 从m_mem的哪个位置开始拷贝
     * num: 拷贝多少个m_mem的元素
     * idx_to: 拷贝到data的哪里
     */
    inline void copyFromHw(std::vector<ec::Float>& data, size_t idx_from, size_t num, size_t idx_to)
    {
      Measurement::inc(Measurement::MeasType::initCopy_VecHw);

      if (idx_from + num > m_mem.size())
        throw std::runtime_error("Size of memory of VecHw to copy from does not fit the elements to be copied.");

      if (idx_to + num > data.size())
        throw std::runtime_error("Size of vector to copy to does not fit the elements to be copied.");

      for (size_t I = 0; I < num; I++)
      {
        Measurement::inc(Measurement::MeasType::copy4B_VecHw);
        data[idx_to + I] = ec::Float(m_mem[idx_from + I]); // Here is added 1x ec::Float assignment to the score additionaly.
      }
    }

    /* 把m_mem的一部分元素拷贝到m_mem的另一部分（相当于平移）
     * idx_from: 要被拷贝的m_mem元素 的位置
     * idx_to: 拷贝到哪里
     * num: 拷贝多少（默认是32个位置）
     *
     */
    inline void assign32(size_t idx_from, size_t idx_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::assign32_VecHw);

      if (idx_from + num > m_mem.size())
        throw std::runtime_error("Size of memory of VecHw to assign from does not fit the elements to be copied.");
      if (idx_to + num > m_mem.size())
        throw std::runtime_error("Size of memory of VecHw to assign to does not fit the elements to be copied.");
      if (num > 32)
        throw std::runtime_error("assign32 can process only max 32 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[idx_to + I] = m_mem[idx_from + I];
      }
    }

    inline void assign32(size_t idx_from, size_t idx_to)
    {
      assign32(idx_from, idx_to, 32ull);
    }

    /* 将m_mem的一部分和m_mem的另一部分相乘，结果存储在m_mem的另一部分（这里是元素对应相乘，结果还是一个向量）
     * Matlab里这个乘法表示为 C = A.*B （A，B，C是长度相等的向量）
     * left_from: A开始的位置
     * right_from: B开始的位置
     * out_to: C开始的位置
     * num: 向量长度（默认是32）
     */
    inline void mul32(size_t left_from, size_t right_from, size_t out_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::mul32_VecHw);

      if (left_from + num > m_mem.size())
        throw std::runtime_error("Left operand of mul32 is out of bound of memory size.");
      if (right_from + num > m_mem.size())
        throw std::runtime_error("Right operand of mul32 is out of bound of memory size.");
      if (out_to + num > m_mem.size())
        throw std::runtime_error("Result of mul32 is out of bound of memory size.");
      if (num > 32)
        throw std::runtime_error("mul32 can process only max 32 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_to + I] = m_mem[left_from + I] * m_mem[right_from + I];
      }
    }

    inline void mul32(size_t left_from, size_t right_from, size_t out_to)
    {
      mul32(left_from, right_from, out_to, 32ull);
    }


    inline void mul32(size_t left_from, ec::Float right_scalar, size_t out_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::mul32_VecHw);

      if (left_from + num > m_mem.size())
        throw std::runtime_error("Left operand of mul32 is out of bound of memory size.");
      if (out_to + num > m_mem.size())
        throw std::runtime_error("Result of mul32 is out of bound of memory size.");
      if (num > 32)
        throw std::runtime_error("mul32 can process only max 32 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_to + I] = m_mem[left_from + I] * right_scalar.m_value;
      }
    }

    inline void mul32(size_t left_from, ec::Float right_scalar, size_t out_to)
    {
      mul32(left_from, right_scalar, out_to, 32ull);
    }
    /* 将m_mem的一部分和m_mem的另一部分相加，结果存储在m_mem的指定位置
     * 假定加法操作是 C = A + B（A，B，C都是向量）
     * left_from: 左加数起始位置，相当于A向量开始的位置
     * right_from: 右加数起始位置，相当于B向量开始的位置
     * out_to: 存放结果的起始位置，相当于C向量开始的位置
     * num: 向量长度，或者说是执行多少次对应相加，默认为32
    */
    inline void add32(size_t left_from, size_t right_from, size_t out_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::add32_VecHw);

      if (left_from + num > m_mem.size())
        throw std::runtime_error("Left operand of add32 is out of bound of memory size.");
      if (right_from + num > m_mem.size())
        throw std::runtime_error("Right operand of add32 is out of bound of memory size.");
      if (out_to + num > m_mem.size())
        throw std::runtime_error("Result of add32 is out of bound of memory size.");
      if (num > 32)
        throw std::runtime_error("add32 can process only max 32 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_to + I] = m_mem[left_from + I] + m_mem[right_from + I];
      }
    }

    inline void add32(size_t left_from, size_t right_from, size_t out_to)
    {
      add32(left_from, right_from, out_to, 32ull);
    }

    inline void add32(size_t left_from, ec::Float right_scalar, size_t out_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::add32_VecHw);

      if (left_from + num > m_mem.size())
        throw std::runtime_error("Left operand of add32 is out of bound of memory size.");
      if (out_to + num > m_mem.size())
        throw std::runtime_error("Result of add32 is out of bound of memory size.");
      if (num > 32)
        throw std::runtime_error("add32 can process only max 32 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_to + I] = m_mem[left_from + I] + right_scalar.m_value;
      }
    }

    inline void add32(size_t left_from, ec::Float right_scalar, size_t out_to)
    {
      add32(left_from, right_scalar, out_to, 32ull);
    }

    /*
    * ===============================================
    * 执行向量相加，将另一个向量直接加到指定向量上
    * B=B+A (原来的B向量直接被结果覆盖了)
    * in_from: A向量的起始位置，不会被覆盖
    * out_scalar_to: B向量的起始位置，后面会直接被覆盖掉
    * num: 向量长度
    * ===============================================
    * 示例：
    * 执行前: m_mem: [1 1 1 2 2 2 3 3 3]
    * in_from 0: A = [1 1 1]
    * out_scalar_to 4: B = [2 2 3]
    * num 3
    * 相当于第一个向量是223，第二个向量是111
    * 将223加111，再把结果直接覆盖到223上
    * 结果: m_mem: [1 1 1 2 3 3 4 3 3]
    */
    inline void acc32(size_t in_from, size_t out_scalar_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::acc32_VecHw);

      if (in_from + num > m_mem.size())
        throw std::runtime_error("Left operand of acc32 is out of bound of memory size.");
      if (out_scalar_to + 1 > m_mem.size())
        throw std::runtime_error("Result of acc32 is out of bound of memory size.");
      if (num > 32)
        throw std::runtime_error("acc32 can process only max 32 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_scalar_to] += m_mem[in_from + I];
      }
    }

    inline void acc32(size_t in_from, size_t out_scalar_to)
    {
      acc32(in_from, out_scalar_to, 32ull);
    }

/* 将m_mem的一部分求sin，存到m_mem的另一部分
* 类似于 B=sin(A)
* in_from: A开始的位置
* out_to: B开始的位置
* num: 向量A，B的长度（默认是4）
*/
    inline void sin4(size_t in_from, size_t out_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::sin4_VecHw);

      if (in_from + 4 > m_mem.size())
        throw std::runtime_error("Left operand of sin4 is out of bound of memory size.");
      if (out_to + 4 > m_mem.size())
        throw std::runtime_error("Result of sin4 is out of bound of memory size.");
      if (num > 4)
        throw std::runtime_error("sin4 can process only max 4 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_to + I] = sinf(m_mem[in_from + I]);
      }
    }

    inline void sin4(size_t in_from, size_t out_to)
    {
      sin4(in_from, out_to, 4ull);
    }

/* 将m_mem的一部分求cos，存到m_mem的另一部分
* 类似于 B=cos(A)
* in_from: A开始的位置
* out_to: B开始的位置
* num: 向量A，B的长度（默认是4）
*/
    inline void cos4(size_t in_from, size_t out_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::cos4_VecHw);

      if (in_from + 4 > m_mem.size())
        throw std::runtime_error("Left operand of cos4 is out of bound of memory size.");
      if (out_to + 4 > m_mem.size())
        throw std::runtime_error("Result of cos4 is out of bound of memory size.");
      if (num > 4)
        throw std::runtime_error("cos4 can process only max 4 values in parallel.");


      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_to + I] = cosf(m_mem[in_from + I]);
      }
    }

    inline void cos4(size_t in_from, size_t out_to)
    {
      cos4(in_from, out_to, 4ull);
    }

/* 将m_mem的一部分开方，存到m_mem的另一部分
* 类似于 B=sqrt(A)
* in_from: A开始的位置
* out_to: B开始的位置
* num: 向量A，B的长度（默认是4）
*/
    inline void sqrt4(size_t in_from, size_t out_to, size_t num)
    {
      Measurement::inc(Measurement::MeasType::sqrt4_VecHw);

      if (in_from + 4 > m_mem.size())
        throw std::runtime_error("Left operand of sqrt4 is out of bound of memory size.");
      if (out_to + 4 > m_mem.size())
        throw std::runtime_error("Result of sqrt4 is out of bound of memory size.");
      if (num > 4)
        throw std::runtime_error("sqrt4 can process only max 4 values in parallel.");

      for (size_t I = 0; I < num; I++)
      {
        m_mem[out_to + I] = sqrtf(m_mem[in_from + I]);
      }
    }

    inline void sqrt4(size_t in_from, size_t out_to)
    {
      sqrt4(in_from, out_to, 4ull);
    }

  };

} // namespace ec