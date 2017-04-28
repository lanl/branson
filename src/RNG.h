//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RNG.h
 * \author Alex Long
 * \date   April 27 2017
 * \brief  RNG modified from the Draco project
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//

// Copyright from Draco project (https://github.com/lanl/Draco) as follows:

//Copyright (C) 2016-2017 Los Alamos National Security, LLC.
//All rights reserved.

//Copyright 2016.  Los Alamos National Security, LLC. This software was 
//produced under U.S. Government contract DE-AC52-06NA25396 for Los 
//Alamos National Laboratory (LANL), which is operated by Los Alamos
//National Security, LLC for the U.S. Department of Energy. The 
//U.S. Government has rights to use, reproduce, and distribute this
//software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
//LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
//FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
//derivative works, such modified software should be clearly marked, so
//as not to confuse it with the version available from LANL.

//Additionally, redistribution and use in source and binary forms, with
//or without modification, are permitted provided that the following
//conditions are met:

//- Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.

//- Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the 
//  documentation and/or other materials provided with the distribution.

//- Neither the name of Los Alamos National Security, LLC, Los Alamos
//  National Laboratory, LANL, the U.S. Government, nor the names of its 
//  contributors may be used to endorse or promote products derived from
//  this software without specific prior written permission.

//THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND 
//CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
//FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS 
//ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY 
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
//IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
//OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
//ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef RNG_h
#define RNG_h

#include <stdlib.h>
#include "random123/threefry.h"

// Select a particular counter-based random number generator from Random123.
typedef r123::Threefry2x64 CBRNG;

// Counter and key types.
typedef CBRNG::ctr_type ctr_type;
typedef CBRNG::key_type key_type;


#include "random123/features/compilerfeatures.h"
#include <limits>
#include <type_traits>


namespace r123{

#if R123_USE_CXX11_TYPE_TRAITS
using std::make_signed;
using std::make_unsigned;
#else
// Sigh... We could try to find another <type_traits>, e.g., from
// boost or TR1.  Or we can do it ourselves in the r123 namespace.
// It's not clear which will cause less headache...
template <typename T> struct make_signed{};
template <typename T> struct make_unsigned{};
#define R123_MK_SIGNED_UNSIGNED(ST, UT)                 \
template<> struct make_signed<ST>{ typedef ST type; }; \
template<> struct make_signed<UT>{ typedef ST type; }; \
template<> struct make_unsigned<ST>{ typedef UT type; }; \
template<> struct make_unsigned<UT>{ typedef UT type; }

R123_MK_SIGNED_UNSIGNED(int8_t, uint8_t);
R123_MK_SIGNED_UNSIGNED(int16_t, uint16_t);
R123_MK_SIGNED_UNSIGNED(int32_t, uint32_t);
R123_MK_SIGNED_UNSIGNED(int64_t, uint64_t);
#if R123_USE_GNU_UINT128
R123_MK_SIGNED_UNSIGNED(__int128_t, __uint128_t);
#endif
#undef R123_MK_SIGNED_UNSIGNED
#endif

#if defined(__CUDACC__) || defined(_LIBCPP_HAS_NO_CONSTEXPR)
// Amazing! cuda thinks numeric_limits::max() is a __host__ function, so
// we can't use it in a device function.  
//
// The LIBCPP_HAS_NO_CONSTEXP test catches situations where the libc++
// library thinks that the compiler doesn't support constexpr, but we
// think it does.  As a consequence, the library declares
// numeric_limits::max without constexpr.  This workaround should only
// affect a narrow range of compiler/library pairings.
// 
// In both cases, we find max() by computing ~(unsigned)0 right-shifted
// by is_signed.
template <typename T>
R123_CONSTEXPR R123_STATIC_INLINE R123_CUDA_DEVICE T maxTvalue(){
    typedef typename make_unsigned<T>::type uT;
    return (~uT(0)) >> std::numeric_limits<T>::is_signed;
 }
#else
template <typename T>
R123_CONSTEXPR R123_STATIC_INLINE T maxTvalue(){
    return std::numeric_limits<T>::max();
}
#endif

// u01: Input is a W-bit integer (signed or unsigned).  It is cast to
//   a W-bit unsigned integer, multiplied by Ftype(2^-W) and added to
//   Ftype(2^(-W-1)).  A good compiler should optimize it down to an
//   int-to-float conversion followed by a multiply and an add, which
//   might be fused, dependingon the architecture.
//
//  If the input is a uniformly distributed integer, then the
//  result is a uniformly distributed floating point number in [0, 1].
//  The result is never exactly 0.0.  
//  The smallest value returned is 2^-W.
//  Let M be the number of mantissa bits in Ftype.
//  If W>M  then the largest value retured is 1.0.
//  If W<=M then the largest value returned is the largest Ftype less than 1.0.
template <typename Ftype, typename Itype>
R123_CUDA_DEVICE R123_STATIC_INLINE Ftype u01(Itype in){
    typedef typename make_unsigned<Itype>::type Utype;
    R123_CONSTEXPR Ftype factor = Ftype(1.)/(maxTvalue<Utype>() + Ftype(1.));
    R123_CONSTEXPR Ftype halffactor = Ftype(0.5)*factor;
#if R123_UNIFORM_FLOAT_STORE
    volatile Ftype x = Utype(in)*factor; return x+halffactor;
#else
    return Utype(in)*factor + halffactor;
#endif
}

// uneg11: Input is a W-bit integer (signed or unsigned).  It is cast
//    to a W-bit signed integer, multiplied by Ftype(2^-(W-1)) and
//    then added to Ftype(2^(-W-2)).  A good compiler should optimize
//    it down to an int-to-float conversion followed by a multiply and
//    an add, which might be fused, depending on the architecture.
//
//  If the input is a uniformly distributed integer, then the
//  output is a uniformly distributed floating point number in [-1, 1].
//  The result is never exactly 0.0.
//  The smallest absolute value returned is 2^-(W-1)
//  Let M be the number of mantissa bits in Ftype.
//  If W>M  then the largest value retured is 1.0 and the smallest is -1.0.
//  If W<=M then the largest value returned is the largest Ftype less than 1.0
//    and the smallest value returned is the smallest Ftype greater than -1.0.
template <typename Ftype, typename Itype>
R123_CUDA_DEVICE R123_STATIC_INLINE Ftype uneg11(Itype in){
    typedef typename make_signed<Itype>::type Stype;
    R123_CONSTEXPR Ftype factor = Ftype(1.)/(maxTvalue<Stype>() + Ftype(1.));
    R123_CONSTEXPR Ftype halffactor = Ftype(0.5)*factor;
#if R123_UNIFORM_FLOAT_STORE
    volatile Ftype x = Stype(in)*factor; return x+halffactor;
#else
    return Stype(in)*factor + halffactor;
#endif
}

// u01fixedpt:  Return a "fixed point" number in (0,1).  Let:
//   W = width of Itype, e.g., 32 or 64, regardless of signedness.
//   M = mantissa bits of Ftype, e.g., 24, 53 or 64
//   B = min(M, W)
// Then the 2^B-1 possible output values are:
//    2^-B*{1, 3, 5, ..., 2^B - 1}
// The smallest output is: 2^-B
// The largest output is:  1 - 2^-B
// The output is never exactly 0.0, nor 0.5, nor 1.0.
// The 2^(B-1) possible outputs:
//   - are equally likely,
//   - are uniformly spaced by 2^-(B-1),
//   - are balanced around 0.5
template <typename Ftype, typename Itype>
R123_CUDA_DEVICE R123_STATIC_INLINE Ftype u01fixedpt(Itype in){
    typedef typename make_unsigned<Itype>::type Utype;
    R123_CONSTEXPR int excess = std::numeric_limits<Utype>::digits - std::numeric_limits<Ftype>::digits;
    if(excess>=0)
    {

// 2015-09-26 KT - Suppress warnings for the following expressions (see https://rtt.lanl.gov/redmine/issues/416)
//
// Basically, GCC under BullseyeCoverage issues the following warning every time this file is included:
//         
// Counter_RNG.hh:124:65:   required from here
// uniform.hpp:200:48: warning: second operand of conditional expression has no effect [-Wunused-value]
//         R123_CONSTEXPR int ex_nowarn = (excess>=0) ? excess : 0;
//
// Unfortunately, if this expression is simplified (see r7628) some compilers will not compile the code because the RHS of the
// assignment may contain values that are not known at comile time (not constexpr).  We don't want to spend to much time debugging
// this issue because the code is essentially vendor owned (Random123).

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#endif
        
        R123_CONSTEXPR int ex_nowarn = (excess>=0) ? excess : 0;

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

        R123_CONSTEXPR Ftype factor = Ftype(1.)/(Ftype(1.) + ((maxTvalue<Utype>()>>ex_nowarn)));
        return (1 | (Utype(in)>>ex_nowarn)) * factor;
    }
    else
    {
        return u01<Ftype>(in);
    }
}

} // namespace r123


namespace // anonymous
{

// Forward declaration.
class Counter_RNG;

// Select a particular counter-based random number generator from Random123.
typedef r123::Threefry2x64 CBRNG;

// Counter and key types.
typedef CBRNG::ctr_type ctr_type;
typedef CBRNG::key_type key_type;

#define CBRNG_DATA_SIZE 4
//----------------------------------------------------------------------------//
/*! \brief Generate a random double.
 *
 * Given a pointer to RNG state data, this function returns a random double in
 * the open interval (0, 1)---i.e., excluding the endpoints.
 */
static inline double _ran(ctr_type::value_type *const data) {
  CBRNG rng;

  // Assemble a counter from the first two elements in data.
  ctr_type ctr = {{data[0], data[1]}};

  // Assemble a key from the last two elements in data.
  const key_type key = {{data[2], data[3]}};

  // Invoke the counter-based rng.
  const ctr_type result = rng(ctr, key);

  // Increment the counter.
  ctr.incr();

  // Copy the updated counter back into data.
  data[0] = ctr[0];
  data[1] = ctr[1];

  // Convert the first 64 bits of the RNG output into a double-precision value
  // in the open interval (0, 1) and return it.
  return r123::u01fixedpt<double, ctr_type::value_type>(result[0]);
}

} // end anonymous


//===========================================================================//
/*!
 * \class RNG
 * \brief A counter-based random-number generator.
 *
 * RNG provides an interface to a counter-based random number generator
 * from the Random123 library from D. E. Shaw Research
 * (http://www.deshawresearch.com/resources_random123.html).
 *
 * Similarly, Rnd_Control is a friend of RNG because initializing a
 * generator requires access to private data that should not be exposed through
 * the public interface.  Rnd_Control takes no responsibility for instantiating
 * RNGs itself, and since copying RNGs is disabled (via a
 * private copy constructor), an Rnd_Control must be able to initialize a
 * generator that was instantiated outside of its control.
 */
//===========================================================================//
class RNG {

public:
  typedef ctr_type::const_iterator const_iterator;

  /*! \brief Default constructor.
   *
   * This default constructor is invoked when a client wants to create a
   * RNG 
   */
  RNG() {}

  //! Return a random double in the interval (0, 1).
  double generate_random_number() const { return _ran(data); }

  //! Return the stream number.
  uint64_t get_num() const { return data[2]; }

  //! Initialize internal state from a seed and stream number.
  inline void set_seed(const uint32_t seed, const uint64_t streamnum=0);
private:
  mutable ctr_type::value_type data[4];

  //! Private copy constructor.
  RNG(const RNG &);

  //! Private assignment operator.
  RNG &operator=(const RNG &);

};

//---------------------------------------------------------------------------//
// Implementation
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//! \brief Initialize internal state from a seed and stream number.
inline void RNG::set_seed(const uint32_t seed, const uint64_t streamnum) {
  // Low bits of the counter.
  data[0] = 0;

  // High bits of the counter; used for the seed.
  data[1] = static_cast<uint64_t>(seed) << 32;

  // Low bits of the key; used for the stream number.
  data[2] = streamnum;

  // High bits of the key; used as a spawn counter.
  data[3] = 0;
}

#endif
//----------------------------------------------------------------------------//
// end of RNG.h
//----------------------------------------------------------------------------//
