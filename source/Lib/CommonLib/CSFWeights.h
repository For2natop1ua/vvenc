/* -----------------------------------------------------------------------------
The copyright in this software is being made available under the Clear BSD
License, included below. No patent rights, trademark rights and/or
other Intellectual Property Rights other than the copyrights concerning
the Software are granted under this license.

The Clear BSD License

Copyright (c) 2019-2026, Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V. & The VVenC Authors.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted (subject to the limitations in the disclaimer below) provided that
the following conditions are met:

     * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

     * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

     * Neither the name of the copyright holder nor the names of its
     contributors may be used to endorse or promote products derived from this
     software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
THIS LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


------------------------------------------------------------------------------------------- */
/** \file     CSFWeights.h
    \brief    Contrast Sensitivity Function (CSF) weighting table for perceptual quantization
*/

#pragma once

#include "TypeDef.h"

#if JVET_DEV_CSF_SCALING_LIST

#include <algorithm>
#include <cmath>

//! \ingroup CommonLib
//! \{

namespace vvenc {

namespace CSF {

// CSF weights for the 8x8 DCT block, normalized so that entry [0][0] = 16.
// Derived from the Mannos-Sakrison model: higher value = less perceptually
// visible frequency component = coarser quantization applied.
static const int csf8x8[8][8] =
{
  {  16,  16,  16,  19,  22,  26,  32,  40 },
  {  16,  16,  17,  20,  24,  30,  38,  48 },
  {  16,  17,  19,  23,  28,  35,  45,  58 },
  {  19,  20,  23,  28,  34,  43,  56,  72 },
  {  22,  24,  28,  34,  43,  55,  71,  92 },
  {  26,  30,  35,  43,  55,  71,  92, 119 },
  {  32,  38,  45,  56,  71,  92, 119, 155 },
  {  40,  48,  58,  72,  92, 119, 155, 200 }
};

// Bilinear interpolation of the 8x8 CSF table to an arbitrary TU block size.
// \param row   row index in the target block
// \param col   column index in the target block
// \param size  target block dimension (must be one of: 4, 8, 16, 32, 64)
inline int interpCSF( int row, int col, int size )
{
  const double scale = 7.0 / static_cast<double>( size - 1 );
  const double y  = row * scale;
  const double x  = col * scale;
  const int    y0 = static_cast<int>( y );
  const int    x0 = static_cast<int>( x );
  const int    y1 = std::min( y0 + 1, 7 );
  const int    x1 = std::min( x0 + 1, 7 );
  const double dy = y - y0;
  const double dx = x - x0;
  const double v  =   csf8x8[y0][x0] * ( 1.0 - dx ) * ( 1.0 - dy )
                    + csf8x8[y0][x1] *         dx    * ( 1.0 - dy )
                    + csf8x8[y1][x0] * ( 1.0 - dx ) *         dy
                    + csf8x8[y1][x1] *         dx    *         dy;
  return static_cast<int>( std::round( v ) );
}

} // namespace CSF

} // namespace vvenc

//! \}

#endif // JVET_DEV_CSF_SCALING_LIST