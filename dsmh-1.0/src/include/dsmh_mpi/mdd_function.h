/*
 * Copyright (C) 1996-2011 Daniel Waggoner and Tao Zha
 *
 * This free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * If you did not received a copy of the GNU General Public License
 * with this software, see <http://www.gnu.org/licenses/>.
 */

#ifndef __DW_MDD_SWITCH__
#define __DW_MDD_SWITCH__

#include "dw_elliptical.h"

#ifdef __cplusplus
extern "C"
{
#endif


PRECISION ComputeLogMarginalDensity_Mueller(TMatrix proposal, TMatrix posterior,int *in1, int *in2);
PRECISION ComputeLogMarginalDensity_Bridge(TMatrix proposal, TMatrix posterior);

TElliptical* CreateEllipticalFromPosterior_TruncatedGaussian(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi);
TElliptical* CreateEllipticalFromPosterior_Gaussian(TVector R, int dim, TVector center, TMatrix scale);
TElliptical* CreateEllipticalFromPosterior_TruncatedPower(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi);
TElliptical* CreateEllipticalFromPosterior_Power(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_hi);
TElliptical* CreateEllipticalFromPosterior_Uniform(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi);
TElliptical* CreateEllipticalFromPosterior_Step(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi);

  
#ifdef __cplusplus
}
#endif  
  
#endif
