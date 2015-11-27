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

#include "mdd_function.h"
#include "dw_matrix.h"
#include "dw_array.h"
#include "dw_matrix_rand.h"
#include "dw_ascii.h"
#include "dw_rand.h"
#include "dw_math.h"
#include "dw_matrix_sort.h"
#include "dw_elliptical.h"
#include "dw_parse_cmd.h"
#include "dw_std.h"
#include "mdd_constant.hpp"

#include <math.h>
#include <string.h>
#include <time.h>

static PRECISION ComputeLinear(TMatrix proposal, PRECISION log_c, int *intercept);
static PRECISION ComputeInverseLinear(TMatrix posterior, PRECISION log_c, int *intercept);

/*******************************************************************************/
/************************** Bridge Method Alternative **************************/
/*******************************************************************************/
static PRECISION BridgeDifference(TMatrix proposal, TMatrix posterior, PRECISION logr)
{
  int i;
  PRECISION x, r, sum1=MINUS_INFINITY, sum2=MINUS_INFINITY, n1=RowM(posterior), n2=RowM(proposal);

  for (r=n2/n1, i=RowM(proposal)-1; i >= 0; i--)
    if ((x=ElementM(proposal,i,0)+logr-ElementM(proposal,i,1)) < 0)
      sum2=AddLogs(sum2,-log(1+r*exp(x)));
    else
      sum2=AddLogs(sum2,-x-log(exp(-x)+r));

  for (r=n1/n2, i=RowM(posterior)-1; i >= 0; i--)
    if ((x=ElementM(posterior,i,1)-ElementM(posterior,i,0)-logr) < 0)
      sum1=AddLogs(sum1,-log(1+r*exp(x)));
    else
      sum1=AddLogs(sum1,-x-log(exp(-x)+r));
    
  return sum2 - sum1;
}
#define MAX_C 1.0E50
#define TOL   1.0E-7
PRECISION ComputeLogMarginalDensity_Bridge(TMatrix proposal, TMatrix posterior)
{
  PRECISION min_c, max_c, mid_c=0.0, diff;
  int i;

  // Bracket the zero
  if ((diff=BridgeDifference(proposal,posterior,mid_c)) < 0.0)
    {
      max_c=mid_c;
      for (min_c=-1.0; min_c > -MAX_C; max_c=min_c, min_c*=10)
	if ((diff=BridgeDifference(proposal,posterior,min_c)) > 0) break;
      if (min_c <= -MAX_C) return min_c;
    }
  else
    {
      min_c=mid_c;
      for (max_c=1.0; max_c < MAX_C; min_c=max_c, max_c*=10)
	if ((diff=BridgeDifference(proposal,posterior,max_c)) < 0) break;
      if (max_c >= MAX_C) return max_c;
    }

  // Divide and conququer
  diff=BridgeDifference(proposal,posterior,mid_c=(min_c + max_c)/2.0);
  for (i=0; i < 50; i++)
    {
      if (diff > 0)
	min_c=mid_c;
      else
	max_c=mid_c;
      if ((fabs(diff=BridgeDifference(proposal,posterior,mid_c=(min_c + max_c)/2.0)) < TOL)) break;
    }
  return mid_c;
}
#undef MAX_C
#undef TOL
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/************ Mueller's method for computing marginal data density *************/
/*******************************************************************************/
#define MAX_C 1E50
static PRECISION ComputeLinear(TMatrix proposal, PRECISION log_c, int *intercept)
{
  int i;
  PRECISION slope=0.0, tmp;
  for (*intercept=0, i=RowM(proposal)-1; i >= 0; i--)
    if (log_c <= (tmp=ElementM(proposal,i,0) - ElementM(proposal,i,1)))
      {
	(*intercept)++;
	slope+=exp(log_c-tmp);
      }
  return slope;
}
static PRECISION ComputeInverseLinear(TMatrix posterior, PRECISION log_c, int *intercept)
{
  int i;
  PRECISION slope=0.0, tmp;
  for (*intercept=0, i=RowM(posterior)-1; i >= 0; i--)
    if (log_c >= (tmp=ElementM(posterior,i,0)-ElementM(posterior,i,1)))
      {
	(*intercept)++;
	slope+=exp(tmp-log_c);
      }
  return slope;
}
PRECISION ComputeLogMarginalDensity_Mueller(TMatrix proposal, TMatrix posterior,int *in1, int *in2)
{
  PRECISION log_c=0.0, slope1, slope2, N1=1.0/(PRECISION)RowM(proposal), N2=1.0/(PRECISION)RowM(posterior),
    intercept1, intercept2, max_log_c, min_log_c, diff, tmp;
  int min_in1, max_in1, min_in2, max_in2;

  slope1=ComputeLinear(proposal,log_c,in1);
  slope2=ComputeInverseLinear(posterior,log_c,in2);
  diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);

  // Bracket the intersection
  if (diff < 0.0)
    {
      do
	if (log_c < -MAX_C) 
	  return log_c;
	else
	  {
	    max_in1=*in1;
	    max_in2=*in2;
	    max_log_c=log_c;
	    log_c=10*(log_c-1);
	    slope1=ComputeLinear(proposal,log_c,in1);
	    slope2=ComputeInverseLinear(posterior,log_c,in2);
	    diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);
	  }
      while (diff < 0.0);
      min_in1=*in1;
      min_in2=*in2;
      min_log_c=log_c;
    }
  else
    {
      do
	if (log_c > MAX_C) 
	  return log_c;
	else
	  {
	    min_in1=*in1;
	    min_in2=*in2;
	    min_log_c=log_c;
	    log_c=10*(log_c+1);
	    slope1=ComputeLinear(proposal,log_c,in1);
	    slope2=ComputeInverseLinear(posterior,log_c,in2);
	    diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);
	  }
      while (diff >= 0.0);
      max_in1=*in1;
      max_in2=*in2;
      max_log_c=log_c;
    }

  // At this point diff(min_log_c) >= 0 and diff(max_log_c) < 0.
  while ((min_in1 != max_in1) || (min_in2 != max_in2))
    {
      log_c=(min_log_c + max_log_c)/2.0;
      slope1=ComputeLinear(proposal,log_c,in1);
      slope2=ComputeInverseLinear(posterior,log_c,in2);
      diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);
      if (diff > 0)
	{
	  min_in1=*in1;
	  min_in2=*in2;
	  min_log_c=log_c;
	}
      else
	{
	  max_in1=*in1;
	  max_in2=*in2;
	  max_log_c=log_c;
	}
    }

  slope1=N1*ComputeLinear(proposal,min_log_c,in1);
  intercept1=N1*(PRECISION)(*in1);
  slope2=N2*ComputeInverseLinear(posterior,min_log_c,in2);
  intercept2=N2*(PRECISION)(*in2);

  tmp=intercept1-intercept2;
  if (slope1 > 0)
    {
      tmp+=sqrt(tmp*tmp + 4*slope1*slope2);
      if (tmp >= 2.0*slope1)
	tmp=min_log_c + log(tmp) - log(2.0*slope1);
      else
	return -min_log_c;
    }
  else
    {
      if (tmp > 0)
	if (slope2 > tmp)
	  tmp=min_log_c + log(slope2) - log(tmp);
	else
	  return -min_log_c;
      else
	return -max_log_c;
    }
  return (tmp > max_log_c) ? -max_log_c : -tmp;
}
#undef MAX_C
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


TElliptical* CreateEllipticalFromPosterior_TruncatedGaussian(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi)
{
  TElliptical *elliptical=(TElliptical*)NULL;
  PRECISION a, b;
  TMatrix V;
  int i;

  if (dim > 0)
    {
      SortVectorAscending(R,R);

      i=floor(p_lo*(PRECISION)(DimV(R)-1));
      a=(i < DimV(R)-1) ? 0.5*(ElementV(R,i+1) + ElementV(R,i)) : ElementV(R,i);
      i=floor(p_hi*(PRECISION)(DimV(R)-1));
      b=(i < DimV(R)-1) ? 0.5*(ElementV(R,i+1) + ElementV(R,i)) : ElementV(R,i);

      V=ProductTransposeMM((TMatrix)NULL,scale,scale);
      elliptical=CreateElliptical_TruncatedGaussian(dim,center,V,a,b);
      FreeMatrix(V);
    }

  return elliptical;
}

TElliptical* CreateEllipticalFromPosterior_Gaussian(TVector R, int dim, TVector center, TMatrix scale)
{
  TElliptical *elliptical=(TElliptical*)NULL;
  TMatrix variance;
  if (dim > 0)
    {
      elliptical=CreateElliptical_Gaussian(dim,center,variance=ProductTransposeMM((TMatrix)NULL,scale,scale));
      FreeMatrix(variance);
    }
  return elliptical;
}

static PRECISION GetPower_truncatedpower(int dim, PRECISION lower_bound, PRECISION upper_bound)
{
  PRECISION lo=0.0001, hi=1.0, mid;
  int i;
  if (VarianceScale_TruncatedPower(dim,lower_bound,upper_bound,lo) > 1.0)
    return lo;
  for (i=0; i < 20; i++)
    if (VarianceScale_TruncatedPower(dim,lower_bound,upper_bound,hi) < 1.0)
      {
	lo=hi;
	hi*=2.0;
      }
    else
      {
	for (i=0; i < 30; i++)
	  {
	    mid=0.5*(lo+hi);
	    if (VarianceScale_TruncatedPower(dim,lower_bound,upper_bound,mid) < 1.0)
	      lo=mid;
	    else
	      hi=mid;
	  }
	return mid;
      }
  return hi;
}

TElliptical* CreateEllipticalFromPosterior_TruncatedPower(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi)
{
  TElliptical *elliptical=(TElliptical*)NULL;
  PRECISION a, b, power;
  int i;

  if (dim > 0)
    {
      SortVectorAscending(R,R);

      i=floor(p_lo*(PRECISION)(DimV(R)-1));
      a=(i < DimV(R)-1) ? 0.5*(ElementV(R,i+1) + ElementV(R,i)) : ElementV(R,i);
      i=floor(p_hi*(PRECISION)(DimV(R)-1));
      b=(i < DimV(R)-1) ? 0.5*(ElementV(R,i+1) + ElementV(R,i)) : ElementV(R,i);
      power=GetPower_truncatedpower(dim,a,b);
      
      // The code below uses the old method to choosing the power and b.  The new 
      // method chooses the power by trying to match the scale of the variance.
      /* i=floor(0.10*(PRECISION)DimV(R)); */
      /* a=(i > 0) ? 0.5*(sqrt(ElementV(R,i-1)) + sqrt(ElementV(R,i))) : sqrt(ElementV(R,i)); */
      /* i=floor(0.90*(PRECISION)DimV(R)); */
      /* b=(i > 0) ? 0.5*(sqrt(ElementV(R,i-1)) + sqrt(ElementV(R,i))) : sqrt(ElementV(R,i)); */
      /* power=log(0.10/0.90)/log(a/b); */
      /* b/=pow(0.90,1/power); */
      /* i=floor(0.10*(PRECISION)DimV(R)); */
      /* a=(i > 0) ? 0.5*(ElementV(R,i-1) + ElementV(R,i)) : ElementV(R,i); */
      /* i=floor(0.90*(PRECISION)DimV(R)); */
      /* b=(i > 0) ? 0.5*(ElementV(R,i-1) + ElementV(R,i)) : ElementV(R,i); */
      /* power=log(0.10/0.90)/log(a/b); */
      /* b/=pow(0.90,1/power); */

      elliptical=CreateElliptical_TruncatedPower(dim,center,scale,a,b,power);
    }

  return elliptical;
}

TElliptical* CreateEllipticalFromPosterior_Power(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_hi)
{
  TElliptical *elliptical=(TElliptical*)NULL;
  PRECISION b, power;
  int i;

  if (dim > 0)
    {
      SortVectorAscending(R,R);

      i=floor(p_hi*(PRECISION)(DimV(R)-1));
      b=(i < DimV(R)-1) ? 0.5*(ElementV(R,i+1) + ElementV(R,i)) : ElementV(R,i);
      power=GetPower_truncatedpower(dim,0.0,b);

      // The code below uses the old method to choosing the power and b.  The new 
      // method chooses the power by trying to match the scale of the variance.
      /* i=floor(0.20*(PRECISION)DimV(R)); */
      /* PRECISION a = (i > 0) ? 0.5*(sqrt(ElementV(R,i-1)) + sqrt(ElementV(R,i))) : sqrt(ElementV(R,i)); */
      /* i=floor(0.80*(PRECISION)DimV(R)); */
      /* b=(i > 0) ? 0.5*(sqrt(ElementV(R,i-1)) + sqrt(ElementV(R,i))) : sqrt(ElementV(R,i)); */
      /* power=log(0.20/0.80)/log(a/b); */
      /* b/=pow(0.80,1/power); */
      /* i=floor(0.20*(PRECISION)DimV(R)); */
      /* PRECISION a = (i > 0) ? 0.5*(ElementV(R,i-1) + ElementV(R,i)) : ElementV(R,i); */
      /* i=floor(0.80*(PRECISION)DimV(R)); */
      /* b=(i > 0) ? 0.5*(ElementV(R,i-1) + ElementV(R,i)) : ElementV(R,i); */
      /* power=log(0.20/0.80)/log(a/b); */
      /* b/=pow(0.80,1/power); */

      elliptical=CreateElliptical_Power(dim,center,scale,b,power);
    }

  return elliptical;
}

TElliptical* CreateEllipticalFromPosterior_Uniform(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi)
{
  TElliptical *elliptical=(TElliptical*)NULL;
  PRECISION a, b;
  int i;

  if (dim > 0)
    {
      SortVectorAscending(R,R);

      i=floor(p_lo*(PRECISION)(DimV(R)-1));
      a=(i < DimV(R)-1) ? 0.5*(ElementV(R,i+1) + ElementV(R,i)) : ElementV(R,i);
      i=floor(p_hi*(PRECISION)(DimV(R)-1));
      b=(i < DimV(R)-1) ? 0.5*(ElementV(R,i+1) + ElementV(R,i)) : ElementV(R,i);

      elliptical=CreateElliptical_Uniform(dim,center,scale,a,b);
    }

  return elliptical;
}

TElliptical* CreateEllipticalFromPosterior_Step(TVector R, int dim, TVector center, TMatrix scale, PRECISION p_lo, PRECISION p_hi)
{
  TElliptical *elliptical=(TElliptical*)NULL;
  PRECISION *table, inc, x, y;
  int i, j, k, m=30;

  if (dim > 0)
    {
      table=(PRECISION*)dw_malloc((m+1)*sizeof(PRECISION));

      SortVectorAscending(R,R);

      for (k=m; k >= 1; k--)
	{
	  inc=(p_hi-p_lo)/(PRECISION)k;
	  for (x=p_lo, j=0; j <= k; x+=inc, j++)
	    {
	      i=floor(y=x*(PRECISION)(DimV(R)-1));
	      if (i >= DimV(R)-1)
		table[j]=ElementV(R,DimV(R)-1);
	      else
		table[j]=(1.0-y+i)*ElementV(R,i) + (y-i)*ElementV(R,i+1);
	      //table[j]=(i > 0) ? 0.5*(ElementV(R,i-1) + ElementV(R,i)) : ElementV(R,i);
	    }
	  for (j=1; j <= k; j++)
	    if (table[j] <= table[j-1]) break;
	  if (j > k) break;
	}

      if (k <= 0)
	{
	  dw_free(table);
	  printf("Not enough variation in the posterior draws to form proposal\n");
	  dw_exit(0);
	}

      elliptical=CreateElliptical_Step(dim,center,scale,table,k);
      dw_free(table);
    }

  return elliptical;
}
/*******************************************************************************/
/*******************************************************************************/
