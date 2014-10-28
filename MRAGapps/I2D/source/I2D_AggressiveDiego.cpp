
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "I2D_AggressiveDiego.h"

#include <complex>

#ifndef max
#define max(a,b) ((a)>=(b)? (a) : (b))
#endif

namespace AggressiveDiego
{
	inline __m128 worse_division(const __m128 a, const __m128 b)
	{
		//from the "Software Optimization Guide for AMD Family 10h and 12h Processors"
		__m128 x0, x1, x2, x3, y;
		x0 = _mm_rcp_ps (b); 
		x1 = _mm_mul_ps (a, x0); 
		x2 = _mm_mul_ps (b, x0); 
		x3 = _mm_sub_ps (_mm_set1_ps (2.0F), x2);
		
		return _mm_mul_ps (x1, x3);
	}
	
#ifdef __INTEL_COMPILER
	inline __m128 operator+(__m128 a, __m128 b){ return _mm_add_ps(a, b); }
	inline __m128 operator*(__m128 a, __m128 b){ return _mm_mul_ps(a, b); }
	inline __m128 operator-(__m128 a,  __m128 b){ return _mm_sub_ps(a, b); }
	
	inline __m128d operator+(__m128d a, __m128d b){ return _mm_add_pd(a, b); }
	inline __m128d operator*(__m128d a, __m128d b){ return _mm_mul_pd(a, b); }
	inline __m128d operator-(__m128d a,  __m128d b){ return _mm_sub_pd(a, b); }
#endif
	
	template <int c>
	inline __m128 _select(__m128 d)  {abort(); return d;}
	
	template <>
	inline __m128 _select<0>(__m128 d)  { return _mm_shuffle_ps(d,d,  _MM_SHUFFLE(0,0,0,0));}
	template <>
	inline __m128 _select<1>(__m128 d)  { return _mm_shuffle_ps(d,d,  _MM_SHUFFLE(1,1,1,1));}
	template <>
	inline __m128 _select<2>(__m128 d)  { return _mm_shuffle_ps(d,d,  _MM_SHUFFLE(2,2,2,2));}
	template <>
	inline __m128 _select<3>(__m128 d)  { return _mm_shuffle_ps(d,d,  _MM_SHUFFLE(3,3,3,3));}
	
#define _MYSMALL_EPS_ 0
	
	void view(__m128 v, const char * n)
	{
		float data[4];
		
		_mm_storeu_ps(data, v);
		printf("%s: %.2f %.2f %.2f %.2f\n", n, data[0], data[1], data[2], data[3]);
	}
	
	/*********************************************************************************************************************************/
	/********************************                   DIRECT INTERACTIONS                     **************************************/
	/*********************************************************************************************************************************/
	
	struct DirectInteractions
	{
		template<bool safe>
		inline void kernel(const __m128 xd, const __m128 yd, const __m128 xs, const __m128 ys, const __m128 ws, float * const u_dest, float * const v_dest);
		
		inline __m128 get(const float * const ptr, const int s, const float val=-1e7)
		{
			switch(s)
			{
				case 1:
					return _mm_set_ps(val, val, val, ptr[0]);
				case 2:
					return _mm_set_ps(val, val, ptr[1], ptr[0]);
				case 3:
					return _mm_set_ps(val, ptr[2], ptr[1], ptr[0]);
				default:
				{
					abort();
					return _mm_set_ps(val, val, val, val);
				}
			}
		}
		
		inline __m128 get(const float * const ptr)
		{
			return _mm_load_ps(ptr);
		}
		
		template<bool safe>
		void evaluate_tile(const __m128 x0, const __m128 y0, const __m128 h, const __m128 d,
						   const __m128 xs, const __m128 ys, const __m128 ws, float * const u, float * const v)
		{
			for(int dy=0; dy<_BLOCKSIZE_; dy++)
			{
				const __m128 y = _mm_set_ps1(dy)*h + y0;
				const int ybase = dy*_BLOCKSIZE_;
				
				for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
				{
					const __m128 x = (d+_mm_set_ps1(dx))*h + x0;
					
					kernel<safe>(x, y, xs, ys, ws, u+ybase+dx, v+ybase+dx);
				}
			}
		}
		
		template<bool safe>
		void evaluate_tile(const float * const xd, const float * const yd,
						   const __m128 xs, const __m128 ys, const __m128 ws, float * const u, float * const v)
		{
			for(int dy=0; dy<_BLOCKSIZE_; dy++)
			{
				const int ybase = dy*_BLOCKSIZE_;
				const __m128 y =  _mm_load_ps(yd+dy*4);
#pragma ivdep
				for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
				{
					const int offset1 = ybase+dx;
					kernel<safe>(_mm_load_ps(xd+dx), y, xs, ys, ws, u+offset1, v+offset1);
				}
			}
		}
		
		template<bool safe>
		void evaluate(const float x_start, const float y_start, const float h,
					  const float * const xs, const float * const ys, const float * const ws, const int nsources,
					  float * const u, float * const v)
		{
			const int offset = (((unsigned long int)xs) & 0xf)/4;
			const int start_tile = 4*(int)ceil(offset/4.) - offset;
			const int end_tile = 4*(int)ceil((offset + nsources - 3)/4.) - offset;
			
			assert((((unsigned long int)xs) & 0x3) == 0);
			assert((((unsigned long int)&xs[start_tile]) & 0xf) == 0);
			assert((((unsigned long int)&xs[end_tile]) & 0xf) == 0);
			assert(0 <= start_tile);
			assert(start_tile <= end_tile);
			assert(end_tile <= nsources);
			
			const int nrem1 = start_tile;
			const int nrem2 = nsources - end_tile;
			
			const __m128 x0 = _mm_set_ps1(x_start);
			const __m128 y0 = _mm_set_ps1(y_start);
			const __m128 dx = _mm_set_ps1(h);
			const __m128 d = _mm_set_ps(+3,+2,+1,+0);
			
			//evaluate first unaligned particles
			if (nrem1 > 0)
			{
				const __m128 first_xs = get(xs, nrem1);
				const __m128 first_ys = get(ys, nrem1);
				const __m128 first_ws = get(ws, nrem1, 0);
				
				evaluate_tile<safe>(x0, y0, dx, d, first_xs, first_ys, first_ws, u, v);
			}
			
			//evaluate aligned particles
			for(int s=start_tile; s<end_tile; s+=4)
				evaluate_tile<safe>(x0, y0, dx, d, get(xs+s), get(ys+s), get(ws+s), u, v);
			
			//evaluate last unaligned particles
			if (nrem2 > 0)
			{
				const __m128 last_xs = get(xs + end_tile, nrem2);
				const __m128 last_ys = get(ys + end_tile, nrem2);
				const __m128 last_ws = get(ws + end_tile, nrem2, 0);
				
				evaluate_tile<safe>(x0, y0, dx, d, last_xs, last_ys, last_ws, u, v);
			}
		}
		
		template<bool safe>
		void evaluate(const float * const xd, const float * const yd, 
					  const float * const xs, const float * const ys, const float * const ws, const int nsources,
					  float * const u, float * const v)
		{
			const int offset = (((unsigned long int)xs) & 0xf)/4;
			const int start_tile = 4*(int)ceil(offset/4.) - offset;
			const int end_tile = 4*(int)ceil((offset + nsources - 3)/4.) - offset;
			
			assert((((unsigned long int)xs) & 0x3) == 0);
			assert((((unsigned long int)&xs[start_tile]) & 0xf) == 0);
			assert((((unsigned long int)&xs[end_tile]) & 0xf) == 0);
			assert(0 <= start_tile);
			assert(start_tile <= end_tile);
			assert(end_tile <= nsources);
			
			const int nrem1 = start_tile;
			const int nrem2 = nsources - end_tile;
			
			//evaluate first unaligned particles
			if (nrem1 > 0)
			{
				const __m128 first_xs = get(xs, nrem1);
				const __m128 first_ys = get(ys, nrem1);
				const __m128 first_ws = get(ws, nrem1, 0);
				
				evaluate_tile<safe>(xd, yd, first_xs, first_ys, first_ws, u, v);
			}
			
			//evaluate aligned particles
			for(int s=start_tile; s<end_tile; s+=4)
				evaluate_tile<safe>(xd, yd, get(xs+s), get(ys+s), get(ws+s), u, v);
			
			//evaluate last unaligned particles
			if (nrem2 > 0)
			{
				const __m128 last_xs = get(xs + end_tile, nrem2);
				const __m128 last_ys = get(ys + end_tile, nrem2);
				const __m128 last_ws = get(ws + end_tile, nrem2, 0);
				
				evaluate_tile<safe>(xd, yd, last_xs, last_ys, last_ws, u, v);
			}
		}
	};
	
	template<>
	void DirectInteractions::kernel<false>(const __m128 xd, const __m128 yd, const __m128 xs, const __m128 ys, const __m128 ws, float * const u_dest, float * const v_dest)
	{
		const __m128 rx1 = xd - _select<0>(xs);
		const __m128 ry1 = yd - _select<0>(ys);
		const __m128 r2A = rx1*rx1 + ry1*ry1;
		
		const __m128 rx2 = xd - _select<1>(xs);
		const __m128 ry2 = yd - _select<1>(ys);
		const __m128 r2B = rx2*rx2 + ry2*ry2;
		
		const __m128 rx3 = xd - _select<2>(xs);
		const __m128 ry3 = yd - _select<2>(ys);
		const __m128 r2C = rx3*rx3 + ry3*ry3;
		
		const __m128 rx4 = xd - _select<3>(xs);
		const __m128 ry4 = yd - _select<3>(ys);
		const __m128 r2D = rx4*rx4 + ry4*ry4;
		
#ifdef _FMM_NOPRECDIV_KERNELS_
		const __m128 factor1 =  worse_division(_select<0>(ws), r2A);
		const __m128 factor2 =  worse_division(_select<1>(ws), r2B);
		const __m128 factor3 =  worse_division(_select<2>(ws), r2C);
		const __m128 factor4 =  worse_division(_select<3>(ws), r2D);
#else
		const __m128 factor1 =  _mm_div_ps(_select<0>(ws), r2A);
		const __m128 factor2 =  _mm_div_ps(_select<1>(ws), r2B);
		const __m128 factor3 =  _mm_div_ps(_select<2>(ws), r2C);
		const __m128 factor4 =  _mm_div_ps(_select<3>(ws), r2D);
#endif
		
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor1*ry1));
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor2*ry2));
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor3*ry3));
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor4*ry4));
		
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor1*rx1));
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor2*rx2));
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor3*rx3));
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor4*rx4));
	}
	
	template<>
	void DirectInteractions::kernel<true>(const __m128 xd, const __m128 yd, const __m128 xs, const __m128 ys, const __m128 ws, float * const u_dest, float * const v_dest)
	{
		const __m128 rx1 = xd - _select<0>(xs);
		const __m128 ry1 = yd - _select<0>(ys);
		const __m128 r2A = rx1*rx1 + ry1*ry1;
		const __m128 flagA = _mm_cmpnle_ps(r2A, _mm_set_ps1(0));
		
		const __m128 rx2 = xd - _select<1>(xs);
		const __m128 ry2 = yd - _select<1>(ys);
		const __m128 r2B = rx2*rx2 + ry2*ry2;
		const __m128 flagB = _mm_cmpnle_ps(r2B, _mm_set_ps1(0));
		
		const __m128 rx3 = xd - _select<2>(xs);
		const __m128 ry3 = yd - _select<2>(ys);
		const __m128 r2C = rx3*rx3 + ry3*ry3;
		const __m128 flagC = _mm_cmpnle_ps(r2C, _mm_set_ps1(0));
		
		const __m128 rx4 = xd - _select<3>(xs);
		const __m128 ry4 = yd - _select<3>(ys);
		const __m128 r2D = rx4*rx4 + ry4*ry4;
		const __m128 flagD = _mm_cmpnle_ps(r2D, _mm_set_ps1(0));
		
#ifdef _FMM_NOPRECDIV_KERNELS_
		const __m128 factor1 = _mm_and_ps(flagA, worse_division(_select<0>(ws), r2A));
		const __m128 factor2 = _mm_and_ps(flagB, worse_division(_select<1>(ws), r2B));
		const __m128 factor3 = _mm_and_ps(flagC, worse_division(_select<2>(ws), r2C));
		const __m128 factor4 = _mm_and_ps(flagD, worse_division(_select<3>(ws), r2D));
#else
		const __m128 factor1 = _mm_and_ps(flagA, _mm_div_ps(_select<0>(ws), r2A));
		const __m128 factor2 = _mm_and_ps(flagB, _mm_div_ps(_select<1>(ws), r2B));
		const __m128 factor3 = _mm_and_ps(flagC, _mm_div_ps(_select<2>(ws), r2C));
		const __m128 factor4 = _mm_and_ps(flagD, _mm_div_ps(_select<3>(ws), r2D));
#endif
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor1*ry1));
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor2*ry2));
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor3*ry3));
		_mm_store_ps(u_dest, _mm_sub_ps(_mm_load_ps(u_dest), factor4*ry4));
		
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor1*rx1));
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor2*rx2));
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor3*rx3));
		_mm_store_ps(v_dest, _mm_add_ps(_mm_load_ps(v_dest), factor4*rx4));
	}
	
	void direct_interactions_sse(float xd[_BLOCKSIZE_], float (* const yd)[4],
								 const float * const xs, const float * const ys, const float * const ws, const int nsources,
								 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const bool safe)
	{
#if (_BLOCKSIZE_ < 8)
#pragma warning _BLOCKSIZE_ is not big enough.
		printf("aborting, _BLOCKSIZE_ is not greater or equal than 8.");
		abort();
#endif
		
		assert (nsources > 0);
		
		if (( ((unsigned long int)xd) & 0xf) != 0 )
		{
			printf("const float * const xd is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)yd) & 0xf) != 0 )
		{
			printf("const float * const * yd is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)u) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)v) & 0xf) != 0 )
		{
			printf("float * const v is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)xs) & 0xf) != (((unsigned long int)ys) & 0xf ))
		{
			printf("const float * const xs not aligned as ys!\n");
			abort();
		}
		
		if (( ((unsigned long int)xs) & 0xf) != (((unsigned long int)ws) & 0xf ))
		{
			printf("const float * const xs not aligned as ws!\n");
			abort();
		}
		/*
		 #ifndef NDEBUG
		 for(int iy=0; iy<_BLOCKSIZE_; iy++)
		 {
		 for(int ix=0; ix<_BLOCKSIZE_; ix++)
		 {
		 printf("%f ", yd[iy][ix]);
		 }
		 printf("\n");
		 }
		 
		 //exit(0);
		 #endif
		 */
		
		DirectInteractions direct_interactions;
		
		if (safe)
			direct_interactions.evaluate<true>((float*)xd, (float*)yd, xs, ys, ws, nsources, (float *)u, (float *)v);
		else
			direct_interactions.evaluate<false>((float*)xd, (float*)yd, xs, ys, ws, nsources, (float *)u, (float *)v);
		
	}
	
	void direct_interactions_sse(const float x_start, const float y_start, const float h,
								 const float * const xs, const float * const ys, const float * const ws, const int nsources,
								 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const bool safe)
	{
#if (_BLOCKSIZE_ < 8)
#pragma warning _BLOCKSIZE_ is not big enough.
		printf("aborting, _BLOCKSIZE_ is not greater or equal than 8.");
		abort();
#endif
		
		assert (nsources > 0);
		
		if (( ((unsigned long int)u) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)v) & 0xf) != 0 )
		{
			printf("float * const v is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)xs) & 0xf) != (((unsigned long int)ys) & 0xf ))
		{
			printf("const float * const xs not aligned as ys!\n");
			abort();
		}
		
		if (( ((unsigned long int)xs) & 0xf) != (((unsigned long int)ws) & 0xf ))
		{
			printf("const float * const xs not aligned as ws!\n");
			abort();
		}
		
		DirectInteractions direct_interactions;
		
		if (safe)
			direct_interactions.evaluate<true>(x_start, y_start, h, xs, ys, ws, nsources, (float *)u, (float *)v);
		//direct_interactions_cpp(x_start, y_start, h, xs, ys, ws, nsources, u, v, safe);
		else
			direct_interactions.evaluate<false>(x_start, y_start, h, xs, ys, ws, nsources, (float *)u, (float *)v);
		
	}
	
	void direct_interactions_cpp(const float x_start, const float y_start, const float h,
								 const float * const xs, const float * const ys, const float * const ws, const int nsources,
								 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const bool safe)
	{
		assert (nsources > 0);
		
		typedef float Real;
		for(int dy=0; dy<_BLOCKSIZE_; dy++)
			for(int dx=0; dx<_BLOCKSIZE_; dx++)
			{
				Real localu[2] = {0,0};
				
				for(int s=0; s<nsources; s++)
				{
					Real r[2] = {
						x_start + dx*h - xs[s],
						y_start + dy*h - ys[s]
					};
					
					const Real distance_2 = r[0]*r[0] + r[1]*r[1];
					
					if (distance_2==0) continue;
					
					localu[0] -= (ws[s]*r[1])/distance_2;
					localu[1] += (ws[s]*r[0])/distance_2;
				}
				
				u[dy][dx] += localu[0];
				v[dy][dx] += localu[1];
			}
	}
	
	void check_quality(const float x_start, const float y_start, const float h,
					   const float * const xs, const float * const ys, const float * const ws, const int nsources,
					   float (* const ures)[_BLOCKSIZE_], float (* const vres)[_BLOCKSIZE_], const double tol)
	{
		assert (nsources > 0);
		
		double L1[2]= {0, 0};
		double LINF[2] = {0, 0};
		
		for(int dy=0; dy<_BLOCKSIZE_; dy++)
			for(int dx=0; dx<_BLOCKSIZE_; dx++)
			{
				double u=0, v=0;
				
				for(int s=0; s<nsources; s++)
				{
					const double rx = x_start + dx*h - xs[s];
					const double ry = y_start + dy*h - ys[s];
					const double r2 = rx*rx + ry*ry;
					const double safety_value = (r2 >= _MYSMALL_EPS_);
					const double factor = safety_value / (r2*safety_value + (1-safety_value));
					
					u -= ry*ws[s]*factor;
					v += rx*ws[s]*factor;
				}
				
				const double uactual = ures[dy][dx];
				const double vactual = vres[dy][dx];
				
				const double diffu = u-uactual;
				const double diffv = v-vactual;
				
				assert(fabs(diffu) < tol);
				assert(fabs(diffv) < tol);
				
				L1[0] += fabs(diffu);
				L1[1] += fabs(diffv);
				
				LINF[0] = max(LINF[0], fabs(diffu));
				LINF[1] = max(LINF[1], fabs(diffv));
			}
		
		printf("\n====================================================================\n");
		printf("============     QUALITY CHECK (%.1e)   PASSED      =============\n", tol);
		printf("====================================================================\n");
		
		printf("linf-norm discrepancies: u: %e v: %e\n", LINF[0], LINF[1]);
		printf("l1-norm discrepancies: u: %e v: %e\n", L1[0], L1[1]);
		printf("average discrepancies: u: %e v: %e\n", L1[0]/(double)(_BLOCKSIZE_*_BLOCKSIZE_), L1[1]/(double)(_BLOCKSIZE_*_BLOCKSIZE_));
		printf("average interaction discrepancy: u: %e v: %e\n", L1[0]/(double)(nsources*_BLOCKSIZE_*_BLOCKSIZE_), L1[1]/(double)(nsources*_BLOCKSIZE_*_BLOCKSIZE_));
	}
	
	
	/*********************************************************************************************************************************/
	/********************************                   INDIRECT INTERACTIONS                   **************************************/
	/*********************************************************************************************************************************/
	
	struct IndirectInteractions
	{
		template<int s>
		inline __m128d sum_step(const __m128 rexp, const __m128 iexp, const __m128d prod)
		{
			return _mm_addsub_pd(_mm_cvtps_pd(_mm_shuffle_ps(rexp, rexp, _MM_SHUFFLE(s,s,s,s)))*prod,
								 _mm_cvtps_pd(_mm_shuffle_ps(iexp, iexp, _MM_SHUFFLE(s,s,s,s)))*_mm_shuffle_pd(prod, prod, _MM_SHUFFLE2(0,1)));
		}
		
		template<int d>
		inline __m128d prod_step(const __m128 realrp, const __m128 imagrp, const __m128d oldprod)
		{
			return _mm_addsub_pd(_mm_cvtps_pd(_mm_shuffle_ps(realrp, realrp, _MM_SHUFFLE(d,d,d,d)))*oldprod,
								 _mm_cvtps_pd(_mm_setzero_ps() - _mm_shuffle_ps(imagrp, imagrp, _MM_SHUFFLE(d,d,d,d)))*_mm_shuffle_pd(oldprod, oldprod, _MM_SHUFFLE2(0,1)));
		}
		
		inline __m128 rearrange_real(const __m128d s0, const __m128d s1, const __m128d s2, const __m128d s3)
		{
			return _mm_shuffle_ps(_mm_cvtpd_ps( _mm_shuffle_pd(s0, s1, _MM_SHUFFLE2(0,0))),
								  _mm_cvtpd_ps( _mm_shuffle_pd(s2, s3, _MM_SHUFFLE2(0,0))),
								  _MM_SHUFFLE(1,0,1,0));
		}
		
		inline __m128 rearrange_imag(const __m128d s0, const __m128d s1, const __m128d s2, const __m128d s3)
		{
			return _mm_shuffle_ps(_mm_cvtpd_ps( _mm_shuffle_pd(s0, s1, _MM_SHUFFLE2(1,1))),
								  _mm_cvtpd_ps( _mm_shuffle_pd(s2, s3, _MM_SHUFFLE2(1,1))),
								  _MM_SHUFFLE(1,0,1,0));
		}
		
		void compute(const float xstart, const float ystart, const float h_,
					 const float xcenter, const float ycenter,
					 const float * const rexpansions, const float * const iexpansions,
					 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_])
		{
#if (_ORDER_ % 4 != 0)
#pragma warning _ORDER_ is not multiple of 4
			printf("aborting, _ORDER_ is not a multiple of 4.");
			abort();
#endif
			
			const __m128 x0 = _mm_set_ps1(xstart-xcenter);
			const __m128 y0 = _mm_set_ps1(ystart-ycenter);
			const __m128 h = _mm_set_ps1(h_);
			const __m128 d = _mm_set_ps(+3,+2,+1,+0);
			
			for(int dy=0; dy<_BLOCKSIZE_; dy++)
			{
				const __m128 imagrp_ = y0 + _mm_set_ps1(dy)*h;
				
				for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
				{
					const __m128 realrp_ = x0 + (d+_mm_set_ps1(dx))*h;
#ifdef _FMM_NOPRECDIV_KERNELS_					
					const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
					const __m128 realrp = worse_division(realrp_, abs_rp);
					const __m128 imagrp = worse_division(imagrp_, abs_rp);
#else
					const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
					const __m128 realrp = _mm_div_ps(realrp_, abs_rp);
					const __m128 imagrp = _mm_div_ps(imagrp_, abs_rp);
#endif					
					__m128d sum0 = _mm_setzero_pd();
					__m128d sum1 = _mm_setzero_pd();
					__m128d sum2 = _mm_setzero_pd();
					__m128d sum3 = _mm_setzero_pd();
					
					__m128d pro0 = _mm_set_pd(0,1);
					__m128d pro1 = _mm_set_pd(0,1);
					__m128d pro2 = _mm_set_pd(0,1);
					__m128d pro3 = _mm_set_pd(0,1);
					
					for(int e=0; e<_ORDER_; e+=4)
					{
						const __m128 rexp = _mm_load_ps(rexpansions+e);
						const __m128 iexp = _mm_load_ps(iexpansions+e);
						
						sum0 = sum0 + sum_step<0>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<0>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<0>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<0>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
						
						sum0 = sum0 + sum_step<1>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<1>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<1>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<1>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
						
						sum0 = sum0 + sum_step<2>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<2>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<2>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<2>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
						
						sum0 = sum0 + sum_step<3>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<3>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<3>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<3>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
					}
					
					const __m128 realsum = rearrange_real(sum0, sum1, sum2, sum3);
					const __m128 imagsum = rearrange_imag(sum0, sum1, sum2, sum3);
					
					_mm_store_ps(&u[dy][dx], _mm_add_ps(_mm_load_ps(&u[dy][dx]), imagsum*realrp - realsum*imagrp));
					_mm_store_ps(&v[dy][dx], _mm_add_ps(_mm_load_ps(&v[dy][dx]), imagsum*imagrp + realsum*realrp));
				}
			}
		}
		
		void compute(const float * const xd, const float * const yd,
					 const float xcenter_, const float ycenter_,
					 const float * const rexpansions, const float * const iexpansions,
					 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_])
		{
#if (_ORDER_ % 4 != 0)
#pragma warning _ORDER_ is not multiple of 4
			printf("aborting, _ORDER_ is not a multiple of 4.");
			abort();
#endif
			
			const __m128 xcenter = _mm_set_ps1(xcenter_);
			const __m128 ycenter = _mm_set_ps1(ycenter_);
			
			for(int dy=0; dy<_BLOCKSIZE_; dy++)
			{
				const int ybase = _BLOCKSIZE_*dy;
				const __m128 imagrp_ = _mm_load_ps(yd+dy*4) - ycenter;
				
				for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
				{
					const __m128 realrp_ = _mm_load_ps(xd+dx) - xcenter;
					
#ifdef _FMM_NOPRECDIV_KERNELS_
					const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
					const __m128 realrp = worse_division(realrp_, abs_rp);
					const __m128 imagrp = worse_division(imagrp_, abs_rp);
#else					
					const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
					const __m128 realrp = _mm_div_ps(realrp_, abs_rp);
					const __m128 imagrp = _mm_div_ps(imagrp_, abs_rp);
#endif					
					__m128d sum0 = _mm_setzero_pd();
					__m128d sum1 = _mm_setzero_pd();
					__m128d sum2 = _mm_setzero_pd();
					__m128d sum3 = _mm_setzero_pd();
					
					__m128d pro0 = _mm_set_pd(0,1);
					__m128d pro1 = _mm_set_pd(0,1);
					__m128d pro2 = _mm_set_pd(0,1);
					__m128d pro3 = _mm_set_pd(0,1);
					
					for(int e=0; e<_ORDER_; e+=4)
					{
						const __m128 rexp = _mm_load_ps(rexpansions+e);
						const __m128 iexp = _mm_load_ps(iexpansions+e);
						
						sum0 = sum0 + sum_step<0>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<0>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<0>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<0>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
						
						sum0 = sum0 + sum_step<1>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<1>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<1>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<1>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
						
						sum0 = sum0 + sum_step<2>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<2>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<2>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<2>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
						
						sum0 = sum0 + sum_step<3>(rexp, iexp, pro0);
						sum1 = sum1 + sum_step<3>(rexp, iexp, pro1);
						sum2 = sum2 + sum_step<3>(rexp, iexp, pro2);
						sum3 = sum3 + sum_step<3>(rexp, iexp, pro3);
						
						pro0 = prod_step<0>(realrp, imagrp, pro0);
						pro1 = prod_step<1>(realrp, imagrp, pro1);
						pro2 = prod_step<2>(realrp, imagrp, pro2);
						pro3 = prod_step<3>(realrp, imagrp, pro3);
					}
					
					const __m128 realsum = rearrange_real(sum0, sum1, sum2, sum3);
					const __m128 imagsum = rearrange_imag(sum0, sum1, sum2, sum3);
					
					_mm_store_ps(&u[dy][dx], _mm_add_ps(_mm_load_ps(&u[dy][dx]), imagsum*realrp - realsum*imagrp));
					_mm_store_ps(&v[dy][dx], _mm_add_ps(_mm_load_ps(&v[dy][dx]), imagsum*imagrp + realsum*realrp));
				}
			}
		}
	};
	
	void indirect_interactions_sse(const float xstart, const float ystart, const float h_,
								   const float xcenter, const float ycenter,
								   const float * const rexpansions, const float * const iexpansions,
								   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_])
	{
		if (( ((unsigned long int)u) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)v) & 0xf) != 0 )
		{
			printf("float * const v is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)rexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)iexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		
		IndirectInteractions indirect;
		
		indirect.compute(xstart, ystart, h_, xcenter, ycenter, rexpansions, iexpansions, u,v);
		
	}
	
	void indirect_interactions_sse(float xd[_BLOCKSIZE_], float (* const yd)[4],
								   const float xcenter, const float ycenter,
								   const float * const rexpansions, const float * const iexpansions,
								   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_])
	{
		if (( ((unsigned long int)xd) & 0xf) != 0 )
		{
			printf("float * const xd is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)yd) & 0xf) != 0 )
		{
			printf("float * const yd is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)u) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)v) & 0xf) != 0 )
		{
			printf("float * const v is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)rexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)iexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		
		IndirectInteractions indirect;
		
		indirect.compute((float*)xd, (float*)yd, xcenter, ycenter, rexpansions, iexpansions, u,v);
		
	}
	
	void indirect_interactions_ssefloat(const float xstart, const float ystart, const float h_,
										const float xcenter, const float ycenter,
										const float * const rexpansions, const float * const iexpansions,
										float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_])
	{
#if (_ORDER_ % 4 != 0)
#pragma warning _ORDER_ is not multiple of 4
		printf("aborting, _ORDER_ is not a multiple of 4.");
		abort();
#endif
		if (( ((unsigned long int)u) & 0xf ) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)v) & 0xf) != 0 )
		{
			printf("float * const v is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)rexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)iexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		const __m128 x0 = _mm_set_ps1(xstart-xcenter);
		const __m128 y0 = _mm_set_ps1(ystart-ycenter);
		const __m128 h = _mm_set_ps1(h_);
		const __m128 d = _mm_set_ps(+3,+2,+1,+0);
		
		for(int dy=0; dy<_BLOCKSIZE_; dy++)
		{
			const __m128 imagrp_ = y0 + _mm_set_ps1(dy)*h;
			
			for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
			{
				const __m128 realrp_ = x0 + (d+_mm_set_ps1(dx))*h;
#ifdef _FMM_NOPRECDIV_KERNELS_
				const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
				const __m128 realrp = worse_division(realrp_, abs_rp);
				const __m128 imagrp = worse_division(imagrp_, abs_rp);
#else
				const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
				const __m128 realrp = _mm_div_ps(realrp_, abs_rp);
				const __m128 imagrp = _mm_div_ps(imagrp_, abs_rp);
#endif	
				
				__m128 realsum = _mm_setzero_ps();
				__m128 imagsum = _mm_setzero_ps();
				
				__m128 realprod = _mm_set_ps(1,1,1,1);
				__m128 imagprod = _mm_setzero_ps();
				
				for(int e=0; e<_ORDER_; e+=4)
				{
					const __m128 rexps = _mm_load_ps(rexpansions+e);
					const __m128 iexps = _mm_load_ps(iexpansions+e);
					
					realsum = realsum + realprod*_select<0>(rexps) - imagprod*_select<0>(iexps);
					imagsum = imagsum + imagprod*_select<0>(rexps) + realprod*_select<0>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
					realsum = realsum + realprod*_select<1>(rexps) - imagprod*_select<1>(iexps);
					imagsum = imagsum + imagprod*_select<1>(rexps) + realprod*_select<1>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
					realsum = realsum + realprod*_select<2>(rexps) - imagprod*_select<2>(iexps);
					imagsum = imagsum + imagprod*_select<2>(rexps) + realprod*_select<2>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
					realsum = realsum + realprod*_select<3>(rexps) - imagprod*_select<3>(iexps);
					imagsum = imagsum + imagprod*_select<3>(rexps) + realprod*_select<3>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
				}
				
				_mm_store_ps(&u[dy][dx], _mm_add_ps(_mm_load_ps(&u[dy][dx]), imagsum*realrp - realsum*imagrp));
				_mm_store_ps(&v[dy][dx], _mm_add_ps(_mm_load_ps(&v[dy][dx]), imagsum*imagrp + realsum*realrp));
			}
		}
	}
	
	
	void indirect_interactions_ssefloat(float xd[_BLOCKSIZE_], float (* const yd_)[4],
										const float xcenter_, const float ycenter_,
										const float * const rexpansions, const float * const iexpansions,
										float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_])
	{
#if (_ORDER_ % 4 != 0)
#pragma warning _ORDER_ is not multiple of 4
		printf("aborting, _ORDER_ is not a multiple of 4.");
		abort();
#endif
		if (( ((unsigned long int)u) & 0xf ) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)v) & 0xf) != 0 )
		{
			printf("float * const v is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)rexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		if (( ((unsigned long int)iexpansions) & 0xf) != 0 )
		{
			printf("float * const u is not 16-byte aligned!\n");
			abort();
		}
		
		const float * const yd = (float *)yd_;
		const __m128 xcenter = _mm_set_ps1(xcenter_);
		const __m128 ycenter = _mm_set_ps1(ycenter_);
		
		for(int dy=0; dy<_BLOCKSIZE_; dy++)
		{
			const int ybase = _BLOCKSIZE_*dy;
			const __m128 imagrp_ = _mm_load_ps(yd+dy*4) - ycenter;
			
			for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
			{
				const __m128 realrp_ = _mm_load_ps(xd+dx) - xcenter;
				
#ifdef _FMM_NOPRECDIV_KERNELS_
				const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
				const __m128 realrp = worse_division(realrp_, abs_rp);
				const __m128 imagrp = worse_division(imagrp_, abs_rp);
#else
				const __m128 abs_rp = (realrp_*realrp_ + imagrp_*imagrp_);
				const __m128 realrp = _mm_div_ps(realrp_, abs_rp);
				const __m128 imagrp = _mm_div_ps(imagrp_, abs_rp);
#endif			
				__m128 realsum = _mm_setzero_ps();
				__m128 imagsum = _mm_setzero_ps();
				
				__m128 realprod = _mm_set_ps(1,1,1,1);
				__m128 imagprod = _mm_setzero_ps();
				
				for(int e=0; e<_ORDER_; e+=4)
				{
					const __m128 rexps = _mm_load_ps(rexpansions+e);
					const __m128 iexps = _mm_load_ps(iexpansions+e);
					
					realsum = realsum + realprod*_select<0>(rexps) - imagprod*_select<0>(iexps);
					imagsum = imagsum + imagprod*_select<0>(rexps) + realprod*_select<0>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
					realsum = realsum + realprod*_select<1>(rexps) - imagprod*_select<1>(iexps);
					imagsum = imagsum + imagprod*_select<1>(rexps) + realprod*_select<1>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
					realsum = realsum + realprod*_select<2>(rexps) - imagprod*_select<2>(iexps);
					imagsum = imagsum + imagprod*_select<2>(rexps) + realprod*_select<2>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
					realsum = realsum + realprod*_select<3>(rexps) - imagprod*_select<3>(iexps);
					imagsum = imagsum + imagprod*_select<3>(rexps) + realprod*_select<3>(iexps);
					
					{
						const __m128 rold = realprod;
						const __m128 iold = imagprod;
						
						realprod = rold*realrp + iold*imagrp;
						imagprod = iold*realrp - rold*imagrp;
					}
					
				}
				
				_mm_store_ps(&u[dy][dx], _mm_add_ps(_mm_load_ps(&u[dy][dx]), imagsum*realrp - realsum*imagrp));
				_mm_store_ps(&v[dy][dx], _mm_add_ps(_mm_load_ps(&v[dy][dx]), imagsum*imagrp + realsum*realrp));
			}
		}
	}
	
	
	
	void indirect_interactions_cpp(const float xstart, const float ystart, const float h,
								   const float xcenter, const float ycenter,
								   const float * const rexpansions, const float * const iexpansions,
								   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_])
	{
		
		typedef std::complex<float> ExpansionsValueType;
		
		for(int dy=0; dy<_BLOCKSIZE_; dy++)
			for(int dx=0; dx<_BLOCKSIZE_; dx++)
			{
				const double location[2] = {
					xstart + dx*h,
					ystart + dy*h
				};
				
				const ExpansionsValueType rp = ExpansionsValueType(location[0],location[1])-ExpansionsValueType(xcenter, ycenter);
				
				std::complex<double> csum = std::complex<double>(0,0);
				std::complex<double> prod = std::complex<double>(1,0);
				
				for (int n=0;n<_ORDER_;++n)
				{
					csum+=(prod*std::complex<double>(rexpansions[n], iexpansions[n]));
					prod/=(std::complex<double>)rp;
				}
				
				csum*= -ExpansionsValueType(0,1)/rp;
				
				u[dy][dx] += csum.real();
				v[dy][dx] -= csum.imag();
				
			}
		
		/*
		 
		 
		 for(int dy=0; dy<_BLOCKSIZE_; dy++)
		 for(int dx=0; dx<_BLOCKSIZE_; dx++)
		 {
		 const double realrp_ = xstart + dx*h - xcenter;
		 const double imagrp_ = ystart + dy*h - ycenter;
		 const double invabs_rp = 1/(realrp_*realrp_ + imagrp_*imagrp_);
		 const double realrp = realrp_*invabs_rp;
		 const double imagrp = imagrp_*invabs_rp;
		 
		 
		 double realsum = 0, imagsum = 0;
		 double realprod = 1, imagprod = 0;
		 
		 
		 for(int i=0; i<_ORDER_; i++)
		 {
		 const double realval = rexpansions[i];
		 const double imagval = iexpansions[i];
		 
		 
		 realsum += realprod*realval - imagprod*imagval;
		 imagsum += imagprod*realval + realprod*imagval;
		 
		 
		 const double rold = realprod;
		 const double iold = imagprod;
		 
		 
		 realprod = rold*realrp + iold*imagrp;
		 imagprod = iold*realrp - rold*imagrp;
		 }
		 u[dy][dx] += imagsum*realrp - realsum*imagrp;
		 v[dy][dx] += imagsum*imagrp + realsum*realrp;
		 }
		 */
	}
}


#include <complex>
namespace AggressiveDiego
{
	using namespace std;
	
	void check_quality(const float xstart, const float ystart, const float h,
					   const float * const xcenter, const float * const ycenter,
					   const float * const rexpansions, const float * const iexpansions, const int nexpansions,
					   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const double tol)
	{
		double L1[2]= {0, 0};
		double LINF[2] = {0, 0};
		
		typedef std::complex<float> ExpansionsValueType;
		
		for(int dy=0; dy<_BLOCKSIZE_; dy++)
			for(int dx=0; dx<_BLOCKSIZE_; dx++)
			{
				double ucheck = 0, vcheck = 0;
				
				for(int e=0; e<nexpansions; e++)
				{
					const double location[2] = {
						xstart + dx*h,
						ystart + dy*h
					};
					
					const ExpansionsValueType rp = ExpansionsValueType(location[0],location[1])-ExpansionsValueType(xcenter[e], ycenter[e]);
					
					std::complex<double> csum = std::complex<double>(0,0);
					std::complex<double> prod = std::complex<double>(1,0);
					
					for (int n=0;n<_ORDER_;++n)
					{
						csum+=(prod*std::complex<double>(rexpansions[e*_ORDER_+n], iexpansions[e*_ORDER_+n]));
						prod/=(std::complex<double>)rp; 
					}
					
					csum*= -ExpansionsValueType(0,1)/rp;
					
					ucheck += csum.real();
					vcheck -= csum.imag();
				}
				
				const double uactual = u[dy][dx];
				const double vactual = v[dy][dx];
				
				const double diffu = (ucheck - uactual)/(max(1e-3, max(fabs(uactual), fabs(ucheck))));
				const double diffv = (vcheck - vactual)/(max(1e-3, max(fabs(vactual), fabs(vcheck))));
				
				assert(fabs(diffu) < tol);
				assert(fabs(diffv) < tol);
				
				L1[0] += fabs(diffu);
				L1[1] += fabs(diffv);
				
				LINF[0] = max(LINF[0], fabs(diffu));
				LINF[1] = max(LINF[1], fabs(diffv));
			}
		
		printf("\n====================================================================\n");
		printf("============     QUALITY CHECK (%.1e)   PASSED      =============\n", tol);
		printf("====================================================================\n");
		
		printf("relative linf-norm discrepancies: u: %e v: %e\n", LINF[0], LINF[1]);
		printf("relative l1-norm discrepancies: u: %e v: %e\n", L1[0], L1[1]);
		printf("average relative discrepancies: u: %e v: %e\n", L1[0]/(double)(_BLOCKSIZE_*_BLOCKSIZE_), L1[1]/(double)(_BLOCKSIZE_*_BLOCKSIZE_));
		printf("average relative expansion discrepancy: u: %e v: %e\n", L1[0]/(double)(_ORDER_*_BLOCKSIZE_*_BLOCKSIZE_*nexpansions), L1[1]/(double)(_ORDER_*_BLOCKSIZE_*_BLOCKSIZE_*nexpansions));
	}
}
