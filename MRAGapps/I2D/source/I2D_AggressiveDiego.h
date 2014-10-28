/*
 *  AggressiveDiego_KernelFMM.h
 *  I2D_CoreFMM_Aggressive_SSE_Kernel.h
 *
 *  Created by Diego Rossinelli on 12/24/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

namespace  AggressiveDiego {
	
	void direct_interactions_sse(const float x_start, const float y_start, const float h, 
								 const float * const xs, const float * const ys, const float * const ws, const int nsources,
								 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const bool safe);
	
	void direct_interactions_sse(float xd[_BLOCKSIZE_], float (* const yd)[4],
								 const float * const xs, const float * const ys, const float * const ws, const int nsources,
								 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const bool safe);
	
	void direct_interactions_cpp(const float x_start, const float y_start, const float h, 
								 const float * const xs, const float * const ys, const float * const ws, const int nsources,
								 float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const bool safe);
	
	void check_quality( const float x_start, const float y_start, const float h, 
					   const float * const xs, const float * const ys, const float * const ws, const int nsources,
					   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const double tol);
	
	
	void indirect_interactions_sse(const float xstart, const float ystart, const float h,
								   const float xcenter, const float ycenter,
								   const float * const rexpansions, const float * const iexpansions,
								   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_]);
	
	void indirect_interactions_sse(float xd[_BLOCKSIZE_], float (* const yd)[4],
								   const float xcenter, const float ycenter,
								   const float * const rexpansions, const float * const iexpansions,
								   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_]);
	
	void indirect_interactions_ssefloat(const float xstart, const float ystart, const float h,
										const float xcenter, const float ycenter,
										const float * const rexpansions, const float * const iexpansions,
										float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_]);
	
	void indirect_interactions_ssefloat(float xd[_BLOCKSIZE_], float (* const yd_)[4],
										const float xcenter_, const float ycenter_,
										const float * const rexpansions, const float * const iexpansions,
										float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_]);
	
	void indirect_interactions_cpp(const float xstart, const float ystart, const float h,
								   const float xcenter, const float ycenter,
								   const float * const rexpansions, const float * const iexpansions,
								   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_]);
	
	void check_quality(const float xstart, const float ystart, const float h,
					   const float * const xcenter, const float * const ycenter,
					   const float * const rexpansions, const float * const iexpansions, const int nexpansions,
					   float (* const u)[_BLOCKSIZE_], float (* const v)[_BLOCKSIZE_], const double tol);
}
