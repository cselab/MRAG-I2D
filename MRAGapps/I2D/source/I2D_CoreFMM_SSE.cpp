/*
 *  I2D_CoreFMM_SSE.cpp
 *  I2D_ROCKS
 *
 *  Created by Roman Schaerer on 12/26/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

extern  double _THETA;
#define _FMMSILENT

#include <cstdio>
#include <fstream>
#include <string>

#include "I2D_CoreFMM_SSE.h"

#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>

#include "I2D_Headers.h"
#include "I2D_Types.h"
#include "I2D_FMMTypes.h"

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGProfiler.h"

#include "mani-fmm2d/VortexExpansions.h"
#include "mani-fmm2d/hcfmm_box.h"
#include "mani-fmm2d/hcfmm_boxBuilder_serial.h"
#ifndef _MRAG_TBB
#include "mani-fmm2d/hcfmm_evaluator_serial.h"
#else
#include "mani-fmm2d/hcfmm_evaluator_tbb.h"
#endif

#include "I2D_CoreFMM_AggressiveVel.h"
#include "I2D_AggressiveDiego.h"
#include "I2D_CoreFMM_PlanBuilder.h"
#include "I2D_CoreFMM_PlanBuilderWim.h"

#include "I2D_CoreFMM_Plan.h"

#ifdef __INTEL_COMPILER
inline __m128 operator+(__m128 a, __m128 b){ return _mm_add_ps(a, b); }
inline __m128 operator*(__m128 a, __m128 b){ return _mm_mul_ps(a, b); }
inline __m128 operator-(__m128 a,  __m128 b){ return _mm_sub_ps(a, b); }
#endif

typedef tbb::tick_count ClockTime;
typedef HCFMM::Box<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBox;

struct FloatVelocityBlock
{
	float u[2][_BLOCKSIZE_][_BLOCKSIZE_];
	float ypos[_BLOCKSIZE_][4];
	float xpos[_BLOCKSIZE_];
	
	void initialize(const float xstart, const float ystart, const float h_gridpoint)
	{
		//clear u[0] and u[1]
		{
			const __m128 zero = _mm_setzero_ps();
			
			float * const myu = (float *)u[0];
			for(int dy=0; dy<_BLOCKSIZE_; dy++)
				for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
					_mm_store_ps(myu + dx + dy*_BLOCKSIZE_ , zero);
			
			float * const myv = (float *)u[1];
			for(int dy=0; dy<_BLOCKSIZE_; dy++)
				for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
					_mm_store_ps(myv + dx + dy*_BLOCKSIZE_ , zero);
		}
		
		//initialize xpos and ypos
		{
			const __m128 x0 = _mm_set_ps1(xstart);
			const __m128 y0 = _mm_set_ps1(ystart);
			const __m128 h = _mm_set_ps1(h_gridpoint);
			const __m128 d = _mm_set_ps(+3,+2,+1,+0);
			
			for(int dx=0; dx<_BLOCKSIZE_; dx += 4)
				_mm_store_ps(xpos + dx, (d+_mm_set_ps1(dx))*h + x0);
			
			for(int dy=0; dy<_BLOCKSIZE_; dy++)
				_mm_store_ps(ypos[dy], _mm_set_ps1(dy)*h + y0);
		}
	}
	
	static FloatVelocityBlock * allocate()
	{
		cache_aligned_allocator<FloatVelocityBlock> allocator;
		FloatVelocityBlock * ptr = allocator.allocate(1);
		return ptr;
	}
	
	static void deallocate(FloatVelocityBlock *& velblocks)
	{
		cache_aligned_allocator<FloatVelocityBlock> allocator;
		allocator.deallocate(velblocks, 1);
		velblocks = NULL;
	}
};

class DirectInfo {
public:
	DirectInfo () : xs (NULL), ys (NULL), ws (NULL), num_sources (-1), is_overlapping (false) {
	}
	
	DirectInfo (const float* const _xs, const float* const _ys, const float* const _ws, int _num_sources, bool _is_overlapping) :
	xs (_xs), ys (_ys), ws (_ws), num_sources (_num_sources), is_overlapping (_is_overlapping) {
	}
	
	DirectInfo (const DirectInfo& other) :
	xs (other.xs), ys (other.ys), ws (other.ws), num_sources (other.num_sources),
	is_overlapping (other.is_overlapping) {
	}
	
	const DirectInfo& operator=(const DirectInfo& other)
	{
		xs = other.xs;
		ys = other.ys;
		ws = other.ws;
		num_sources = other.num_sources;
		is_overlapping = other.is_overlapping;
		
		return *(this);
	}
	
	const float* xs;
	const float* ys;
	const float* ws;
	int num_sources;
	bool is_overlapping;
};

class IndirectInfo {
public:
	IndirectInfo () : xbox (-1), ybox (-1), real_values (NULL), imag_values (NULL) {
	}
	
	IndirectInfo (float _xbox, float _ybox, const float* const _real_values, const float* const _imag_values) :
	xbox (_xbox), ybox (_ybox), real_values (_real_values), imag_values (_imag_values) {
	}
	
	IndirectInfo (const IndirectInfo& other) :
	xbox (other.xbox), ybox (other.ybox), real_values (other.real_values), imag_values (other.imag_values) {
	}
	
	const IndirectInfo& operator=(const IndirectInfo& other)
	{
		xbox = other.xbox;
		ybox = other.ybox;
		real_values = other.real_values;
		imag_values = other.imag_values;
		
		return *(this);
	}
	
	float xbox, ybox;
	const float* real_values;
	const float* imag_values;
};

class TargetInfo {
public:
	TargetInfo () : xstart (-1), ystart (-1), h (-1) {
	}
	
	TargetInfo (float _xstart, float _ystart, float _h) : xstart (_xstart), ystart (_ystart), h (_h) {
	}
	
	float xstart, ystart, h;
};

typedef std::vector <DirectInfo, tbb::scalable_allocator<DirectInfo > > DirectInfoVector;
typedef std::vector <IndirectInfo, tbb::scalable_allocator<IndirectInfo > > IndirectInfoVector;
typedef std::vector <std::pair <int,int>,  tbb::scalable_allocator<std::pair <int,int> > > IntervalVector;
typedef std::vector <IntervalVector, tbb::scalable_allocator <IntervalVector> > IntervalVectors;
typedef std::vector <std::vector <bool, tbb::scalable_allocator <bool> >, tbb::scalable_allocator <std::vector <bool, tbb::scalable_allocator <bool> > > > OverlappingVectors;
typedef std::list <std::pair <int,int>, tbb::scalable_allocator <std::pair<int,int> > > IntervalList;

class PlanPerBlock {
public:
	PlanPerBlock () : m_velocity_block (NULL) {}
	
	DirectInfoVector m_direct_infos;
	IndirectInfoVector m_indirect_infos;
	
	VelocityBlock* m_velocity_block;
};

template< typename T >
class range
{
public:
    typedef T value_type;
	
    range( T const & center ) : min_( center ), max_( center ) {}
    range( T const & min, T const & max )
	: min_( min ), max_( max ) {}
    T min() const { return min_; }
    T max() const { return max_; }
private:
    T min_;
    T max_;
};

template <typename T>
struct left_of_range : public std::binary_function< range<T>, range<T>, bool >
{
    bool operator()( range<T> const & lhs, range<T> const & rhs ) const
    {
        return lhs.min() < rhs.min()
		&& lhs.max() < rhs.min();
    }
};

typedef std::map< range<int>, int, left_of_range<int> > TranslationMap;

template <typename T>
class IntervalComparisonByStartPoint {
public:
	bool operator() (const std::pair <T,T>& interval_a, const std::pair <T,T>& interval_b) const {
		return (interval_a.first < interval_b.first);
	}
};


class TimeInfo {
public:
	TimeInfo () : direct_reduction_time (0), indirect_reduction_time (0),
	direct_total_gflop (0), indirect_total_gflop (0),
	direct_avg_gflops(0), indirect_avg_gflops(0), 
	direct_kernel_calls (0), indirect_kernel_calls (0){}
	
	TimeInfo (double _direct_reduction_time, double _indirect_reduction_time,
			  double _direct_total_gflop, double _indirect_total_gflop,
			  double _direct_avg_gflops, double _indirect_avg_gflops,
			  unsigned int _direct_kernel_calls, unsigned int _indirect_kernel_calls) :
	direct_reduction_time (_direct_reduction_time), indirect_reduction_time (_indirect_reduction_time),
	direct_total_gflop (_direct_total_gflop), indirect_total_gflop (_indirect_total_gflop),
	direct_avg_gflops(_direct_avg_gflops), indirect_avg_gflops(_indirect_avg_gflops), 
	direct_kernel_calls (_direct_kernel_calls), indirect_kernel_calls (_indirect_kernel_calls) {}
	
	double direct_reduction_time, indirect_reduction_time;
	double direct_total_gflop, indirect_total_gflop;
	double direct_avg_gflops, indirect_avg_gflops; 
	double direct_kernel_calls, indirect_kernel_calls;
};

typedef std::vector <TimeInfo> TimeInfos;

class SSE_Plan : public Plan {
public:
	
	SSE_Plan (int _num_target_blocks, BlockInfo* _target_blocks, int num_source_particles, tBox* root_node, bool _b_verbose) :
	m_num_target_blocks (_num_target_blocks), m_target_blocks (_target_blocks),
	m_plan_per_block (_num_target_blocks), m_num_source_particles (num_source_particles),
	m_direct_overlapping_intervals (_num_target_blocks), m_direct_non_overlapping_intervals (_num_target_blocks),
	m_root_source_particles (&root_node->vparticles[0]),
	m_indirect_interactions (_num_target_blocks),
	m_direct_buffer_size (0), m_indirect_buffer_size (0), b_verbose (_b_verbose),
	m_xs (NULL), m_ys (NULL), m_ws (NULL), m_real_values (NULL), m_imag_values (NULL) {
		
	}
	
	~SSE_Plan () {
		if (m_direct_buffer_size > 0) {
			assert (m_xs != NULL && m_ws != NULL && m_ws != NULL);
			_mm_free (m_xs);
			_mm_free (m_ys);
			_mm_free (m_ws);
		}
		if (m_indirect_buffer_size > 0) {
			assert (m_real_values != NULL && m_imag_values != NULL);
			_mm_free (m_real_values);
			_mm_free (m_imag_values);
		}
	}
	
	void merge_direct_intervals (int target_block) {
		IntervalVector* intervals;
		
		for (int i=0; i<=1; ++i) {
			if (i == 0)
				intervals = &m_direct_non_overlapping_intervals [target_block];
			else
				intervals = &m_direct_overlapping_intervals [target_block];
			
			std::sort (intervals->begin (), intervals->end (), IntervalComparisonByStartPoint <int> ());
			
			if (!intervals->empty ()) {
				IntervalList merged_intervals;
				IntervalList::iterator merge_iter;
				IntervalVector::iterator all_iter;
				
				merged_intervals.push_back (*intervals->begin ());
				merge_iter = merged_intervals.begin ();
				all_iter = intervals->begin (); ++all_iter;
				while (all_iter != intervals->end ()) {
					assert (merge_iter != merged_intervals.end());
					if (all_iter->first <= merge_iter->second+1) {
						merge_iter->second = std::max (merge_iter->second, all_iter->second);
						++all_iter;
					} else {
						merged_intervals.push_back (*all_iter);
						++all_iter;
						++merge_iter;
					}
				}
				
				if (merged_intervals.size () != intervals->size ()) {
					intervals->resize (merged_intervals.size ());
					
					merge_iter = merged_intervals.begin ();
					for (int i=0; merge_iter != merged_intervals.end (); ++merge_iter, ++i) {
						(*intervals) [i] = *merge_iter;
					}
				}
			}
		}
	}
	
	void extract_direct_data () 
	{
		IntervalVector all_intervals;
		
#ifdef _I2D_MPI_
		
		int size = 0;
		for (int i=0; i<m_num_target_blocks; ++i) 
			size += m_direct_overlapping_intervals [i].size ()+m_direct_non_overlapping_intervals [i].size ();
		
		all_intervals.resize(size);
		
		if (all_intervals.size () == 0) 
		{
			printf("Akamon now.\n");
			abort();
			return;
		}
		
		// copy data from m_direct_intervals into all_intervals
		{
			int index = 0;
			for (int i=0; i<m_num_target_blocks; ++i) {
				for (int j=0;j<m_direct_overlapping_intervals [i].size (); ++j) {
					all_intervals [index+j] = m_direct_overlapping_intervals [i][j];
				}
				index += m_direct_overlapping_intervals [i].size ();
				for (int j=0;j<m_direct_non_overlapping_intervals [i].size (); ++j) {
					all_intervals [index+j] = m_direct_non_overlapping_intervals [i][j];
				}
				index += m_direct_non_overlapping_intervals [i].size ();
			}
		}
		
		// perform parallel sort of intervals by sorting according to their starting value
		tbb::parallel_sort (all_intervals.begin (), all_intervals.end (), IntervalComparisonByStartPoint <int> ());
		
#else
		all_intervals.push_back(pair<int, int>(0, m_num_source_particles-1));
#endif
		
		IntervalVector merged_intervals;
		IntervalVector::iterator all_iter = all_intervals.begin();
		
		merged_intervals.push_back (*all_iter);
		++all_iter;
		
		for( ; all_iter!=all_intervals.end(); ++all_iter)
		{
			pair<int,int>& interval = merged_intervals.back();
			
			if (all_iter->first <= interval.second+1) 
				interval.second = std::max (interval.second, all_iter->second);
			else 
				merged_intervals.push_back (*all_iter);
		}
		
		// construct translation map
		{
			int offset = 0;
			int last_endpoint = 0;
			for (IntervalVector::iterator it = merged_intervals.begin (); it != merged_intervals.end (); ++it) 
			{
				offset += it->first-last_endpoint;
				last_endpoint = it->second;
				
				m_translation_map [range<int>(it->first, it->second)] = offset;
			}
			assert (m_translation_map.size () > 0);
			
			m_direct_buffer_size = last_endpoint - offset + 1;
			assert (m_direct_buffer_size > 0);
		}
		
		//allocate memory
		assert (m_xs == NULL && m_ys == NULL && m_ws == NULL);
		m_xs = (float*)_mm_malloc(sizeof(float)*m_direct_buffer_size, 16);
		m_ys = (float*)_mm_malloc(sizeof(float)*m_direct_buffer_size, 16);
		m_ws = (float*)_mm_malloc(sizeof(float)*m_direct_buffer_size, 16);
		
		//copy data
		TranslationMap::const_iterator map_iter;
		for (map_iter = m_translation_map.begin (); map_iter != m_translation_map.end (); ++map_iter) {
			
			const int old_start = map_iter->first.min ();
			const int old_end = map_iter->first.max ();
			const int new_start = old_start - map_iter->second;
			const int new_end = old_end - map_iter->second;
			
			assert (old_start <= old_end && old_start >= 0 && old_end < m_num_source_particles);
			assert (new_start >= 0 && new_end < m_direct_buffer_size);
			
			for (int i = new_start, j = old_start; i<= old_end; ++i, ++j) 
			{
				m_xs [i] = (m_root_source_particles+j)->x[0];
				m_ys [i] = (m_root_source_particles+j)->x[1];
				m_ws [i] = (m_root_source_particles+j)->w[0];
			}
		}
	}
	
	void make_direct_plan () 
	{
		for (int i=0; i<m_num_target_blocks; ++i) 
		{
			assert(m_plan_per_block [i].m_direct_infos.size()==0);
			
			for (unsigned int j=0; j<m_direct_overlapping_intervals [i].size (); ++j)
			{
				TranslationMap::const_iterator map_iter = m_translation_map.find (m_direct_overlapping_intervals [i][j].first);
				
				assert (map_iter != m_translation_map.end ());
				assert(map_iter->first.min() <= m_direct_overlapping_intervals [i][j].first);
				assert(map_iter->first.max() >= m_direct_overlapping_intervals [i][j].first);
				assert(map_iter->first.min() <= m_direct_overlapping_intervals [i][j].second);
				assert(map_iter->first.max() >= m_direct_overlapping_intervals [i][j].second);
				
				const int index = m_direct_overlapping_intervals [i][j].first - map_iter->second;
				
				const int num_particles = m_direct_overlapping_intervals [i][j].second - m_direct_overlapping_intervals [i][j].first + 1;
				assert (num_particles > 0);
				
				assert (index >= 0 && index < m_direct_buffer_size);
				assert (index + num_particles - 1 < m_direct_buffer_size);
				
				const float* const xsp = &m_xs [index];
				const float* const ysp = &m_ys [index];
				const float* const wsp = &m_ws [index];
				
				DirectInfo direct_info (xsp, ysp, wsp, num_particles, true);
				
				m_plan_per_block[i].m_direct_infos.push_back(direct_info);
			}
			
			for (unsigned int j=0; j<m_direct_non_overlapping_intervals [i].size (); ++j)
			{
				TranslationMap::const_iterator map_iter = m_translation_map.find (m_direct_non_overlapping_intervals [i][j].first);
				
				assert (map_iter != m_translation_map.end ());
				assert(map_iter->first.min() <= m_direct_non_overlapping_intervals [i][j].first);
				assert(map_iter->first.max() >= m_direct_non_overlapping_intervals [i][j].first);
				assert(map_iter->first.min() <= m_direct_non_overlapping_intervals [i][j].second);
				assert(map_iter->first.max() >= m_direct_non_overlapping_intervals [i][j].second);
				
				const int index = m_direct_non_overlapping_intervals [i][j].first - map_iter->second;
				
				const int num_particles = m_direct_non_overlapping_intervals [i][j].second - m_direct_non_overlapping_intervals [i][j].first + 1;
				
				assert (num_particles > 0);
				assert (index >= 0 && index < m_direct_buffer_size);
				assert (index + num_particles - 1 < m_direct_buffer_size);
				
				const float* const xsp = &m_xs [index];
				const float* const ysp = &m_ys [index];
				const float* const wsp = &m_ws [index];
				
				DirectInfo direct_info (xsp, ysp, wsp, num_particles, true);
				
				m_plan_per_block[i].m_direct_infos.push_back(direct_info);
			}
		}
	}
	
	void addDirectInteraction (const tBox* const source_box, int target_id, bool is_overlapping){
		assert (source_box->vparticles[0].x[0] == m_root_source_particles [&source_box->vparticles[0] - m_root_source_particles].x[0]);
		
		int left_index = (&source_box->vparticles[0]-m_root_source_particles);
		int right_index = (&source_box->vparticles[source_box->nParticles-1]-m_root_source_particles);
		if (is_overlapping) {
			assert (target_id < m_direct_overlapping_intervals.size ());
			m_direct_overlapping_intervals [target_id].push_back (std::pair <int,int> (left_index, right_index));
		}
		else {
			assert (target_id < m_direct_non_overlapping_intervals.size ());
			m_direct_non_overlapping_intervals [target_id].push_back (std::pair <int,int> (left_index, right_index));
		}
		
	}
	
	void addIndirectInteraction (const tBox* const source_box, int target_id) 
	{
		assert (target_id < m_indirect_interactions.size ());
		m_indirect_interactions [target_id].push_back (source_box);
	}
	
	struct CreateExpansionPlan
	{
		const float * const expansion_reals;
		const float * const expansion_imags;
		
		const std::vector <std::vector <const tBox*> >& indirectinfo;
		const std::map <int, std::pair <const tBox*, int> >& id2index;
		
		std::vector <PlanPerBlock>& m_plan_per_block;
		
		CreateExpansionPlan(const CreateExpansionPlan& other):
		expansion_reals(other.expansion_reals), expansion_imags(other.expansion_imags),
		indirectinfo(other.indirectinfo), id2index(other.id2index), m_plan_per_block(other.m_plan_per_block){}
		
		CreateExpansionPlan(const float * const expansion_reals_, const float * const expansion_imags_, 
							const std::vector <std::vector <const tBox*> >& indirectinfo_, const  std::map <int, std::pair <const tBox*, int> >& id2index_,
							std::vector <PlanPerBlock>& m_plan_per_block_): 
		expansion_reals(expansion_reals_), expansion_imags(expansion_imags_),
		indirectinfo(indirectinfo_),id2index(id2index_),
		m_plan_per_block(m_plan_per_block_){}
		
		void operator()(blocked_range<int> range) const
		{
			for (int i=range.begin();i<range.end();++i) 
			{
				const std::vector <const tBox*>& myinfo = indirectinfo[i];
				
				for(std::vector <const tBox*>::const_iterator it = myinfo.begin(); it != myinfo.end(); ++it)
				{
					assert( *it != NULL);
					const tBox& box = **it;
					
					const int boxID = box.getId();
					std::map <int, std::pair <const tBox*, int> >::const_iterator itIndex = id2index.find (boxID);
					assert (itIndex != id2index.end ());
					
					const int index = itIndex->second.second;
					
					IndirectInfo indirect_info (box.expansions.Center[0],
												box.expansions.Center[1],
												&expansion_reals [index],  &expansion_imags [index]);
					
					m_plan_per_block[i].m_indirect_infos.push_back(indirect_info);
				}
			}
		}
	};
	
	void make_indirect_plan () 
	{
		for (int i=0;i<m_num_target_blocks;++i) {
			m_plan_per_block [i].m_velocity_block = (VelocityBlock*)m_target_blocks[i].ptrBlock;
			assert (m_plan_per_block [i].m_velocity_block != NULL);
		}
		
		std::map <int, std::pair <const tBox*, int> > m_indirect_source_id2index;
		
		{
			int ctr = 0;
			for (int i=0; i<m_num_target_blocks; ++i) 
			{
				const std::vector <const tBox*>::const_iterator START = m_indirect_interactions[i].begin();
				const std::vector <const tBox*>::const_iterator END = m_indirect_interactions[i].end();
				
				for (std::vector <const tBox*>::const_iterator it = START; it != END; ++it) 
				{
					const tBox* source_box = *it;
					
					std::map <int, std::pair <const tBox*, int> >::iterator map_it = m_indirect_source_id2index.find (source_box->getId());
					
					if (map_it == m_indirect_source_id2index.end ()) 
					{
						m_indirect_source_id2index.insert (std::pair <int, std::pair <const tBox*, int> > (source_box->id, std::pair <const tBox*, int> (source_box, ctr*_ORDER_)));
						ctr++;
					} 
					else 
					{
						assert (map_it->second.first == source_box);
						assert (map_it->second.second < ctr*_ORDER_);
					}
				}
			}
			
			m_indirect_buffer_size = m_indirect_source_id2index.size()*_ORDER_;
			assert(m_indirect_buffer_size == ctr*_ORDER_);
		}
		
		if (m_indirect_buffer_size == 0)
			return;
		
		assert (m_real_values == NULL && m_imag_values == NULL);
		
		m_real_values = (float*)_mm_malloc(sizeof(float)*m_indirect_buffer_size, 16);
		m_imag_values = (float*)_mm_malloc(sizeof(float)*m_indirect_buffer_size, 16);
		
		// loop over all elements in m_indirect_interactions
		CreateExpansionPlan create_indirect_plan(m_real_values, m_imag_values, m_indirect_interactions, m_indirect_source_id2index, m_plan_per_block);
		tbb::parallel_for (tbb::blocked_range<int> (0,m_num_target_blocks), create_indirect_plan, auto_partitioner ());
		
		// loop over m_indirect_source_id2index
		for (std::map <int, std::pair <const tBox*, int> >::const_iterator it_id2index = m_indirect_source_id2index.begin(); it_id2index != m_indirect_source_id2index.end (); ++it_id2index) 
		{
			const tBox * source_box = it_id2index->second.first;
			const int start_index = it_id2index->second.second;
			
			assert (start_index >= 0);
			assert(start_index+_ORDER_ <= m_indirect_source_id2index.size()*_ORDER_);
			
			for (int j=0; j<_ORDER_; j++) 
			{
				m_real_values [j+start_index] = source_box->expansions.values[0][j].real(); 
				m_imag_values [j+start_index] = source_box->expansions.values[0][j].imag(); 
			}
		}
	}
	
	Real m_scale_factor;
	
	const int m_num_target_blocks;
	BlockInfo* m_target_blocks;
	const int m_num_source_particles;
	VelocitySourceParticle* const m_root_source_particles;
	
	
	int m_direct_buffer_size, m_indirect_buffer_size;
	float* m_xs, *m_ys, *m_ws;
	float* m_real_values, *m_imag_values;
	
	std::vector <std::vector <const tBox*> > m_indirect_interactions;
	
	std::vector <PlanPerBlock> m_plan_per_block;
	IntervalVectors m_direct_overlapping_intervals, m_direct_non_overlapping_intervals;
	TranslationMap m_translation_map;
	bool b_verbose;
};

struct ReduceDirectSourceContributions {
public:
	
	void operator () (tbb::blocked_range <int> range) 
	{
		ClockTime start = tick_count::now ();
		
		double gflop = 0;
		
		for (int i=range.begin (); i!=range.end(); ++i) 
		{
			assert (source_infos[i].xs != NULL);
			assert (source_infos[i].ys != NULL);
			assert (source_infos[i].ws != NULL);
			
			if (source_infos[i].num_sources < 4) 
			{
				AggressiveDiego::direct_interactions_cpp (target_info.xstart, target_info.ystart, target_info.h,
														  source_infos[i].xs, source_infos[i].ys, source_infos[i].ws, source_infos[i].num_sources,
														  velocity_block->u[0], velocity_block->u[1], source_infos[i].is_overlapping);
				
			} 
			else 
			{
				AggressiveDiego::direct_interactions_sse (velocity_block->xpos, velocity_block->ypos,
														  source_infos[i].xs, source_infos[i].ys, source_infos[i].ws, source_infos[i].num_sources,
														  velocity_block->u[0], velocity_block->u[1], source_infos[i].is_overlapping);
			}
			
			gflop += source_infos[i].num_sources * _BLOCKSIZE_ * _BLOCKSIZE_ * 14.*1e-9;
		}
		
		ClockTime end = tick_count::now ();		
		const double time = (end-start).seconds();
		elapsed_time += (end-start).seconds();
		kernel_calls += range.size();
		total_gflop += gflop;
		avg_gflops += gflop/time;
	}
	
	void _setup()
	{
		velocity_block = FloatVelocityBlock::allocate();
		assert (velocity_block != NULL);
		velocity_block->initialize(target_info.xstart, target_info.ystart, target_info.h);
	}
	
	ReduceDirectSourceContributions (const DirectInfoVector& _source_infos, const TargetInfo& _target_info) :
	velocity_block (NULL), source_infos (_source_infos), target_info (_target_info), elapsed_time (0), total_gflop (0), kernel_calls (0), avg_gflops(0)
	{
		_setup();
	}
	
	
	~ReduceDirectSourceContributions () 
	{
		FloatVelocityBlock::deallocate(velocity_block);
	}
	
	ReduceDirectSourceContributions (ReduceDirectSourceContributions& other, tbb::split) :
	velocity_block (NULL), source_infos (other.source_infos), target_info (other.target_info), elapsed_time (0), total_gflop (0), kernel_calls (0), avg_gflops(0)
	{
		_setup();
	}
	
	void join (ReduceDirectSourceContributions& other) 
	{
		for (int d=0; d<2; ++d) {
			for (int y=0; y<_BLOCKSIZE_; ++y) {
				for (int x=0; x<_BLOCKSIZE_; ++x) {
					velocity_block->u[d][y][x] += other.velocity_block->u[d][y][x];
				}
			}
		}
		
		elapsed_time += other.elapsed_time;
		total_gflop += other.total_gflop;
		avg_gflops += other.avg_gflops;
		kernel_calls += other.kernel_calls;
	}
	
	const DirectInfoVector & source_infos;
	const TargetInfo target_info;
	FloatVelocityBlock * velocity_block;
	double elapsed_time, total_gflop, avg_gflops;
	unsigned int kernel_calls;
};

struct ReduceIndirectSourceContributions {
public:
	
	void operator () (tbb::blocked_range <int> range) 
	{
		ClockTime start = tick_count::now ();
		
		for (int i=range.begin (); i!=range.end(); ++i) 
		{
			assert (source_infos[i].real_values != NULL);
			assert (source_infos[i].imag_values != NULL);
			
#ifdef _FMM_MIXEDPREC_KERNELS_
			AggressiveDiego::indirect_interactions_sse (velocity_block->xpos, velocity_block->ypos,
														source_infos[i].xbox, source_infos[i].ybox, source_infos[i].real_values, source_infos[i].imag_values, velocity_block->u[0], velocity_block->u[1]);
#else
			AggressiveDiego::indirect_interactions_ssefloat (velocity_block->xpos, velocity_block->ypos,
															 source_infos[i].xbox, source_infos[i].ybox, source_infos[i].real_values, source_infos[i].imag_values, velocity_block->u[0], velocity_block->u[1]);
			
			//AggressiveDiego::indirect_interactions_sse (target_info.xstart, target_info.ystart, target_info.h,
			//											source_infos[i].xbox, source_infos[i].ybox, source_infos[i].real_values, source_infos[i].imag_values, velocity_block->u[0], velocity_block->u[1]);
#endif
			//	AggressiveDiego::indirect_interactions_cpp (target_info.xstart, target_info.ystart, target_info.h,
			//			source_infos[i].xbox, source_infos[i].ybox, source_infos[i].real_values, source_infos[i].imag_values, velocity_block->u[0], velocity_block->u[1]);
		}
		
		ClockTime end = tick_count::now ();		
		const double time = (end-start).seconds();
		
		const double gflop = 1e-9*(double)(_BLOCKSIZE_ * _BLOCKSIZE_ * (23 + 14*_ORDER_))*range.size();
		total_gflop += gflop;
		kernel_calls += range.size();
		elapsed_time += (end-start).seconds();
		avg_gflops += gflop/time;
	}
	
	void _setup()
	{
		velocity_block = FloatVelocityBlock::allocate();
		assert (velocity_block != NULL);
		velocity_block->initialize(target_info.xstart, target_info.ystart, target_info.h);
	}
	
	ReduceIndirectSourceContributions (const IndirectInfoVector& _source_infos, const TargetInfo& _target_info) :
	velocity_block (NULL), source_infos (_source_infos), target_info (_target_info), elapsed_time (0), total_gflop (0), kernel_calls (0), avg_gflops(0) {
		_setup();
	}
	
	~ReduceIndirectSourceContributions () {
		FloatVelocityBlock::deallocate(velocity_block);
	}
	
	ReduceIndirectSourceContributions (const ReduceIndirectSourceContributions& other, tbb::split) :
	velocity_block (NULL), source_infos (other.source_infos), target_info (other.target_info), elapsed_time (0), total_gflop (0), kernel_calls (0), avg_gflops(0) {
		_setup();
	}
	
	void join (ReduceIndirectSourceContributions& other) {
		for (int d=0; d<2; ++d) {
			for (int y=0; y<_BLOCKSIZE_; ++y) {
				for (int x=0; x<_BLOCKSIZE_; ++x) {
					velocity_block->u[d][y][x] += (double)other.velocity_block->u[d][y][x];
				}
			}
		}
		
		elapsed_time += other.elapsed_time;
		total_gflop += other.total_gflop;
		avg_gflops += other.avg_gflops;
		kernel_calls += other.kernel_calls;
	}
	
	const IndirectInfoVector& source_infos;
	const TargetInfo target_info;
	FloatVelocityBlock * velocity_block;
	double elapsed_time, total_gflop, avg_gflops;
	unsigned int kernel_calls;
};

template< int stage >
class ComputeTargetBlocks {
public:
	void operator () (tbb::blocked_range<int>& range) const {
		
		for (int i=range.begin(); i != range.end(); ++i) {
			assert (i >= 0 && i < m_plan.m_num_target_blocks);
			
			double culo[2] ;
			m_plan.m_target_blocks[i].pos(culo, 0, 0);
			
			TargetInfo target_info (culo[0], culo[1], m_plan.m_target_blocks [i].h [0]);
			
			VelocityBlock* const velocity_block = m_plan.m_plan_per_block [i].m_velocity_block;
			assert(velocity_block != NULL);
			
			switch (stage)
			{
				case 0:
				{
					ReduceDirectSourceContributions reduce_direct (m_plan.m_plan_per_block[i].m_direct_infos, target_info);
					const int ndirect = m_plan.m_plan_per_block [i].m_direct_infos.size ();
					parallel_reduce (tbb::blocked_range <int> (0,ndirect), reduce_direct, auto_partitioner ());
					
					for (int d=0; d<2; ++d) 
						for (int y=0; y<_BLOCKSIZE_; ++y) 
							for (int x=0; x<_BLOCKSIZE_; ++x) 
								velocity_block->u [d][y][x] = m_plan.m_scale_factor * reduce_direct.velocity_block->u [d][y][x];
					
					(*time_infos)[i].direct_reduction_time = reduce_direct.elapsed_time;
					(*time_infos)[i].direct_total_gflop = reduce_direct.total_gflop;
					(*time_infos)[i].direct_avg_gflops = reduce_direct.avg_gflops;
					(*time_infos)[i].direct_kernel_calls = reduce_direct.kernel_calls;
				}
					break;
					
				case 1:
				{
					ReduceIndirectSourceContributions reduce_indirect (m_plan.m_plan_per_block[i].m_indirect_infos, target_info);
					const int nindirect = m_plan.m_plan_per_block [i].m_indirect_infos.size ();
					parallel_reduce (tbb::blocked_range <int> (0,nindirect), reduce_indirect, auto_partitioner ());
					
					for (int d=0; d<2; ++d) 
						for (int y=0; y<_BLOCKSIZE_; ++y) 
							for (int x=0; x<_BLOCKSIZE_; ++x) 
								velocity_block->u [d][y][x] += m_plan.m_scale_factor * reduce_indirect.velocity_block->u [d][y][x];
					
					(*time_infos)[i].indirect_reduction_time = reduce_indirect.elapsed_time;
					(*time_infos)[i].indirect_total_gflop = reduce_indirect.total_gflop;
					(*time_infos)[i].indirect_avg_gflops = reduce_indirect.avg_gflops;
					(*time_infos)[i].indirect_kernel_calls = reduce_indirect.kernel_calls;
				}
					break;
			}
		}
	}
	
	ComputeTargetBlocks (ComputeTargetBlocks& other, tbb::split) :
	m_plan (other.m_plan), b_verbose(other.b_verbose), time_infos (other.time_infos) {
	}
	
	ComputeTargetBlocks (const ComputeTargetBlocks& other) :
	m_plan (other.m_plan), b_verbose(other.b_verbose), time_infos (other.time_infos) {
	}
	
	~ComputeTargetBlocks () {
		
	}
	
	ComputeTargetBlocks (const SSE_Plan& _plan, bool _b_verbose, TimeInfos* const _time_infos) :
	m_plan (_plan), b_verbose (_b_verbose), time_infos (_time_infos) {
	}
	
	TimeInfos* const time_infos;
	const SSE_Plan& m_plan;
	bool b_verbose;
};


class VelocityEvaluatorSSE {
public:
	
	typedef HCFMM::Box<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBox;
	typedef HCFMM::boxBuilder_serial<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBoxBuilder;
	
	BlockInfo * const m_target_blocks;
	const int m_num_target_blocks;
	tBox * const m_root_node;
	const Real inv_scaling;
	SSE_Plan m_plan;
	const int m_num_source_particles;
	bool b_verbose;
	TimeInfos m_time_infos;
	
	VelocityEvaluatorSSE (BlockInfo* target_blocks, int num_target_blocks, tBox * root_node,
						  const Real inv_scaling, int num_particles, bool _b_verbose) :
	m_target_blocks (target_blocks), m_num_target_blocks (num_target_blocks),
	m_plan (num_target_blocks, target_blocks, num_particles, root_node, _b_verbose),
	m_root_node (root_node), inv_scaling(inv_scaling), m_num_source_particles (num_particles),
	b_verbose (_b_verbose), m_time_infos (num_target_blocks) {
		assert (m_target_blocks != NULL);
		assert (m_num_target_blocks >= 0);
		assert (m_root_node != NULL);
		
		m_plan.m_scale_factor = 1./(2.0*M_PI)*inv_scaling;
		
		for (int i=0;i<num_target_blocks;++i) {
			assert (target_blocks[i].ptrBlock != NULL);
		}
	}
	
	void extract_interaction_data () {
//		PlanBuilder plan_builder (&m_plan, m_root_node, m_target_blocks, m_num_target_blocks);
		PlanBuilderWim plan_builder (&m_plan, m_root_node, m_target_blocks, m_num_target_blocks);
		plan_builder.run ();
		
	}
	
	void extract_direct_data () {
		m_plan.extract_direct_data ();
	}
	void make_direct_plan () {
		m_plan.make_direct_plan ();
	}
	void make_indirect_plan () {
		m_plan.make_indirect_plan ();
	}
	
	void run () {
		tbb::parallel_for (tbb::blocked_range<int> (0,m_num_target_blocks), ComputeTargetBlocks<0> (m_plan, b_verbose, &m_time_infos), auto_partitioner ());
		tbb::parallel_for (tbb::blocked_range<int> (0,m_num_target_blocks), ComputeTargetBlocks<1> (m_plan, b_verbose, &m_time_infos), auto_partitioner ());
	}
	
};

void I2D_CoreFMM_SSE::solve(const Real theta, const Real inv_scaling, BlockInfo * dest, const int nblocks, VelocitySourceParticle * srcparticles, const int nparticles) {
	
	typedef HCFMM::Box<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBox;
	typedef HCFMM::boxBuilder_serial<_VortexExpansions<VelocitySourceParticle, _ORDER_>,_FMM_MAX_LEVEL_> tBoxBuilder;
	
	Profiler profiler;
	
	_THETA = theta;
	
	tBox * rootBox= new tBox;
	ClockTime start_before_tree = tick_count::now ();
	profiler.push_start("tree");
	tBoxBuilder::buildBoxes(srcparticles, nparticles, rootBox);
	profiler.pop_stop();
	
	ClockTime start_before_expansions = tick_count::now ();
	profiler.push_start("expansions");
	tBoxBuilder::generateExpansions(rootBox);
	profiler.pop_stop();
	
	ClockTime start_before_plan = tick_count::now ();
	VelocityEvaluatorSSE evaluator(dest, nblocks, rootBox, inv_scaling, nparticles, b_verbose);
	profiler.push_start("extract interaction data");
	evaluator.extract_interaction_data ();
	profiler.pop_stop();
	profiler.push_start("extract direct data");
	evaluator.extract_direct_data ();
	profiler.pop_stop();
	profiler.push_start("make direct plan");
	evaluator.make_direct_plan ();
	profiler.pop_stop();
	profiler.push_start("make indirect plan");
	evaluator.make_indirect_plan ();
	profiler.pop_stop();
	profiler.push_start("evaluations");
	ClockTime start_before_run = tick_count::now ();
	evaluator.run ();
	ClockTime end = tick_count::now ();
	profiler.pop_stop();
	
	if(b_verbose || timestamp % 5 == 0) profiler.printSummary();
	
	if((b_verbose && timestamp % 5 == 0) || timestamp % 10 == 0)
	{
		double total_direct_time = 0, total_indirect_time = 0;
		double total_direct_particles = 0, total_indirect_boxes = 0;
		double num_direct_kernel_calls = 0, num_indirect_kernel_calls = 0;
		double direct_total_gflop = 0, indirect_total_gflop = 0;	
		double direct_avg_gflops = 0, indirect_avg_gflops = 0;
		
		unsigned int num_direct_cpp_evals = 0;
		unsigned int measured_direct_kernel_calls = 0, measured_indirect_kernel_calls = 0;
		
		for (int i=0;i<evaluator.m_num_target_blocks;++i) 
		{
			total_direct_time += evaluator.m_time_infos[i].direct_reduction_time;
			total_indirect_time += evaluator.m_time_infos[i].indirect_reduction_time;
			total_indirect_boxes += evaluator.m_plan.m_plan_per_block[i].m_indirect_infos.size();
			
			for (unsigned int j=0; j<evaluator.m_plan.m_plan_per_block [i].m_direct_infos.size (); ++j)
			{
				total_direct_particles += evaluator.m_plan.m_plan_per_block [i].m_direct_infos [j].num_sources;
				
				if (evaluator.m_plan.m_plan_per_block [i].m_direct_infos [j].num_sources < 4)
					num_direct_cpp_evals += 1;
			}
			
			num_direct_kernel_calls += evaluator.m_plan.m_plan_per_block [i].m_direct_infos.size ();
			num_indirect_kernel_calls += evaluator.m_plan.m_plan_per_block [i].m_indirect_infos.size ();
			
			direct_total_gflop += evaluator.m_time_infos [i].direct_total_gflop;
			indirect_total_gflop += evaluator.m_time_infos [i].indirect_total_gflop;
			direct_avg_gflops += evaluator.m_time_infos [i].direct_avg_gflops;
			indirect_avg_gflops += evaluator.m_time_infos [i].indirect_avg_gflops;
			
			measured_direct_kernel_calls += evaluator.m_time_infos [i].direct_kernel_calls;
			measured_indirect_kernel_calls += evaluator.m_time_infos [i].indirect_kernel_calls;
		}
		
		direct_avg_gflops /= measured_direct_kernel_calls;
		indirect_avg_gflops /= measured_indirect_kernel_calls;
		
		const double total_time = total_direct_time + total_indirect_time;
		const double total_gflop = direct_total_gflop + indirect_total_gflop;
		
		const double average_num_particles_per_box = total_direct_particles/num_direct_kernel_calls;
		
		const double direct_average_gflops = direct_total_gflop/total_direct_time;
		const double indirect_average_gflops = indirect_total_gflop/total_indirect_time;
		
		double num_effective_interactions = (double)evaluator.m_plan.m_direct_buffer_size*(double)evaluator.m_num_target_blocks*(double)_BLOCKSIZE_*(double)_BLOCKSIZE_;
		const double total_wallclock_time = (end-start_before_plan).seconds ();
		const double evaluation_wallclock_time = (end-start_before_run).seconds ();
		const double plan_wallclock_time = (start_before_run-start_before_plan).seconds ();
		const double expansions_wallclock_time = (start_before_plan-start_before_expansions).seconds ();
		//const double fmm_wallclock_time = (end-start_before_tree).seconds ();
		const double tree_wallclock_time = (start_before_expansions-start_before_tree).seconds ();
		
		std::cout << "\n\n///////////////////// SSE Velocity Solver Information /////////////////////\n";
		
		std::cout.precision (4);
		std::cout << "\nGeneral Information\n";
		std::cout << "\n\tNumber of direct interacting source particles: " << evaluator.m_plan.m_direct_buffer_size << "\n";
		std::cout << "\tNumber of target blocks: " << evaluator.m_num_target_blocks << "\n";
		std::cout << "\tEffective number of interactions: " << setprecision (6) << scientific << num_effective_interactions << "\n";
		std::cout << "\tActual number of interactions: " << (total_direct_particles+total_indirect_boxes)*_BLOCKSIZE_*_BLOCKSIZE_ << ", (direct: " << fixed << setprecision (2) << 100*total_direct_particles/(total_direct_particles+total_indirect_boxes) << "%, indirect: " << 100*total_indirect_boxes/(total_direct_particles+total_indirect_boxes) << "%)\n";
		std::cout << "\tNumber of direct c++ kernel calls: " << num_direct_cpp_evals << scientific << " (" << (double)num_direct_cpp_evals/num_direct_kernel_calls << " % of total direct kernel calls ( " << num_direct_kernel_calls << ") )\n";
		std::cout << "\tAverage number of source particles per box: " << average_num_particles_per_box << "\n";
		
		std::cout.precision (2);
		std::cout << "\nMemory\n";
		std::cout << "\n\tTotal memory for direct interactions: " << scientific << (double)3*evaluator.m_plan.m_direct_buffer_size*sizeof(float)/(1024*1024) << " MB\n";
		std::cout << "\tTotal memory for indirect interactions: " << (double)2*evaluator.m_plan.m_indirect_buffer_size*sizeof(float)/(1024*1024) << " MB\n";
		std::cout << "\nElapsed time\n";
		std::cout << "\n\tTotal elapsed wall-clock time: " << total_wallclock_time << "\n";
		std::cout << "\tTotal CPU time: " << total_time << " [s],\t(direct: " << fixed << setprecision (1) << 100*total_direct_time/total_time << "%, indirect: " << 100*total_indirect_time/total_time << "%)\n";
		std::cout.precision (2);
		std::cout << "\tAverage CPU time per target block: " << scientific << total_time/evaluator.m_num_target_blocks << " [s]\n";
		std::cout << "\tAverage direct SSE kernel time per source particle: " <<  total_direct_time/total_direct_particles << " [s]\n";
		std::cout << "\tAverage indirect SSE kernel time per source box: " << total_indirect_time/total_indirect_boxes << " [s]\n";
		std::cout << "\nFloating point operations\n";
		std::cout << "\n\tTotal number of operations: " << total_gflop << " [GFLOP], (direct: " << fixed << setprecision (1) << 100*direct_total_gflop/total_gflop << "%, indirect: " << 100*indirect_total_gflop/total_gflop << "%)\n";
		std::cout.precision (2);
		std::cout << "\tAverage number of operations per target block: " << scientific << total_gflop/evaluator.m_num_target_blocks << " [GFLOP]\n";
		std::cout << "\nPerformance\n";
		std::cout << "\tDirect SSE kernel (per call): " << fixed << setprecision (3) << direct_avg_gflops << " [GFLOP/s]\n";
		std::cout << "\tIndirect SSE kernel (per call): " <<  indirect_avg_gflops << " [GFLOP/s]\n";
		std::cout << "\n\tDirect SSE kernel (average): " << fixed << setprecision (3) << direct_average_gflops << " [GFLOP/s]\n";
		std::cout << "\tIndirect SSE kernel (average): " << indirect_average_gflops << " [GFLOP/s]\n";
		
		std::cout << "\nFloating point operations per second wall-clock time\n";
		std::cout << "\n\tPlan + Evaluations: " << total_gflop/total_wallclock_time << " [GFLOP/s]\n";
		std::cout << "\n\tEvaluations only: " << total_gflop/evaluation_wallclock_time << " [GFLOP/s]\n";
		
		std::fstream file;
		
		//Write data for break even plot
		std::string file_name;
		
		file.open ("measurements.txt",  std::fstream::app | std::fstream::out);
		if (file.fail ()) {
			std::cerr << "\nFailed to open measurements.txt file. Aborting now.\n";
			abort ();
		}
		file << setprecision (6) << scientific << num_effective_interactions << "\t" << total_gflop << "\t" << \
		evaluation_wallclock_time << "\t" << plan_wallclock_time << "\t" << expansions_wallclock_time << "\t" << \
		tree_wallclock_time << "\n";
		file.close ();
		if (file.fail ()) {
			std::cerr << "\nFailed to close measurements.txt file. Aborting now.\n";
			abort ();
		}
		
		std::cout << "\n///////////////////////////////////////////////////////////////////////////\n\n";
		
	}
	
	delete rootBox;
	
	timestamp++;
}
