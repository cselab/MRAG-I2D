/*
 *  MRAGRefiner_BlackList.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/21/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <set>

using namespace std;

#include "MRAGRefiner.h"

namespace MRAG
{
	struct RefinerInspector 
	{
		const int max_jump;
		const int max_level;
		
		const GridNode * requested_node;
		set<int> & black_list;
		
		set<const GridNode *> marked;
		set<const GridNode *> to_refine;
		
		const NeighborhoodType& neighborhood;
		const HierarchyType& hierarchy;
		
		RefinerInspector(const GridNode * initial_node, const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, set<int> & black_list, const int max_jump, const int max_level): 
		requested_node(initial_node), black_list(black_list), hierarchy(hierarchy), neighborhood(neighborhood), max_jump(max_jump), max_level(max_level) {}
		
		void inspect()
		{
			assert(requested_node->level < max_level || max_level == -1);
			
			if (!_step_inspect(requested_node))
				to_refine.clear();
			
			//printf("NODE: %d %d %d %d, SIZE: %d\n", requested_node->index[0], 
			//requested_node->index[1], requested_node->index[2], requested_node->level, to_refine.size());
		}
		
		bool _step_inspect(const GridNode * node)
		{
			const bool bInsideBlackList = ( black_list.find(node->blockID)!=black_list.end());
			if (bInsideBlackList) return false;
			
			marked.insert(node);
			to_refine.insert(node);
			
			NeighborhoodType::const_iterator it = neighborhood.find(const_cast<GridNode*>(node));
			assert(it!=neighborhood.end());
			vector<GridNode *> neighbors = it->second;
			
			for (vector<GridNode *>::iterator it = neighbors.begin(); it!=neighbors.end(); ++it)
			{
				const GridNode * neighbor = _findRealNode(hierarchy, neighborhood, *it);
				
				const bool bMarked = ( marked.find(neighbor)!=marked.end() );
				if (bMarked) continue;
				
				const bool bToRefine = abs(node->level+1 - neighbor->level) > max_jump;
				
				if (bToRefine)
				{
					if (neighbor->level	>= max_level && max_level>=0)
						return false;
					
					if (!_step_inspect(neighbor))
						return false;				
				}
			}
			
			return true;
		}
		
		static GridNode* _findRealNode(const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, GridNode* node)
		{
			//this function finds the "true" grid node when the input grid node is a ghost
			
			NeighborhoodType::const_iterator itNode = neighborhood.find(const_cast<GridNode*>(node));
			const bool bNotGhost = (itNode != neighborhood.end());
			
			if (bNotGhost) 
				return node;
			else 
			{
				//look for the parent's children. the one with the same blockID is the one.				
				HierarchyType::const_iterator itParent = hierarchy.find(node->parent);
				
				assert(itParent != hierarchy.end());
				
				vector<GridNode *> children = itParent->second;
				
				int iChild;
				for(iChild =0; iChild<children.size(); iChild++)
					if (children[iChild]->blockID == node->blockID) break;
				
				assert(children[iChild]->blockID == node->blockID);
				
				return children[iChild];
				
			}
		}
	};
	
	class Refiner_BlackList: public Refiner
	{
		typedef set<const GridNode *> SetOfNodes;
		
		set<int> * black_list;
		
		bool _check(const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, set<int> black_list, SetOfNodes& to_refine)
		{
			for(SetOfNodes::iterator it=to_refine.begin(); it!=to_refine.end(); ++it)
			{		
				assert((*it)->level + 1 <= m_nMaxLevel || m_nMaxLevel==-1);
				assert(black_list.find((*it)->blockID) == black_list.end());
				
				NeighborhoodType::const_iterator itNeighbors = neighborhood.find(const_cast<GridNode*>(*it));
				assert(itNeighbors!=neighborhood.end());
				vector<GridNode *> neighbors = itNeighbors->second;
				
				for (vector<GridNode *>::iterator itNeighbor = neighbors.begin(); itNeighbor!=neighbors.end(); ++itNeighbor)
				{
					const GridNode * neighbor = RefinerInspector::_findRealNode(hierarchy, neighborhood, *itNeighbor);
					
					const bool bToRefine = to_refine.find(neighbor) != to_refine.end();
					
					if (!bToRefine)
						assert(abs(neighbor->level - ((*it)->level+1)) <= m_nMaxLevelJump);
					else
						assert(abs(neighbor->level+1 - ((*it)->level+1)) <= m_nMaxLevelJump);
				}
			}
			
			return true;
		}
		
	public:
		
		Refiner_BlackList(int nMaxLevelJump=1,  int nMaxLevel=-1): Refiner(nMaxLevelJump, nMaxLevel), black_list(NULL)
		{
		}
		
		void set_blacklist(set<int> * black_list)
		{
			this->black_list = black_list;
		}
		
		virtual RefinementPlan* createPlan(const HierarchyType& hierarchy, const NeighborhoodType& neighborhood, const bool vProcessingDirections[3], vector<NodeToRefine>& vRefinements)
		{			
			if(black_list == NULL)
				return Refiner::createPlan(hierarchy, neighborhood, vProcessingDirections, vRefinements);
			else 
			{
				
				SetOfNodes requested_nodes;
				
				for(HierarchyType::const_iterator it=hierarchy.begin(); it!=hierarchy.end(); it++)
				{
					const GridNode * node = it->first;
					
					if (node != NULL && !node->isEmpty && node->shouldBeRefined && (m_nMaxLevel==-1? true: node->level < m_nMaxLevel)) requested_nodes.insert(node);
				}
				
				SetOfNodes final_list;
				
				for(SetOfNodes::iterator it = requested_nodes.begin(); it!=requested_nodes.end(); ++it)
				{
					RefinerInspector inspector(*it, hierarchy, neighborhood, *black_list, this->m_nMaxLevelJump, this->m_nMaxLevel);
					inspector.inspect();
					
					final_list.insert(inspector.to_refine.begin(), inspector.to_refine.end());
				}
				
				//printf("SIZE=%d\n", final_list.size());
				assert(_check(hierarchy, neighborhood, *black_list, final_list));
				
				vRefinements.clear();
				RefinementPlan * plan = new RefinementPlan;
				
				{
					int nNewChildren = 0;
					for(SetOfNodes::const_iterator it = final_list.begin();  it!=final_list.end(); ++it)
					{
						const GridNode * node = *it;
						
						SingleRefinement& singleRefinement = *(plan->createEntry());
						singleRefinement.block_info_source.blockID = node->blockID;
						singleRefinement.block_info_source.index[0] = node->index[0];
						singleRefinement.block_info_source.index[1] = node->index[1];
						singleRefinement.block_info_source.index[2] = node->index[2];
						singleRefinement.block_info_source.level = node->level;
						
						for(int i=0; i<8; i++)
						{
							const int offset[3] = {i&1, (i>>1)&1, (i>>2)&1};
							
							if ((offset[0]==0 || vProcessingDirections[0]) && 
								(offset[1]==0 || vProcessingDirections[1]) &&
								(offset[2]==0 || vProcessingDirections[2]))
							{
								RefinementPlanNode& child = singleRefinement.createEntry();
								child.level = node->level + 1;
								
								child.index[0] = node->index[0]*2 + offset[0];
								child.index[1] = node->index[1]*2 + offset[1];
								child.index[2] = node->index[2]*2 + offset[2];
								
								child.relative_index[0] = offset[0];
								child.relative_index[1] = offset[1];
								child.relative_index[2] = offset[2];
								
								nNewChildren++;
							}
						}
						
						vRefinements.push_back(NodeToRefine(node, vRefinements.size()));
					}
					
					plan->nNewBlocks = nNewChildren;
				}
				
				black_list = NULL;
				
				return plan;
			}

		}
	};
}
