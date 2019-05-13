//
//  Author: Harry Kai-Ho Chan
//  Email: khchanak@cse.ust.hk
//


/*
 *  Implementation of the approximate algorithms proposed in the paper
 *  "Efficient Processing of Spatial Group Keyword Queries"
 *  by Xin Cao et al.
 */
#ifndef __coskq__cao_alg_new__
#define __coskq__cao_alg_new__

#include "unified.h"
#include "cao_alg.h"

obj_set_t* Cao_Exact_new( int cost_tag, query_t* q, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap);

void enumrateBestGroup( int cost_tag, obj_set_t* S, obj_t* pivot, B_KEY_TYPE dist, B_KEY_TYPE& curCost, obj_set_t*& curGroup, query_t* q);

//void cao_search( int cost_tag, query_t* q, obj_set_t* selectedSet, obj_set_t* candidateSet, B_KEY_TYPE pairDist, B_KEY_TYPE furDist, int startId, B_KEY_TYPE& curCost, obj_set_t*& curGroup);

void cao_search( int& cost_tag, query_t*& q, obj_set_t*& selectedSet, obj_set_t*& candidateSet, B_KEY_TYPE& pairDist, B_KEY_TYPE& furDist, int& startId, B_KEY_TYPE& curCost, obj_set_t*& curGroup);

obj_set_t* Cao_Appro2_new( query_t* q, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap);

B_KEY_TYPE findMostInfreqKey(psi_t* psi_v, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap);

void add_obj_set_entry_sorted( obj_t* obj_v, obj_set_t* obj_set_v, B_KEY_TYPE dist);

bool is_covered_obj_set( obj_set_t* obj_set_v, psi_t* psi_v);


#endif /* defined(__coskq__cao_alg_new__) */
