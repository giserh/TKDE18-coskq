//
//  Author: Chan Kai Ho
//  Email: khchanak@ust.hk
//

#ifndef GENERAL_H
#define GENERAL_H

#include "cao_alg.h"
#include "cao_alg_new.h"

#define Cost_MaxMax			1
#define Cost_MaxMax2		2
#define Cost_Sum			3
#define Cost_Max			4
#define Cost_MinMax			5
#define Cost_MinMax2		6
#define Cost_SumMax			7
#define Cost_SumMax2		8

//=================================================================================
//=================================================================================
//=======  Approximate Algorithm   ================================================
//=================================================================================
//=================================================================================

obj_set_t* Unified_A( query_t* q, int cost_tag, B_KEY_TYPE cost_c_global);

obj_set_t* Gen_ConstructGreedyFeasibleSet( obj_t* o, query_t* q, int cost_tag, obj_set_t* O_t, B_KEY_TYPE pdist_max);

//=================================================================================
//=================================================================================
//=======  Exact Algorithm  =======================================================
//=================================================================================
//=================================================================================

obj_set_t* Unified_Approach( query_t* q, int s_tag, int cost_tag, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap);

obj_set_t* Unified_E( query_t* q, int cost_tag, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap, B_KEY_TYPE cost_c_global);

void LoopObjPair(int cost_tag, query_t* q, obj_t* obj_v1, obj_t* obj_v2, B_KEY_TYPE d_l, obj_set_t*& S_c, B_KEY_TYPE& cost_c);

bool EnumerateSubset( int cost_tag, obj_set_t* O_t, psi_t* psi,  tri_t* tri_v, query_t* q, obj_set_t*& S_a, B_KEY_TYPE& cost_c, B_KEY_TYPE d);

bool EnumerateSubset_sub(int cost_tag, bst_t* IF_v, obj_set_t* S, obj_t* o, query_t* q, B_KEY_TYPE x, B_KEY_TYPE cost, obj_set_t*& curSet,  B_KEY_TYPE& curCost);

//=================================================================================
//=================================================================================
//=======================Some useful functions=====================================
//=================================================================================
//=================================================================================

void psi_insert( psi_t* psi_v,k_node_t* k_head);

psi_t* key_intersection( k_node_t* k_head1, k_node_t* k_head2);

bool is_contain_key(psi_t* psi_v, B_KEY_TYPE key);

psi_t* psi_exclusion( psi_t* psi_v1, obj_t* obj_v);

int is_relevant_node(node_t* node_v, psi_t* psi_v);

int is_relevant_obj( obj_t* obj_v, psi_t* psi_v);

int is_relevant_node(node_t* node_v, obj_t* obj_v);

psi_t* node_intersection(node_t* node_v, psi_t* psi_v);

int number_intersection( k_node_t* k_head1, k_node_t* k_head2);

psi_t* uncover_keyword( query_t* q, obj_set_t* obj_set_v);

obj_set_t* retrieve_obj_key( KEY_TYPE key);

obj_set_t* range_query_sorted_dist( disk_t* disk_v, query_t* q, loc_t* loc_v);

void range_query_sorted_dist_sub( node_t* node_v, disk_t* disk_v, obj_set_t* &obj_set_v, query_t* q, loc_t* loc_v);

obj_set_t* range_query_sorted_dist( disk_t* disk_v1, disk_t* disk_v2, query_t* q, loc_t* loc_v);

void range_query_sorted_dist_sub( node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t* &obj_set_v, query_t* q, loc_t* loc_v);

void retrieve_sub_tree_sorted_dist( node_t* node_v, obj_set_t* &obj_set_v, query_t* q, loc_t* loc_v);

obj_set_t* range_query_sorted_dist( query_t* q, loc_t* loc_v);

//====

//same as before, no need to overload
//void retrieve_sub_tree( node_t* node_v, obj_set_t* &obj_set_v, query_t* q);


void range_query_sub( node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t* &obj_set_v, query_t* q);

obj_set_t* range_query( disk_t* disk_v1, disk_t* disk_v2, query_t* q);

query_t** gen_query_set4( int query_n, int n, data_t* data_v);

bool is_old2( psi_t* psi_v, B_KEY_TYPE key);

b_heap_t* heap_sort_obj_set( obj_set_t* obj_set_v, loc_t* loc_v2);

b_heap_t* heap_sort_obj_set2( obj_set_t* obj_set_v, loc_t* loc_v2, B_KEY_TYPE dist_max);

b_heap_t* heap_sort_obj_set3( obj_set_t* obj_set_v, query_t* q, B_KEY_TYPE& maxDist);

obj_t* const_NN_key2( loc_t* loc_v, KEY_TYPE key, disk_t* disk_v1, disk_t* disk_v2);

bool check_dist_constraint( obj_set_t* obj_set_v, obj_t* obj_v, B_KEY_TYPE d);

void update_heap_obj(bst_t* heap, bst_node_t* bst_node_v, psi_t* psi_v, query_t* q, bst_node_list_t* bst_node_list_v, psi_t* psi_excluded);

//=======

void retrieve_sub_tree_dominant( node_t* node_v, obj_set_t* &obj_set_v, query_t* q);
void range_query_dominant_sub( node_t* node_v, disk_t* disk_v, obj_set_t* &obj_set_v, query_t* q);
obj_set_t* range_query_dominant( disk_t* disk_v, query_t* q);

void add_obj_set_entry_domanint( obj_t* obj_v, obj_set_t* obj_set_v, query_t* q);
bool check_dominant(node_t* node_v, obj_set_t* obj_set_v, query_t* q);
bool isDominated(node_t* node_v, obj_t* obj_v2, psi_t* psi_v);

obj_set_t* range_query_dominant( disk_t* disk_v1, disk_t* disk_v2, query_t* q);
void range_query_dominant_sub( node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t* &obj_set_v, query_t* q);

bool check_cost(obj_t* obj_v, obj_set_t* S, B_KEY_TYPE pdist_max);

//---

bool isDominated(obj_t* obj_v, obj_t* obj_v2, psi_t* psi_v);

void remove_obj_set_entry( obj_set_t* obj_set_v,  obj_t* obj_v );

void release_IF( bst_t* T);

void release_IF_sub( bst_node_t* x);

#endif




