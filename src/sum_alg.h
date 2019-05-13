//
//  Author: Harry Kai-Ho Chan
//  Email: khchanak@cse.ust.hk
//

#ifndef SUM_ALG_H
#define SUM_ALG_H

#include "unified.h"

//================ sum-exact-new ==========================

obj_set_t* Sum_Exact_new( query_t* q);

//================ sum-exact ==========================

obj_set_t* Sum_Exact( query_t* q);

//================ sum-greedy ==========================

obj_set_t* Sum_Appro( query_t* q);


//====================================================
//=============bit_set for sum exact==================
//====================================================

typedef struct bit_node
{
	BIT_TYPE                  bit_string;
	struct bit_node*     next;
}	bit_node_t;

//The structure for storing a set of sets of keywords.
typedef struct bit_set
{
	int			bit_n;
	bit_node_t*      p_head;//dummy head used
}	bit_set_t;

bit_node_t* alloc_bit_node();

bit_node_t* copy_bit_node(bit_node_t* bit_node_v);

bit_set_t* alloc_bit_set();

bool bit_set_find(bit_set_t* bit_set_v, bit_node_t* bit_v);

void bit_set_insert(bit_set_t* bit_set_v, bit_node_t* bit_v);

bool bit_set_remove(bit_set_t* bit_set_v, bit_node_t* bit_v);

void bit_set_append(bit_set_t* bit_set_v1, bit_set_t* bit_set_v2);

void find_all_bit_set_sub(bit_set_t* bit_set_v, int k, BIT_TYPE bit);

bit_set_t* find_all_bit_set(bit_node_t* bit_node_v);

bit_node_t* psi_to_bit_node(B_KEY_TYPE a[],int size, psi_t* psi_v);

void release_bit_set( bit_set_t* bit_set_v);

void print_bit_node( bit_node_t* bit_node_v);

void print_bit_set( bit_set_t* bit_set_v);

unsigned int map_to_integer(bit_node_t* bit_node_v);

#endif
