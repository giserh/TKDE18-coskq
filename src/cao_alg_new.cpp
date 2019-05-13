//
//  cao_alg_new.cpp
//  coskq
//
//  Created by Harry on 15/9/15.
//  Copyright (c) 2015 Harry. All rights reserved.
//

#include "cao_alg_new.h"

/*
 *	The implementation of the "Cao-Exact-new" algorithm.
 *
 *  handle MAXMAX and MINMAX cost
 */
obj_set_t* Cao_Exact_new(int cost_tag, query_t* q, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap)
{
    int i, rear_u, top_u, rear_w, top_w, j;
    KEY_TYPE t_inf;
    B_KEY_TYPE curCost, dist;
    b_heap_t *U, *W;
    obj_set_t *curGroup, *S;
    void *e, *e_w;
    node_t *node_v, *node_v_w;
    obj_t *obj_v, *obj_v_w;
    bst_node_t* bst_node_v;
    BIT_TYPE p_list, p_list_w;
    loc_t* loc_v;
    B_KEY_TYPE d_o1_o2;
    B_KEY_TYPE key;
    loc_t* loc_v_temp;

    U = alloc_b_heap(INI_HEAP_SIZE);

    rear_u = 1;
    U->obj_arr[rear_u].element = (void*)IRTree_v.root;
    U->obj_arr[rear_u].e_tag = 1;
    U->obj_arr[rear_u].key = calc_minDist_node(IRTree_v.root, q->loc_v);

    b_h_insert(U, rear_u++);

    if (cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) //MinMax or MinMax2
        curGroup = Cao_Appro1(q);
    else //MaxSum or Dia or SumMax
        curGroup = Cao_Appro2_new(q, keyfreq_hashmap);

    /*s*/
    stat_v.n_k++;
    /*s*/

    if (curGroup == NULL) {
        release_b_heap(U);
        return NULL;
    }

    curCost = comp_cost(cost_tag, curGroup, q);

    //new :find the most infrequent keyword
    t_inf = findMostInfreqKey(q->psi_v, keyfreq_hashmap);

    //Best-first search process.
    while (!b_h_is_empty(U)) {
        top_u = b_h_get_top(U);
        e = U->obj_arr[top_u].element;

        /*t/
         if( e == NULL)
         printf( "");
         /*t*/

        if (U->obj_arr[top_u].e_tag == 1) {
            //e is an node_t*.
            node_v = (node_t*)e;

            //Check the distance.
            dist = calc_minDist_node(node_v, q->loc_v);

            //new: pruning
            if (cost_tag == Cost_MinMax2) {
                if (dist >= 2 * curCost)
                    break;
            } else if (dist >= curCost)
                break;

            bst_node_v = bst_search(node_v->bst_v, t_inf);
            if (bst_node_v == NULL)
                continue;

            p_list = bst_node_v->p_list;
            for (i = 0; i < node_v->num; i++) {
                //Check the c_key keyword.
                if (!get_k_bit(p_list, i))
                    continue;

                /*t/
                 if( node_v->child[ i] == NULL)
                 printf( "");
                 /*t*/

                /*t/
                 if( rear == 9999)
                 printf( "");
                 /*t*/

                U->obj_arr[rear_u].element = node_v->child[i];
                if (node_v->level > 0) {
                    //node_v is an inner-node.
                    U->obj_arr[rear_u].e_tag = 1;
                    U->obj_arr[rear_u].key = calc_minDist(node_v->MBRs[i], q->loc_v);
                } else {
                    //node_v is a leaf-node.
                    U->obj_arr[rear_u].e_tag = 2;
                    U->obj_arr[rear_u].key = calc_minDist(node_v->MBRs[i], q->loc_v);
                }

                //Enqueue.
                b_h_insert(U, rear_u++);
            } //
        } else {

            //e is an obj_t*.
            obj_v = (obj_t*)e;
            //printf("e:%d\n",obj_v->id);
            loc_v = get_obj_loc(obj_v);

            //new: pruning
            dist = calc_dist_loc(q->loc_v, loc_v);
            if (cost_tag == Cost_MinMax2) {
                if (dist >= 2 * curCost)
                    break;
            } else if (dist >= curCost)
                break;

            S = alloc_obj_set();

            //========
            W = alloc_b_heap(INI_HEAP_SIZE);

            rear_w = 1;
            W->obj_arr[rear_w].element = (void*)IRTree_v.root;
            W->obj_arr[rear_w].e_tag = 1;
            W->obj_arr[rear_w].key = 0;

            b_h_insert(W, rear_w++);
            while (!b_h_is_empty(W)) {

                top_w = b_h_get_top(W);
                e_w = W->obj_arr[top_w].element;

                if (cost_tag == Cost_MinMax2) {
                    if (W->obj_arr[top_w].key > 2 * curCost)
                        break;
                } else {
                    if (W->obj_arr[top_w].key > curCost)
                        break;
                }

                if (W->obj_arr[top_w].e_tag == 1) {
                    //e is an node_t*.
                    node_v_w = (node_t*)e_w;
                    p_list_w = is_relevant_node(node_v_w, q);
                    for (j = 0; j < node_v_w->num; j++) {
                        if (!get_k_bit(p_list_w, j))
                            continue;

                        if (node_v_w->level > 0) {
                            node_t* child_node = (node_t*)node_v_w->child[j];

                            d_o1_o2 = calc_minDist_node(child_node, loc_v);
                            key = calc_minDist_node(child_node, q->loc_v);

                            //region checking:

                            if (cost_tag == Cost_MaxMax) //MaxSum
                            {
                                // > E is checked by the key and heap
                                if (d_o1_o2 > curCost - dist) //C_oi^diam
                                    continue;
                                key = fmax(dist, key) + d_o1_o2;

                            } else if (cost_tag == Cost_MaxMax2) //Dia
                            {
                                //no checking
                                key = fmax(fmax(dist, key), d_o1_o2);
                            } else if (cost_tag == Cost_MinMax) //MinMax
                            {
                                if (key > curCost || d_o1_o2 > curCost) //both circle are checked
                                    continue;
                            } else if (cost_tag == Cost_SumMax) //SumMax
                            {
                                if (d_o1_o2 > curCost - dist) //C_oi^diam
                                    continue;
                                key = dist + key + d_o1_o2;

                            } else if (cost_tag == Cost_MinMax2) //MinMax2
                            {
                                if (d_o1_o2 > curCost)
                                    continue;
                            }

                            W->obj_arr[rear_w].element = child_node;
                            W->obj_arr[rear_w].e_tag = 1;
                            W->obj_arr[rear_w].key = key;
                            //Enqueue.
                            b_h_insert(W, rear_w++);

                        } else {
                            obj_t* child_obj = (obj_t*)node_v_w->child[j];

                            loc_v_temp = get_obj_loc(child_obj);
                            d_o1_o2 = calc_dist_loc(loc_v, loc_v_temp);
                            key = calc_dist_loc(loc_v_temp, q->loc_v);
                            release_loc(loc_v_temp);

                            if (cost_tag == Cost_MaxMax) //MaxSum
                            {
                                if (d_o1_o2 > curCost - dist)
                                    continue;
                                key = fmax(dist, key) + d_o1_o2;
                            } else if (cost_tag == Cost_MaxMax2) //Dia
                            {
                                key = fmax(fmax(dist, key), d_o1_o2);
                            } else if (cost_tag == Cost_MinMax) //MinMax
                            {
                                if (key > curCost || d_o1_o2 > curCost)
                                    continue;
                            } else if (cost_tag == Cost_SumMax) //SumMax
                            {
                                if (d_o1_o2 > curCost - dist)
                                    continue;
                                key = dist + key + d_o1_o2;
                            } else if (cost_tag == Cost_MinMax2) // MinMax2
                            {
                                if (d_o1_o2 > curCost)
                                    continue;
                            }

                            W->obj_arr[rear_w].element = child_obj;
                            W->obj_arr[rear_w].e_tag = 2;
                            W->obj_arr[rear_w].key = key;

                            //Enqueue.
                            b_h_insert(W, rear_w++);
                        }
                    }
                } else {
                    //e is an obj_t*.
                    obj_v_w = (obj_t*)e_w;

                    loc_v_temp = get_obj_loc(obj_v_w);
                    B_KEY_TYPE d_o_q = calc_dist_loc(loc_v_temp, q->loc_v);
                    release_loc(loc_v_temp);

                    add_obj_set_entry_sorted(obj_v_w, S, d_o_q);
                }
            }
            release_b_heap(W);
            //            printf("S:\n");
            //            print_obj_set(S, stdout);
            enumrateBestGroup(cost_tag, S, obj_v, dist, curCost, curGroup, q);
            release_obj_set(S);
            //========
            release_loc(loc_v);
        } //else
    } //while

    release_b_heap(U);

    return curGroup;
}

void enumrateBestGroup(int cost_tag, obj_set_t* S, obj_t* pivot, B_KEY_TYPE d_p_q, B_KEY_TYPE& curCost, obj_set_t*& curGroup, query_t* q)
{
    obj_node_t* obj_node_v;
    loc_t* loc_v2;
    B_KEY_TYPE pairDist, furDist;
    obj_set_t *candidateSet, *selectedSet;

    candidateSet = alloc_obj_set();
    obj_node_v = S->head->next;

    //  printf("pivot:%d\t\n",pivot->id);
    //    print_obj_set(S, stdout);

    psi_t* psi_temp = alloc_psi();
    copy_k_list(psi_temp->k_head, pivot->k_head);
    psi_t* psi_v = psi_exclusion(q->psi_v, psi_temp);

    if ((cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) && !is_covered_obj_set(S, psi_v))
        goto E;

    while (obj_node_v != NULL) {

        //  printf("obj_node_v:%d\t key:%f\n",obj_node_v->obj_v->id, obj_node_v->dist);

        //pivot must in selected set
        if (obj_node_v->obj_v == pivot) {
            obj_node_v = obj_node_v->next;
            continue;
        }
        if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_SumMax) //MaxSum or Dia or SumMax
            add_obj_set_entry(obj_node_v->obj_v, candidateSet);

        selectedSet = alloc_obj_set();
        add_obj_set_entry(pivot, selectedSet);
        add_obj_set_entry(obj_node_v->obj_v, selectedSet);
        pairDist = calc_dist_obj(pivot, obj_node_v->obj_v);

        loc_v2 = get_obj_loc(obj_node_v->obj_v);

        if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2) //MaxSum or Dia
            furDist = fmax(d_p_q, calc_dist_loc(loc_v2, q->loc_v));
        else if (cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) //MinMax or MinMax2
            furDist = fmin(d_p_q, calc_dist_loc(loc_v2, q->loc_v));
        else if (cost_tag == Cost_SumMax) //SumMax
            furDist = d_p_q + calc_dist_loc(loc_v2, q->loc_v);

        release_loc(loc_v2);

        if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_SumMax) //MaxSum or Dia or SumMax
        {
            if (is_covered_obj_set(candidateSet, psi_v)) {
                // printf("cao_search...\n");
                int temp = -1;
                cao_search(cost_tag, q, selectedSet, candidateSet, pairDist, furDist, temp, curCost, curGroup);
            }
            //optimal set updating: done in the function
        } else //MinMax or MinMax2
        {
            int temp = -1;
            cao_search(cost_tag, q, selectedSet, S, pairDist, furDist, temp, curCost, curGroup);
        }
        //      printf("end search \t\t:%f\n",stat_v.memory_v / ( 1024 * 1024));

        release_obj_set(selectedSet);
        obj_node_v = obj_node_v->next;
    }
E:
    release_obj_set(candidateSet);
    release_psi(psi_temp);
    release_psi(psi_v);

    //    printf("end enum \t:%f\n",stat_v.memory_v / ( 1024 * 1024));
}

// furDist=
//  max d(o,q) (MaxSum, Dia)
//  min d(o,q) (MinMax, MinMax2)
//  sum d(o,q) (SumMax)

void cao_search(int& cost_tag, query_t*& q, obj_set_t*& selectedSet, obj_set_t*& candidateSet, B_KEY_TYPE& pairDist, B_KEY_TYPE& furDist, int& startId, B_KEY_TYPE& curCost, obj_set_t*& curGroup)
{

    obj_set_t* nextCandSet;
    psi_t *leftKeywords, *psi_v, *psi_v1;
    obj_node_t *o_c, *o_n, *o_s;
    k_node_t* k_node_v;
    B_KEY_TYPE selectedDiam, selectedDiam2, cost, furDist2;
    loc_t* loc_v;

    // printf("\t 1 \t:%f\n",stat_v.memory_v / ( 1024 * 1024));

    if (is_covered_obj_set(selectedSet, q)) {
        //                B_KEY_TYPE cost_temp = comp_cost(selectedSet, q);
        //                printf("cost:%f\t f+p:%f f:%f \t p:%f \n",cost_temp, furDist+pairDist, furDist, pairDist);
        if (cost_tag == Cost_MaxMax || cost_tag == Cost_MinMax) //MaxSum or MinMax
            cost = furDist + pairDist;
        else if (cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax2)
            cost = fmax(furDist, pairDist);
        else if (cost_tag == Cost_SumMax)
            cost = furDist + pairDist;

        if (cost < curCost) {
            release_obj_set(curGroup);
            curGroup = copy_obj_set(selectedSet);
            curCost = cost;
        }
        return;
    }

    // printf("\t 2 \t:%f\n",stat_v.memory_v / ( 1024 * 1024));

    // printf("1\n");
    nextCandSet = alloc_obj_set();
    leftKeywords = alloc_psi();
    //------------------------------------------------------------------------
    psi_v = uncover_keyword(q, selectedSet);
    //for each obj in cand set
    o_c = candidateSet->head->next;
    while (o_c != NULL) {

        //===
        if (cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) { //MinMax
            loc_v = get_obj_loc(o_c->obj_v);
            if (calc_dist_loc(loc_v, q->loc_v) < furDist) {
                release_loc(loc_v);
                o_c = o_c->next;
                continue;
            }
            release_loc(loc_v);
        }
        //===

        //textual pruning
        if (!is_relevant_obj(o_c->obj_v, psi_v)) {
            //        if (number_intersection(psi_v->k_head, o_c->obj_v->k_head)==0){
            o_c = o_c->next;
            continue;
        }

        //avoid duplicate enum.
        if (o_c->obj_v->id < startId) {
            o_c = o_c->next;
            continue;
        }

        //--------------------------------
        //update max d(o1,o2) in the selected set
        //for each obj in the selected set
        //compute the dist between o_c
        selectedDiam = 0.0;
        o_s = selectedSet->head->next;
        while (o_s != NULL) {
            selectedDiam = fmax(selectedDiam, calc_dist_obj(o_s->obj_v, o_c->obj_v));
            o_s = o_s->next;
        }

        if (((cost_tag == Cost_MaxMax || cost_tag == Cost_MinMax || cost_tag == Cost_SumMax) && selectedDiam + furDist > curCost)
            || ((cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax2) && fmax(selectedDiam, furDist) > curCost)) {
            o_c = o_c->next;
            continue;
        }
        //---------------------------------

        add_obj_set_entry(o_c->obj_v, nextCandSet);
        nextCandSet->head->next->dist = selectedDiam;

        //for each keyword in o_c
        k_node_v = o_c->obj_v->k_head->next;
        while (k_node_v != NULL) {
            if (is_contained_query(q, k_node_v->key)) {
                add_psi_entry(leftKeywords, k_node_v->key);
            }
            k_node_v = k_node_v->next;
        }

        o_c = o_c->next;

    } //end while

    //------------------------------------------------------------------------
    //  printf("\t 3 \t:%f\n",stat_v.memory_v / ( 1024 * 1024));

    psi_v1 = psi_exclusion(psi_v, leftKeywords);

    //    print_k_list(psi_v->k_head, stdout);
    //    print_k_list(leftKeywords->k_head, stdout);
    //    print_k_list(psi_v1->k_head, stdout);

    // printf("\t 4 \t:%f\n",stat_v.memory_v / ( 1024 * 1024));

    release_psi(psi_v);
    release_psi(leftKeywords);

    if (psi_v1->key_n != 0) {
        release_psi(psi_v1);
        //        printf("\t 51 \t:%f\n",stat_v.memory_v / ( 1024 * 1024));
        goto E;
    }
    release_psi(psi_v1);

    //    printf("\t 52 \t:%f\n",stat_v.memory_v / ( 1024 * 1024));

    //------------------------------------------------------------------------
    o_n = nextCandSet->head->next;
    while (o_n != NULL) {

        add_obj_set_entry(o_n->obj_v, selectedSet);
        furDist2 = furDist;
        //-------------
        if (cost_tag == Cost_SumMax) { //SumMax
            loc_t* loc_n = get_obj_loc(o_n->obj_v);
            furDist2 = furDist + calc_dist_loc(loc_n, q->loc_v);
            release_loc(loc_n);
            if (furDist2 + fmax(pairDist, o_n->dist) > curCost) {
                remove_obj_set_entry(selectedSet);
                o_n = o_n->next;
                continue;
            }
        }
        //--------------

        if (selectedSet->obj_n <= q->psi_v->key_n) {
            selectedDiam2 = fmax(pairDist, o_n->dist);
            cao_search(cost_tag, q, selectedSet, nextCandSet, selectedDiam2, furDist2, o_n->obj_v->id, curCost, curGroup);
        }
        remove_obj_set_entry(selectedSet);
        o_n = o_n->next;
    }

E:
    release_obj_set(nextCandSet);
    //    printf("\t 6 \t:%f\n",stat_v.memory_v / ( 1024 * 1024));

    return;
}

/*
 *	The implementation of the "Cao-Appro2-new" algorithm.
 *
 */
obj_set_t* Cao_Appro2_new(query_t* q, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap)
{
    int i, rear, top;
    KEY_TYPE c_key;
    B_KEY_TYPE costV, costV_1, dist;
    b_heap_t* U;
    obj_set_t *V, *V_1, *V_tmp;
    void* e;
    node_t* node_v;
    obj_t* obj_v;
    bst_node_t* bst_node_v;
    BIT_TYPE p_list;
    loc_t* loc_v;
    query_t* q_new;

    //printf( "Cao-Appro2:\n");

    U = alloc_b_heap(INI_HEAP_SIZE);

    rear = 1;
    U->obj_arr[rear].element = (void*)IRTree_v.root;
    U->obj_arr[rear].e_tag = 1;
    U->obj_arr[rear].key = calc_minDist_node(IRTree_v.root, q->loc_v);

    b_h_insert(U, rear++);

    V = Cao_Appro1(q);

    /*s*/
    stat_v.n_k++;
    /*s*/

    if (V == NULL) {
        release_b_heap(U);
        return NULL;
    }

    costV = comp_cost(cost_tag, V, q);

    //new :find the most infrequent keyword
    c_key = findMostInfreqKey(q->psi_v, keyfreq_hashmap);

    //Best-first search process.
    while (!b_h_is_empty(U)) {
        top = b_h_get_top(U);
        e = U->obj_arr[top].element;

        /*t/
         if( e == NULL)
         printf( "");
         /*t*/

        if (U->obj_arr[top].e_tag == 1) {
            //e is an node_t*.
            node_v = (node_t*)e;

            //Check the distance.
            dist = calc_minDist_node(node_v, q->loc_v);

            //new: pruning
            if (dist > costV) {
                break;
            }

            bst_node_v = bst_search(node_v->bst_v, c_key);
            if (bst_node_v == NULL)
                continue;

            p_list = bst_node_v->p_list;
            for (i = 0; i < node_v->num; i++) {
                //Check the c_key keyword.
                if (!get_k_bit(p_list, i))
                    continue;

                U->obj_arr[rear].element = node_v->child[i];
                if (node_v->level > 0) {
                    //node_v is an inner-node.
                    U->obj_arr[rear].e_tag = 1;
                    U->obj_arr[rear].key = calc_minDist(node_v->MBRs[i], q->loc_v);
                } else {
                    //node_v is a leaf-node.
                    U->obj_arr[rear].e_tag = 2;
                    U->obj_arr[rear].key = calc_minDist(node_v->MBRs[i], q->loc_v);
                }

                //Enqueue.
                b_h_insert(U, rear++);
            } //
        } else {
            //e is an obj_t*.
            obj_v = (obj_t*)e;
            loc_v = get_obj_loc(obj_v);

            //new: pruning
            dist = calc_dist_loc(q->loc_v, loc_v);
            if (dist > costV) {
                break;
            }

            //Construct a new query instance.
            q_new = alloc_query();
            q_new->loc_v = loc_v;

            q_new->psi_v = alloc_psi();
            copy_k_list(q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

            //Solve the new query.
            V_1 = Cao_Appro1(q_new);
            //add_obj_set_entry( obj_v, V_1);

            if (V_1 == NULL) {
                release_query(q_new);
                continue;
            }

            costV_1 = comp_cost(cost_tag, V_1, q);

            if (costV_1 < costV) {
                costV = costV_1;
                V_tmp = V;
                V = V_1;

                release_obj_set(V_tmp);
            } else
                release_obj_set(V_1);

            release_query(q_new);

            /*s*/
            stat_v.n_k++;
            /*s*/
        } //else
    } //while

    release_b_heap(U);

    return V;
}

//=========================================================================
//=========================================================================
//=========================================================================

B_KEY_TYPE findMostInfreqKey(psi_t* psi_v, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap)
{

    B_KEY_TYPE mostInfreqKey = -1;
    B_KEY_TYPE count = INFINITY;
    k_node_t* k_node_v = psi_v->k_head->next;
    std::unordered_map<KEY_TYPE, KEY_TYPE>::const_iterator got;

    while (k_node_v != NULL) {
        got = keyfreq_hashmap->find(k_node_v->key);

        if (got == keyfreq_hashmap->end())
            exit(-1); //not found means error and exit

        //find min
        if (count > got->second) {
            mostInfreqKey = got->first;
            count = got->second;
        }

        k_node_v = k_node_v->next;
    }

    // printf("key:%f \t cnt:%f\n",mostInfreqKey, count);

    return mostInfreqKey;
}

/*
 *	Add an object entry @obj_v to @obj_set_v.
 *  the linked list is sorted by dist, from smallest to largest
 */
void add_obj_set_entry_sorted(obj_t* obj_v, obj_set_t* obj_set_v, B_KEY_TYPE dist)
{
    obj_node_t *obj_node_v, *obj_node_new;

    obj_node_v = obj_set_v->head;
    while (obj_node_v->next != NULL) {
        if (dist < obj_node_v->next->dist)
            break;
        obj_node_v = obj_node_v->next;
    }

    //obj_node_v pointing to the node before insert position
    //-----
    obj_node_new = (obj_node_t*)malloc(sizeof(obj_node_t));
    memset(obj_node_new, 0, sizeof(obj_node_t));

    /*s*/
    stat_v.memory_v += sizeof(obj_node_t);
    if (stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/

    obj_node_new->obj_v = obj_v;
    obj_node_new->dist = dist;

    obj_node_new->next = obj_node_v->next;
    obj_node_v->next = obj_node_new;

    obj_set_v->obj_n++;
    //-----
}

/*
 *	Check whether the keywords in psi_v are covered by a set of objs in @obj_set_v.
 */
bool is_covered_obj_set(obj_set_t* obj_set_v, psi_t* psi_v)
{
    KEY_TYPE key;
    k_node_t* k_node_iter;
    obj_node_t* obj_node_iter;

    k_node_iter = psi_v->k_head->next;
    while (k_node_iter != NULL) {
        key = k_node_iter->key;

        obj_node_iter = obj_set_v->head->next;
        while (obj_node_iter != NULL) {
            if (has_key_obj(obj_node_iter->obj_v, key))
                break;

            obj_node_iter = obj_node_iter->next;
        }

        if (obj_node_iter == NULL)
            return false;

        k_node_iter = k_node_iter->next;
    }

    return true;
}
