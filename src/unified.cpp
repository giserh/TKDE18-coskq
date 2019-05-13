//
//  general cost
//

#include "unified.h"

//=========  Approximate Algorithm ======================

/*
 * The approximate algorithm "Unified-A"
 *
 */
obj_set_t* Unified_A(query_t* q, int cost_tag, B_KEY_TYPE cost_c_global)
{
    int top;
    B_KEY_TYPE LB, cost_c, cost, dist;
    obj_set_t* S_c; //the current best solution.
    obj_set_t* S; //the newly constructed feasible set.
    obj_set_t* R; //region R.
    disk_t *disk_u, *disk_l, *disk_u2; //the outer and inner disks.
    obj_set_t* O_t;
    obj_t *o = NULL, *o_next, *o_prev;
    loc_t* loc_v;
    b_heap_t* R_heap;

    obj_node_t* obj_node_v;

    B_KEY_TYPE pdist_max = INFINITY;

    //Compute the LB and UB of cost(S*, q).c_1;
    S_c = comp_bounds(cost_tag, q, LB, cost_c);
    ///LB and UB are updated in comp_bounds function

    cost_c = fmin(cost_c, cost_c_global);

    if (S_c == NULL)
        return NULL;

    if (cost_tag == Cost_Max)
        return S_c;

    ///---------------------------------------------------
    //Initialize region R.
    disk_u = alloc_disk(IRTree_v.dim);
    disk_l = alloc_disk(IRTree_v.dim);
    disk_u2 = alloc_disk(IRTree_v.dim);

    ///disk_u & disk_l is centered at loc_v with r=UB & LB

    set_disk(disk_u, q->loc_v, cost_c);
    set_disk(disk_l, q->loc_v, LB);
    set_disk(disk_u2, q->loc_v, cost_c * 2);

    //R storing the objects to be iterated
    // O_t is use to maintain the objects to be picked
    // updated in each iteration
    if (cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) {
        //R storing the objects to be iterated
        R = range_query_sorted_dist(disk_u, q, q->loc_v);
        if (cost_tag == Cost_MinMax)
            //    O_t = range_query_sorted_dist( disk_u, q, q->loc_v);
            O_t = copy_obj_set(R);
        else
            O_t = range_query_sorted_dist(disk_u2, q, q->loc_v);
    } else {
        //Direct range query with the range be a "ring".
        if (cost_tag == Cost_Sum) {
            R = range_query_dominant(disk_u, q);
            O_t = range_query_dominant(disk_l, q);
        } else {
            R = range_query_sorted_dist(disk_u, q, q->loc_v);
            O_t = range_query(disk_l, q);
        }
        ///exclude the objects in region_l that are in the boundary of disk_l
        refine_region(O_t, disk_l);
        ///exclude the objects in region_u that are inside disk_l
        obj_exclusion_disk(R, disk_l);
    }

    ///--------------------------------------------------
    /*s/
     printf( "#cands: %i\n", R->obj_n);
     /*s*/
    /*t/
     print_obj_set( R, stdout);
     /*t*/

    //Pre-checking.
    if (R->obj_n == 0)
        goto E;

    //-while there exist unprocessed relevant object o in R(S)
    obj_node_v = R->head->next;

    while (obj_node_v != NULL) {

        o = obj_node_v->obj_v;

        loc_v = get_obj_loc(o);
        dist = calc_dist_loc(loc_v, q->loc_v);
        release_loc(loc_v);

        ///r_max
        if (dist > cost_c)
            break;

        //   printf("cost_c:%0.3lf\tdist:%0.3lf\n", cost_c, dist);

        if (cost_tag == Cost_MaxMax || cost_tag == Cost_MinMax) //MaxSum or MinMax
            pdist_max = cost_c - dist;
        else if (cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax2) //Dia or MinMax2
            pdist_max = cost_c;

        S = Gen_ConstructGreedyFeasibleSet(o, q, cost_tag, O_t, pdist_max);

        if (cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2)
            remove_obj_set_entry(O_t, o);
        else
            add_obj_set_entry(o, O_t);

        if (!S) {
            obj_node_v = obj_node_v->next;
            continue;
        }
        cost = comp_cost(cost_tag, S, q);

        if (cost <= cost_c) {
            release_obj_set(S_c);
            S_c = S;
            S = NULL;
            cost_c = cost;
        } else {
            release_obj_set(S);
            S = NULL;
        }

        obj_node_v = obj_node_v->next;
    } //while
//Release the resource.
//release_b_heap( R_heap);

E:
    release_obj_set(O_t);

    remove_identical_obj(S_c);
    release_disk(disk_u);
    release_disk(disk_u2);
    release_disk(disk_l);
    release_obj_set(R);
    return S_c;
}

/*
 *	The implementation of the "ConstructGreedyFeasibleSet-Appro" procedure in the paper.
 *
 */
obj_set_t* Gen_ConstructGreedyFeasibleSet(obj_t* o, query_t* q, int cost_tag, obj_set_t* O_t, B_KEY_TYPE pdist_max)
{
    obj_node_t* obj_node_v;
    obj_set_t* S;
    psi_t* psi_v;
    k_node_t* k_node_v;
    loc_t* loc_v;

    bst_t* heap;
    bst_node_t* bst_node_v;

    psi_t* psi_excluded;

    //Obtain the "un-covered" keywords by S.
    k_node_v = key_exclusion(q->psi_v->k_head, o->k_head);
    psi_v = const_psi(k_node_v);

    if (psi_v->key_n == 0) {
        release_psi(psi_v);
        S = alloc_obj_set();
        //Include object o in S.
        add_obj_set_entry(o, S);

        return S;
    }

    S = alloc_obj_set();
    //Include object o in S.
    add_obj_set_entry(o, S);

    loc_v = get_obj_loc(o);

    b_heap_t* R_heap;
    obj_t *o_cur, *o_next, *obj_pick;
    int top;
    if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) //MaxSum or Dia or MinMax or MinMax2
    {
        //Sort the objects in R by their distances to o.
        /// only need to consider objects with d(o,o') < pdist_max
        //        R_heap = heap_sort_obj_set( O_t, loc_v);
        R_heap = heap_sort_obj_set2(O_t, loc_v, pdist_max);
        top = b_h_get_top(R_heap);
        o_next = R_heap->obj_arr[top].obj_v;

    } else { //Sum or SumMax
        heap = bst_ini();

        obj_node_v = O_t->head->next;
        while (obj_node_v != NULL) {

            bst_node_v = (bst_node_t*)malloc(sizeof(bst_node_t));
            memset(bst_node_v, 0, sizeof(bst_node_t));

            /*s*/
            stat_v.memory_v += sizeof(bst_node_t);
            /*s*/

            bst_node_v->obj_v1 = obj_node_v->obj_v;
            bst_node_v->key = calc_minDist(obj_node_v->obj_v->MBR, q->loc_v) / ((double)is_relevant_obj(obj_node_v->obj_v, psi_v));

            bst_insert(heap, bst_node_v);

            obj_node_v = obj_node_v->next;
        }
    }
    //--

    ///while there are keyword not yet covered
    while (psi_v->key_n > 0) {
        //========================================
        //============== PickObject ==============
        //========================================

        obj_pick = NULL;

        if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) //MaxSum or Dia or MinMax
        {
            while (true) {
                //find out next relevant obj
                if (o_next == NULL)
                    break;

                o_cur = o_next;
                top = b_h_get_top(R_heap);
                if (top == 0)
                    o_next = NULL;
                else
                    o_next = R_heap->obj_arr[top].obj_v;

                if (is_relevant_obj(o_cur, psi_v)) {
                    obj_pick = o_cur;
                    break;
                }
            }
        } else if (cost_tag == Cost_Sum || cost_tag == Cost_SumMax) //Sum
        {
            bst_node_v = bst_get_min(heap->root);

            if (bst_node_v != NULL) {
                obj_pick = bst_node_v->obj_v1;
                bst_delete(heap, bst_node_v);
            }
        }
        ///if no object is found, we cannot cover all keywords
        if (obj_pick == NULL) {
            release_obj_set(S);
            S = NULL;
            break;
        }

        //S=S \cup o'
        // \psi = \psi - o'\psi
        add_obj_set_entry(obj_pick, S);
        psi_excluded = psi_exclusion(psi_v, obj_pick);

        //update heap
        if (cost_tag == Cost_Sum || cost_tag == Cost_SumMax) {

            //list storing the nodes to be reinserted
            bst_node_list_t* bst_node_list_v = (bst_node_list_t*)malloc(sizeof(bst_node_list_t));
            memset(bst_node_list_v, 0, sizeof(bst_node_list_t));

            //  printf("1. node_n:%d\n", heap->node_n);
            update_heap_obj(heap, heap->root, psi_v, q, bst_node_list_v, psi_excluded);
            //   printf("2. node_n:%d\n", heap->node_n);

            bst_node_list_t* bst_node_list_iter;
            bst_node_list_iter = bst_node_list_v->next;
            while (bst_node_list_iter != NULL) {
                //        printf("a");
                bst_insert(heap, bst_node_list_iter->bst_node_v);
                bst_node_list_iter = bst_node_list_iter->next;
            }

            //     printf("\n3. node_n:%d\n\n", heap->node_n);

            release_bst_node_list(bst_node_list_v);
        }

        release_psi(psi_excluded);

    } ///end while

    if (psi_v->key_n != 0) {
        release_obj_set(S);
        S = NULL;
    }
    //Release the memory.
    release_loc(loc_v);
    release_psi(psi_v);

    if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) //MaxSum or Dia or MinMax
        release_b_heap(R_heap);
    else {
        bst_release(heap);
    }
    return S;
}

//=================================================================================
//================================   Main Algorithm  ==============================
//=================================================================================
/*
 * The Unfified Approach in the paper.
 *
 * @s_tag
 *  Unified-E = 1
 *  Unified-A = 2
 */
obj_set_t* Unified_Approach(query_t* q, int s_tag, int cost_tag, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap)
{
    obj_set_t* S_c; //the current best solution.
    if (s_tag == 1)
        S_c = Unified_E(q, cost_tag, keyfreq_hashmap, INFINITY);
    else
        S_c = Unified_A(q, cost_tag, INFINITY);
    return S_c;
}

/*
 * The exact algorithm "Unified-E"
 *
 */
obj_set_t* Unified_E(query_t* q, int cost_tag, std::unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap, B_KEY_TYPE cost_c_global)
{
    int rear, top;
    B_KEY_TYPE LB, curCost, cost, d_o1_o2, key, d_o1_q, d_o2_q;
    obj_set_t* curSet; //the current best solution.
    obj_set_t* S; //the newly constructed feasible set.
    obj_set_t* R; //region R.
    disk_t* disk_u; //the outer and inner disks.
    obj_set_t* region_u;
    //obj_set_t* region_l;
    loc_t *loc_v, *loc_v2;
    obj_node_t *obj_node_v, *obj_node_v2, *obj_node_tmp;

    b_heap_t* U;
    obj_t *obj_v, *obj_v2;
    psi_t *psi_v, *psi_v2;

    KEY_TYPE t_inf;
    obj_set_t* obj_set_inf;
    bool t_inf_flag;

    //Compute the LB and UB of cost(S*, q);
    curSet = comp_bounds(cost_tag, q, LB, curCost);
    ///LB and UB are updated in comp_bounds function

    curCost = fmin(curCost, cost_c_global);

    if (curSet == NULL)
        return NULL;
    ///S_a is NN
    /*t/
     printf("LB:%0.5lf\n",LB);
     printf("UB:%0.5lf\n",UB);
     /*t*/
    if (cost_tag == Cost_Max)
        return curSet;
    ///---------------------------------------------------
    //Initialize region R.
    //An alternative implementation is possible here.
    //Direct range query with the range be a "ring".
    disk_u = alloc_disk(IRTree_v.dim);

    ///disk_u & disk_l is centered at loc_v with r=UB & LB
    if (cost_tag == Cost_MinMax2)
        set_disk(disk_u, q->loc_v, 2 * curCost);
    else
        set_disk(disk_u, q->loc_v, curCost);

    ///region_u & region_l store the set of objects inside disk_u & disk_l
    if (cost_tag == Cost_Sum)
        region_u = range_query_dominant(disk_u, q);
    else
        region_u = range_query(disk_u, q);

    R = region_u;

    //Pre-checking.
    if (R->obj_n == 0)
        goto E;

    U = alloc_b_heap(INI_HEAP_SIZE);
    rear = 1;

    //-----------------------------------------------------------
    // inintialize the obj obj pair

    //---
    //find the most infrequent keyword
    t_inf = findMostInfreqKey(q->psi_v, keyfreq_hashmap);
    obj_set_inf = alloc_obj_set();
    obj_node_tmp = R->head->next;
    while (obj_node_tmp != NULL) {
        if (has_key_obj(obj_node_tmp->obj_v, t_inf))
            add_obj_set_entry(obj_node_tmp->obj_v, obj_set_inf);
        obj_node_tmp = obj_node_tmp->next;
    }
    //printf("t_inf:%f\t obj_n:%d\n",t_inf, obj_set_inf->obj_n);
    //---

    //printf("|R|:%d\n", R->obj_n);

    //    printf("query:%d\n", q->psi_v->key_n);
    //    print_k_list(q->psi_v->k_head, stdout);

    obj_node_v = R->head->next;
    while (obj_node_v != NULL) {
        obj_node_v->psi_v = key_intersection(q->psi_v->k_head, obj_node_v->obj_v->k_head);

        obj_node_v = obj_node_v->next;
    }

    //   printf("=====\n");
    obj_node_v = R->head->next;
    while (obj_node_v->next != NULL) {

        loc_v = get_obj_loc(obj_node_v->obj_v);
        d_o1_q = calc_dist_loc(loc_v, q->loc_v);

        obj_node_v2 = obj_node_v->next;
        while (obj_node_v2 != NULL) {

            //--pair pruning---
            //if o1 cover all keywords that o2 cover, no need to have this pair
            //similarly, if o2 cover all keywords that o1 cover, no need to have this pair
            int cnt = number_intersection(obj_node_v->psi_v->k_head, obj_node_v2->psi_v->k_head);

            //            print_k_list(obj_node_v->psi_v->k_head, stdout);
            //            print_k_list(obj_node_v2->psi_v->k_head, stdout);
            //            printf("cnt:%d\n",cnt);

            if (cnt == obj_node_v->psi_v->key_n || cnt == obj_node_v2->psi_v->key_n) {
                obj_node_v2 = obj_node_v2->next;
                continue;
            }
            //--

            loc_v2 = get_obj_loc(obj_node_v2->obj_v);
            d_o2_q = calc_dist_loc(loc_v2, q->loc_v);
            d_o1_o2 = calc_dist_loc(loc_v, loc_v2);

            //check whether the area within disk D(q,d_l)
            //if yes, the area does not contain feasible set
            if (fmin(d_o1_q, d_o2_q) + d_o1_o2 < LB) {

                release_loc(loc_v2);
                obj_node_v2 = obj_node_v2->next;
                continue;
            }
            //--

            //textual pruning: check whether the area contain t_inf
            t_inf_flag = false;
            obj_node_tmp = obj_set_inf->head->next;
            while (obj_node_tmp != NULL) {
                if (calc_dist_obj(obj_node_v->obj_v, obj_node_tmp->obj_v) <= d_o1_o2 && calc_dist_obj(obj_node_v2->obj_v, obj_node_tmp->obj_v) <= d_o1_o2) {
                    t_inf_flag = true;
                    break;
                }
                obj_node_tmp = obj_node_tmp->next;
            }

            if (!t_inf_flag) {
                release_loc(loc_v2);
                obj_node_v2 = obj_node_v2->next;
                continue;
            }
            //--

            //key = lower bound of the cost of the set contain o1 o2
            if (cost_tag == Cost_MaxMax) //MaxSum
                key = d_o1_o2 + fmax(fmax(d_o1_q, d_o2_q), LB);
            else if (cost_tag == Cost_MaxMax2) //Dia
                key = fmax(d_o1_o2, fmax(fmax(d_o1_q, d_o2_q), LB));
            else if (cost_tag == Cost_Sum) //Sum
                key = d_o1_q + d_o2_q;
            else if (cost_tag == Cost_MinMax) //MinMax
                key = fmax(fmax(d_o1_q, d_o2_q), fmax(LB, d_o1_o2));
            else if (cost_tag == Cost_SumMax) //SumMax
                key = d_o1_o2 + d_o1_q + d_o2_q;
            else if (cost_tag == Cost_MinMax2) //MinMax2
                key = fmax(d_o1_o2, fmax(d_o1_q, d_o2_q) - d_o1_o2);
            /*
             else if (cost_tag==12)//SumMax2
             key = fmax( d_o1_q + d_o2_q, d_o1_o2);
             */
            else
                key = d_o1_o2;

            if ((cost_tag == Cost_MaxMax && d_o1_o2 < curCost - LB && d_o1_o2 > LB - fmin(d_o1_q, d_o2_q) && key < curCost)
                || (cost_tag == Cost_MaxMax2 && d_o1_o2 < curCost && d_o1_o2 >= LB - fmin(d_o1_q, d_o2_q) && key < curCost)
                || (cost_tag == Cost_Sum && d_o1_o2 < curCost - LB && d_o1_o2 > LB - fmin(d_o1_q, d_o2_q) && key < curCost)
                || (cost_tag == Cost_MinMax && d_o1_o2 < curCost && d_o1_o2 > LB - fmin(d_o1_q, d_o2_q) && key < curCost)
                || (cost_tag == Cost_SumMax && 2 * d_o1_o2 < curCost && d_o1_o2 < curCost - LB && d_o1_o2 > LB - fmin(d_o1_q, d_o2_q) && key < curCost)
                || (cost_tag == Cost_MinMax2 && d_o1_o2 < curCost && d_o1_o2 > LB - fmin(d_o1_q, d_o2_q) && key < curCost)
                //||(cost_tag==12 && d_o1_o2 < cost_c && key < cost_c)
                ) {

                U->obj_arr[rear].element = obj_node_v->obj_v;
                U->obj_arr[rear].element2 = obj_node_v2->obj_v;
                U->obj_arr[rear].key = key;

                b_h_insert(U, rear++);
            }

            release_loc(loc_v2);
            obj_node_v2 = obj_node_v2->next;
        }
        release_loc(loc_v);

        obj_node_v = obj_node_v->next;
    }

    //  printf("qwqw\n");

    //-----------------------------------------------------------
    //for each pair (o1, e2)
    //note that we enforce o1 o2 cannot be the same object
    while (!b_h_is_empty(U)) {
        //  printf("==========\n");
        top = b_h_get_top(U);
        obj_v = (obj_t*)U->obj_arr[top].element;
        obj_v2 = (obj_t*)U->obj_arr[top].element2;

        //cnt++;
        //printf("top:%d\te:%d\te2:%d \t cnt:%d\n",top, obj_v->id ,obj_v2->id, cnt);
        //checking UB

        if (U->obj_arr[top].key >= curCost)
            break;

        //Pre-checking (for the boundary case that |S| < 3).
        S = alloc_obj_set();
        add_obj_set_entry(obj_v, S);
        add_obj_set_entry(obj_v2, S);

        //---
        if (is_covered_obj_set(S, q)) {
            cost = comp_cost(cost_tag, S, q);
            if (cost < curCost) {
                release_obj_set(curSet);
                curSet = S;
                curCost = cost;
            } else
                release_obj_set(S);
        } else {
            release_obj_set(S);

            LoopObjPair(cost_tag, q, obj_v, obj_v2, LB, curSet, curCost);
        }

    } //while

    release_b_heap(U);
// printf("cnt:%d\n",cnt);
E:
    remove_identical_obj(curSet);
    release_disk(disk_u);

    obj_node_v = R->head->next;
    while (obj_node_v != NULL) {
        release_psi(obj_node_v->psi_v);
        obj_node_v = obj_node_v->next;
    }
    release_obj_set(R);
    //release_obj_set( region_l);
    return curSet;
}

/*
 * The "LoopObjPair" algorithm
 * S_c and cost_c store the current best set and cost, respectively.
 *
 */
void LoopObjPair(int cost_tag, query_t* q, obj_t* obj_v1, obj_t* obj_v2, B_KEY_TYPE d_l, obj_set_t*& curSet, B_KEY_TYPE& curCost)
{
    tri_t tri_v;
    bool feasibleSetFlag = false;
    obj_set_t *S, *O_t, *R;
    loc_t *loc_v, *loc_v1, *loc_v2;
    disk_t *disk_v1, *disk_v2;
    B_KEY_TYPE dist, d_o1_q, d_o2_q, d_o1_o2, cost, r_min, r_max, r_max2;
    b_heap_t* R_heap;
    obj_t *o, *o_next;
    int top;
    psi_t *psi_v, *psi_v2, *psi_v3;
    k_node_t *k_head_1, *k_head_2, *k_head_3;
    obj_t* o_q;

    loc_v1 = get_obj_loc(obj_v1);
    loc_v2 = get_obj_loc(obj_v2);

    d_o1_q = calc_dist_loc(loc_v1, q->loc_v);
    d_o2_q = calc_dist_loc(loc_v2, q->loc_v);

    d_o1_o2 = calc_dist_loc(loc_v1, loc_v2);

    //------------------------------------------------

    //pruning: r_min and r_max for d(o,q)
    if (cost_tag == Cost_MaxMax) { //MaxSum
        r_min = fmax(d_l, fmax(d_o1_q, d_o2_q));
        r_max = curCost - d_o1_o2;
    } else if (cost_tag == Cost_MaxMax2) { //Dia
        r_min = fmax(d_o1_o2, d_l);
        r_max = curCost;
    } else if (cost_tag == Cost_Sum) { //Sum
        r_min = fmax(d_l, fmax(d_o1_q, d_o2_q));
        r_max = curCost - d_o1_q - d_o2_q;
    } else if (cost_tag == Cost_MinMax) { //Min+Max
        r_min = fmax(d_l - d_o1_o2, 0);
        r_max = fmin(curCost - d_o1_o2, fmin(d_o1_q, d_o2_q));
    } else if (cost_tag == Cost_SumMax) { //Sum+Max
        r_min = fmax(d_l, fmax(d_o1_q, d_o2_q));
        r_max = curCost - d_o1_o2;
    } else if (cost_tag == Cost_MinMax2) { //MinMax2
        r_min = fmax(d_l - d_o1_o2, 0);
        r_max = fmin(d_o1_q, d_o2_q);
        /*
         }else if (cost_tag==12){    //SumMax2
         r_min = fmax(d_l, fmax(d_o1_q,d_o2_q));
         r_max = cost_c;*/
    } else {
        r_min = 0;
        r_max = INFINITY;
    }
    //------------------------------------------------
    disk_v1 = const_disk(loc_v1, d_o1_o2);
    disk_v2 = const_disk(loc_v2, d_o1_o2);

    if (cost_tag == Cost_Sum)
        R = range_query_dominant(disk_v1, disk_v2, q);
    else
        R = range_query(disk_v1, disk_v2, q);

    release_loc(loc_v1);
    release_loc(loc_v2);

    k_head_1 = key_exclusion(q->psi_v->k_head, obj_v1->k_head);
    k_head_2 = key_exclusion(k_head_1, obj_v2->k_head);
    psi_v = const_psi(k_head_2);
    release_k_list(k_head_1);

    //        printf("%d\t%d\n",psi_v->key_n,psi_v2->key_n);
    //        print_k_list(psi_v->k_head, stdout);
    //        print_k_list(psi_v2->k_head, stdout);

    //O_t is the object set to be enumerated
    //online update: insert o into O_t after each iteration
    //For MinMax and MinMax2: remove o from O_t after each iteration
    //-------
    //psi_v2 is used to keep check whether O_t can cover all keywords
    //online update: when an object is inserted to O_t, psi_v2 is updated
    //enumerate only after psi_v2 can cover all keywords
    //For MinMax and MinMax2: freq-- for the keywords contained in the removed object
    //if freq of a keyword in psi_v2 ==0 --> cannot cover all keywords
    if (cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) { //MinMax or MinMax2
        //        O_t = copy_obj_set(R);
        O_t = R;
        psi_v2 = alloc_psi();
        copy_k_list(psi_v2->k_head, psi_v->k_head);
        psi_v2->key_n = psi_v->key_n;

        //--update freq in each node in psi_v2
        k_node_t* k_node_2 = psi_v2->k_head->next;
        while (k_node_2 != NULL) {
            k_node_2->freq = 0;
            obj_node_t* obj_node_v = O_t->head->next;
            while (obj_node_v != NULL) {
                if (has_key_obj(obj_node_v->obj_v, k_node_2->key))
                    k_node_2->freq++;
                obj_node_v = obj_node_v->next;
            }
            k_node_2 = k_node_2->next;
        }

        //----------------------------------------
    } else {
        O_t = alloc_obj_set();
        psi_v2 = alloc_psi();
        copy_k_list(psi_v2->k_head, psi_v->k_head);
        psi_v2->key_n = psi_v->key_n;
    }

    //Check whether R covers the keywords in the query.
    if (!FeasibilityCheck(R, psi_v))
        goto E;

    //------------------------------------------------
    //pre-checking: Dia
    if (cost_tag == Cost_MaxMax2) {

        //o1 o2 are the distance owner group
        if ((d_o1_o2 >= fmax(d_o1_q, d_o2_q))) {
            //Create o_q.
            o_q = (obj_t*)malloc(sizeof(obj_t));
            alloc_obj(o_q, q->loc_v->dim);
            for (int j = 0; j < q->loc_v->dim; j++) {
                o_q->MBR[j].min = q->loc_v->coord[j];
                o_q->MBR[j].max = o_q->MBR[j].min;
            }
            tri_v.o = o_q;
            tri_v.o_1 = obj_v1;
            tri_v.o_2 = obj_v2;

            feasibleSetFlag = EnumerateSubset(cost_tag, R, psi_v, &tri_v, q, curSet, curCost, d_o1_o2);

            if (feasibleSetFlag) {
                goto E;
            }
        }
        //* Dia: in the following, o, o_q are the distance owner group
    }

    //------------------------------------------------
    //Sort the objects in R by their distances to q.
    R_heap = heap_sort_obj_set(R, q);
    top = b_h_get_top(R_heap);
    o_next = R_heap->obj_arr[top].obj_v;

    tri_v.o_1 = obj_v1;
    tri_v.o_2 = obj_v2;
    while (true) {

        if (o_next == NULL)
            break;

        o = o_next;
        /*t/
         printf("\t o:%d\n",o->id);
         /*t*/
        top = b_h_get_top(R_heap);
        if (top == 0)
            o_next = NULL;
        else
            o_next = R_heap->obj_arr[top].obj_v;

        loc_v = get_obj_loc(o);
        dist = calc_dist_loc(loc_v, q->loc_v);
        release_loc(loc_v);

        //r_min
        if (dist < r_min) {
            if (!(cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2)) {
                add_obj_set_entry(o, O_t);
                if (psi_v2->key_n > 0)
                    psi_exclusion(psi_v2, o);
            }
            continue;
        }

        //r_max
        if (dist > r_max)
            break;
        //*/

        tri_v.o = o;

        /*t/
         printf("o,o_1, o_2: %d,%d,%d\n",tri_v.o->id,tri_v.o_1->id,tri_v.o_2->id);
         /*t*/

        //psi_v3 to store q.\psi - (o1.\psi + o2.psi + o.\psi) for enum
        k_head_3 = key_exclusion(k_head_2, tri_v.o->k_head);
        psi_v3 = const_psi(k_head_3);

        //3 objects cover all query keywords
        if (psi_v3->key_n == 0) {
            S = const_obj_set(&tri_v);
            cost = comp_cost(cost_tag, S, q);

            if (cost <= curCost) {
                release_obj_set(curSet);
                curSet = S;
                S = NULL;
                curCost = cost;
            } else
                release_obj_set(S);

            release_psi(psi_v3);
            break;
        }

        //--
        if (cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) {
            feasibleSetFlag = EnumerateSubset(cost_tag, O_t, psi_v3, &tri_v, q, curSet, curCost, d_o1_o2);
            remove_obj_set_entry(O_t, obj_v1);

            k_node_t* k_node_2 = psi_v2->k_head->next;
            while (k_node_2 != NULL) {
                if (has_key_obj(obj_v1, k_node_2->key)) {
                    k_node_2->freq--;
                    if (k_node_2->freq == 0) {
                        release_psi(psi_v3);
                        goto E2;
                    }
                }

                k_node_2 = k_node_2->next;
            }

        } else {
            if (psi_v2->key_n > 0)
                psi_exclusion(psi_v2, o);

            if (psi_v2->key_n == 0) //if O_t is a feasible set
            {
                // printf("dist:%f\td_o1_o2:%f\n",dist,d_o1_o2);
                //start enumeration
                if (cost_tag == Cost_MaxMax2)
                    feasibleSetFlag = EnumerateSubset(cost_tag, O_t, psi_v3, &tri_v, q, curSet, curCost, dist);
                else
                    feasibleSetFlag = EnumerateSubset(cost_tag, O_t, psi_v3, &tri_v, q, curSet, curCost, d_o1_o2);
            }
            add_obj_set_entry(o, O_t);
        }

        release_psi(psi_v3);

        if (feasibleSetFlag)
            break;

    } ///end while

E2:
    release_b_heap(R_heap);

E:
    release_disk(disk_v1);
    release_disk(disk_v2);
    if (!(cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2))
        release_obj_set(O_t);
    release_obj_set(R);
    release_psi(psi_v);
    release_psi(psi_v2);
    return;
}

/*
 * "EnumerateSubset".
 * S_c and cost_c store the current best set and cost, respectively.
 //Sum or SumMax must reutrn false
 */
bool EnumerateSubset(int cost_tag, obj_set_t* O_t, psi_t* psi, tri_t* tri_v, query_t* q, obj_set_t*& curSet, B_KEY_TYPE& curCost, B_KEY_TYPE d)
{

    bst_t* IF_v;
    obj_set_t* S;
    bool feasibleSetFlag = false;
    B_KEY_TYPE cost;

    //Construct an IF structure.
    IF_v = const_IF(O_t, psi);

    //Initialize the S_0.
    S = alloc_obj_set();
    add_obj_set_entry(tri_v->o, S);

    if (tri_v->o != tri_v->o_1)
        add_obj_set_entry(tri_v->o_1, S);
    if (tri_v->o != tri_v->o_2)
        add_obj_set_entry(tri_v->o_2, S);

    cost = comp_cost(cost_tag, S, q);

    //printf("EnumerateSubset_sub\n");
    //Invoke the sub-procedure "recursively".
    feasibleSetFlag = EnumerateSubset_sub(cost_tag, IF_v, S, tri_v->o, q, d, cost, curSet, curCost);
    //Release the resources.
    //bst_release( IF_v);

    release_IF(IF_v);
    release_obj_set(S);

    return feasibleSetFlag;
}

/*
 *	The sub-procedure of function "EnumerateSubset".
 *  curSet and curCost store the current best set and cost, respectively.
 *  x = d(o1,o2)
 *  cost = cost of the current set (updateable by cost_new when Sum, SumMax)
 */
bool EnumerateSubset_sub(int cost_tag, bst_t* IF_v, obj_set_t* S, obj_t* o, query_t* q, B_KEY_TYPE x, B_KEY_TYPE cost, obj_set_t*& curSet, B_KEY_TYPE& curCost)
{
    obj_t* obj_v;
    bst_node_t* bst_node_v;
    obj_node_t* obj_node_v;
    bst_node_list_t* bst_node_list_v;
    bool tag = false;
    loc_t* loc_v;
    B_KEY_TYPE cost_new = cost;

    if (IF_v->node_n == 0) {
        return true; // T/F does not matter here
    }

    bst_node_v = IF_v->root;
    obj_node_v = bst_node_v->p_list_obj->head->next;

    while (obj_node_v != NULL) { //Pick an object.
        obj_v = obj_node_v->obj_v;

        //Distance constraint checking.
        //MaxSum or Dia or MinMax or SumMax or MinMax2
        if ((cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax || cost_tag == Cost_SumMax || cost_tag == Cost_MinMax2) && !check_dist_constraint(S, obj_v, o, x)) {
            obj_node_v = obj_node_v->next;
            continue;
        }

        if (cost_tag == Cost_Sum || cost_tag == Cost_SumMax) {
            loc_v = get_obj_loc(obj_v);
            cost_new = cost + calc_dist_loc(loc_v, q->loc_v);
            release_loc(loc_v);
            if (cost_new > curCost) {
                obj_node_v = obj_node_v->next;
                continue;
            }
        }

        //Update the IF_v.
        bst_node_list_v = update_IF_obj(IF_v, obj_v);

        //Update the S_0.
        //obj_v is added at the first place of S_0.
        add_obj_set_entry(obj_v, S);

        //Recursive call
        if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) //MaxSum or Dia or MinMax or MinMax2
        {
            tag = EnumerateSubset_sub(cost_tag, IF_v, S, o, q, x, cost, curSet, curCost);
        } else if (cost_tag == Cost_Sum || cost_tag == Cost_SumMax) //Sum or SumMax
        {
            EnumerateSubset_sub(cost_tag, IF_v, S, o, q, x, cost_new, curSet, curCost);
        }

        //Checking.
        if (IF_v->node_n == 0) //== if(FeasibilityCheck(S_0, psi)){
        {
            //update current best
            if (cost_new < curCost) {
                release_obj_set(curSet);
                curSet = copy_obj_set(S);
                curCost = cost_new;
            }
            if (cost_tag == Cost_MaxMax || cost_tag == Cost_MaxMax2 || cost_tag == Cost_MinMax || cost_tag == Cost_MinMax2) //MaxSum or Dia or MinMax or MinMax2
                tag = true; //OK to return
        }

        //Restore the S_0.
        remove_obj_set_entry(S);

        //Restore the IF_v.
        restore_IF_bst_node_list(IF_v, bst_node_list_v);
        release_bst_node_list(bst_node_list_v);

        if (tag)
            return true;

        //Try the next object candidate.
        obj_node_v = obj_node_v->next;
    }
    return false;
}
//=====================================================
//================Some useful functions================
//=====================================================

/*
 *	Exclude the keywords that occur in @obj_v from @psi_v1.
 *
 *	Return the resulting keywords in @psi_v1.
 
 * modified from psi_exclusion
 * 20160714 added - return the keyword(s) excluded
 */
psi_t* psi_exclusion(psi_t* psi_v1, obj_t* obj_v)
{

    psi_t* psi_excluded = alloc_psi();
    k_node_t* k_node_excluded = psi_excluded->k_head; //pointer always point to the end of link list

    int tag;
    k_node_t *k_node_prev, *k_node_v1, *k_node_v2;

    k_node_prev = psi_v1->k_head;
    k_node_v1 = psi_v1->k_head->next;
    while (k_node_v1 != NULL) {
        tag = 0;
        k_node_v2 = obj_v->k_head->next;
        while (k_node_v2 != NULL) {
            if (k_node_v2->key == k_node_v1->key) {
                tag = 1;
                break;
            }

            k_node_v2 = k_node_v2->next;
        }
        if (tag == 1) {
            //The current keyword should be deleted from psi_v1

            //--
            add_keyword_entry(k_node_excluded, k_node_v1->key); //pointed moved to new end inside the function
            psi_excluded->key_n++;
            //--

            k_node_prev->next = k_node_v1->next;
            k_node_v1->next = NULL;
            free(k_node_v1);
            psi_v1->key_n--;

            k_node_v1 = k_node_prev->next;
        } else {
            k_node_v1 = k_node_v1->next;
            k_node_prev = k_node_prev->next;
        }
    }

    return psi_excluded;
}

/*
 *	Construct a psi_t structure based on a list of keywords @k_head.
 *  Deep copy performed.
 */
void psi_insert(psi_t* psi_v, k_node_t* k_head)
{
    k_node_t* k_node_v;

    k_node_v = k_head->next;
    while (k_node_v != NULL) {
        add_keyword_entry(psi_v->k_head, k_node_v->key);
        psi_v->key_n++;

        k_node_v = k_node_v->next;
    }

    return;
}

/******
 *	Find the keywords that occur in @k_head2 and @k_head1.
 *  modified from "number_intersection"
 ******/
psi_t* key_intersection(k_node_t* k_head1, k_node_t* k_head2)
{
    k_node_t *k_node_v1, *k_node_v2;
    psi_t* psi_v;

    psi_v = alloc_psi();
    k_node_t* k_node_v = psi_v->k_head;

    k_node_v1 = k_head1->next;
    while (k_node_v1 != NULL) {
        //printf("k_node_v1->key: %.0lf \n",k_node_v1->key);

        k_node_v2 = k_head2->next;
        while (k_node_v2 != NULL) {
            //  printf("k_node_v2->key: %.0lf \n",k_node_v2->key);

            if (k_node_v2->key == k_node_v1->key) {
				add_keyword_entry(k_node_v, k_node_v2->key);
                psi_v->key_n++;
				break;
            }
            k_node_v2 = k_node_v2->next;
        }
        k_node_v1 = k_node_v1->next;
    }

    return psi_v;
}

//==================================
bool is_contain_key(psi_t* psi_v, B_KEY_TYPE key){
	k_node_t* k_node_v = psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		if( k_node_v->key == key)
			return true;
		k_node_v = k_node_v->next;
	}
	return false;
	//---
}


/*
 *	Check whether a node @node_v is "relevant" or not.
 *  modified from is_relevant_node( node_v, q)
 */
int is_relevant_node(node_t* node_v, psi_t* psi_v)
{
    int cnt = 0;
    k_node_t* k_node_v;

    k_node_v = psi_v->k_head->next;
    while (k_node_v != NULL) {
        if (has_key_node(node_v, k_node_v->key))
            cnt++;

        k_node_v = k_node_v->next;
    }

    return cnt;
}

/*
 *	Check whether an object @obj_v is "relevant" to the query @q.
 *	That is, whether @obj_v contains a keyword in the query @q.
 *  modified from is_relevant_obj( obj_v, q)
 */
int is_relevant_obj(obj_t* obj_v, psi_t* psi_v)
{
    int cnt = 0;
    k_node_t* k_node_v;

    k_node_v = psi_v->k_head->next;
    while (k_node_v != NULL) {
        if (has_key_obj(obj_v, k_node_v->key))
            cnt++;

        k_node_v = k_node_v->next;
    }

    return cnt;
}

/*
 *	Check whether a node @node_v is "relevant" or not.
 *  modified from is_relevant_node( node_v, q)
 */
int is_relevant_node(node_t* node_v, obj_t* obj_v)
{
    int cnt = 0;
    k_node_t* k_node_v;

    k_node_v = obj_v->k_head->next;
    while (k_node_v != NULL) {
        if (has_key_node(node_v, k_node_v->key))
            cnt++;

        k_node_v = k_node_v->next;
    }

    return cnt;
}

/*
 *	return the list of keywords exist in both node_v and psi_v
 *  modified from is_relevant_node( node_v, q)
 */
psi_t* node_intersection(node_t* node_v, psi_t* psi_v)
{
    psi_t* psi_rtn;
    k_node_t* k_node_v;
    B_KEY_TYPE key;

    psi_rtn = alloc_psi();
    k_node_v = psi_v->k_head->next;
    while (k_node_v != NULL) {
        key = k_node_v->key;
        if (has_key_node(node_v, key))
            add_psi_entry(psi_rtn, key);

        k_node_v = k_node_v->next;
    }

    return psi_rtn;
}

/******
 *	Count the number of keywords that occur in @k_head2 and @k_head1.
 *  modified from "key_exclusion"
 ******/
/* note that there MUST be NO dulipcate in the list
 */
int number_intersection(k_node_t* k_head1, k_node_t* k_head2)
{
    int count = 0;
    k_node_t *k_node_v1, *k_node_v2;
    KEY_TYPE key;

    k_node_v1 = k_head1->next;
    while (k_node_v1 != NULL) {
        key = k_node_v1->key;
        k_node_v2 = k_head2->next;
        while (k_node_v2 != NULL) {
            if (k_node_v2->key == key) {
                count++;
                break;
            }
            k_node_v2 = k_node_v2->next;
        }
        k_node_v1 = k_node_v1->next;
    }

    return count;
}

//remove the keywords in the object in the object set obj_set_v
psi_t* uncover_keyword(query_t* q, obj_set_t* obj_set_v)
{

    psi_t* psi_v;
    obj_node_t* obj_node_v;
    k_node_t *k_node_v, *k_node_v2;
    bool flag;

    psi_v = alloc_psi();

    //for each query keyword
    k_node_v = q->psi_v->k_head->next;
    while (k_node_v != NULL) {

        flag = false;

        //for each object
        obj_node_v = obj_set_v->head->next;
        while (obj_node_v != NULL) {

            //for each keyword in the object
            k_node_v2 = obj_node_v->obj_v->k_head->next;
            while (k_node_v2 != NULL) {

                if (k_node_v->key == k_node_v2->key) {
                    flag = true;
                    break;
                }
                k_node_v2 = k_node_v2->next;
            }
            if (flag)
                break;
            obj_node_v = obj_node_v->next;
        }
        //if no object contain this keyword
        if (!flag) {
            add_psi_entry(psi_v, k_node_v->key);
        }

        k_node_v = k_node_v->next;
    }

    return psi_v;
}

/*
 *	Retrieve all the objects that contain keyword @key.
 *  IR-tree is used
 */
obj_set_t* retrieve_obj_key(KEY_TYPE key)
{
    obj_set_t* obj_set_v;
    b_heap_t* U;
    int rear_u, top_u, i;
    void* e;
    node_t* node_v;
    bst_node_t* bst_node_v;
    BIT_TYPE p_list;
    obj_t* obj_v;

    obj_set_v = alloc_obj_set();

    U = alloc_b_heap(INI_HEAP_SIZE);

    rear_u = 1;
    U->obj_arr[rear_u].element = (void*)IRTree_v.root;
    U->obj_arr[rear_u].e_tag = 1;
    U->obj_arr[rear_u].key = 0;

    b_h_insert(U, rear_u++);
    while (!b_h_is_empty(U)) {
        top_u = b_h_get_top(U);
        e = U->obj_arr[top_u].element;

        if (U->obj_arr[top_u].e_tag == 1) {
            //e is an node_t*.
            node_v = (node_t*)e;

            bst_node_v = bst_search(node_v->bst_v, key);
            if (bst_node_v == NULL)
                continue;

            p_list = bst_node_v->p_list;
            for (i = 0; i < node_v->num; i++) {
                //Check the c_key keyword.
                if (!get_k_bit(p_list, i))
                    continue;

                U->obj_arr[rear_u].element = node_v->child[i];
                if (node_v->level > 0) {
                    //node_v is an inner-node.
                    U->obj_arr[rear_u].e_tag = 1;
                    U->obj_arr[rear_u].key = 0;
                } else {
                    //node_v is a leaf-node.
                    U->obj_arr[rear_u].e_tag = 2;
                    U->obj_arr[rear_u].key = 0;
                }

                //Enqueue.
                b_h_insert(U, rear_u++);
            } //
        } else {
            //e is an obj_t*.
            obj_v = (obj_t*)e;
            add_obj_set_entry(obj_v, obj_set_v);
        }
    }

    release_b_heap(U);

    return obj_set_v;
}

//===========================================================================================
/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 *
 *  1. objects are sorted by distance to loc_v, closest to farthest
 *  2. obj->dist storing the distance
 */
obj_set_t* range_query_sorted_dist(disk_t* disk_v, query_t* q, loc_t* loc_v)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    range_query_sorted_dist_sub(IRTree_v.root, disk_v, obj_set_v, q, loc_v);

    return obj_set_v;
}

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_sorted_dist_sub(node_t* node_v, disk_t* disk_v, obj_set_t*& obj_set_v, query_t* q, loc_t* loc_v)
{
    int i;
    BIT_TYPE p_list;
    range* MBR;

    if (node_v->parent == NULL)
        MBR = get_MBR_node(node_v, IRTree_v.dim);
    else
        MBR = node_v->parent->MBRs[node_v->loc];

    //No overlapping.
    if (!is_overlap(MBR, disk_v))
        return;

    //Enclosed entrely.
    if (is_enclosed(MBR, disk_v)) {
        retrieve_sub_tree_sorted_dist(node_v, obj_set_v, q, loc_v);
        if (node_v->parent == NULL) {
            free(MBR);

            /*s*/
            stat_v.memory_v -= IRTree_v.dim * sizeof(range);
            /*s*/
        }

        return;
    }

    //The remaining cases.
    if (node_v->level == 0) //node_v is a leaf-node.
    {
        ///for each object inside the leaf node
        for (i = 0; i < node_v->num; i++) {
            ///if inside the disk and relevant, then add into obj_set_v
            if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v) && is_relevant_obj((obj_t*)(node_v->child[i]), q)) {
                add_obj_set_entry_sorted((obj_t*)(node_v->child[i]), obj_set_v, calc_minDist(((obj_t*)(node_v->child[i]))->MBR, loc_v));
            }
        }
    } else //node_v is an inner-node.
    {
        ///retrieve a list of relevant childean
        p_list = is_relevant_node(node_v, q);

        ///recursive call for each child in the list
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                range_query_sorted_dist_sub((node_t*)(node_v->child[i]), disk_v, obj_set_v, q, loc_v);
        }
    }

    if (node_v->parent == NULL) {
        free(MBR);
        /*s*/
        stat_v.memory_v -= IRTree_v.dim * sizeof(range);
        /*s*/
    }
}

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */
void retrieve_sub_tree_sorted_dist(node_t* node_v, obj_set_t*& obj_set_v, query_t* q, loc_t* loc_v)
{
    int i;
    BIT_TYPE p_list;

    if (node_v->level == 0) {
        //node_v is a leaf-node.
        //Retrieve all its objects.
        for (i = 0; i < node_v->num; i++) {
            if (is_relevant_obj((obj_t*)(node_v->child[i]), q))
                add_obj_set_entry_sorted((obj_t*)(node_v->child[i]), obj_set_v, calc_minDist(((obj_t*)(node_v->child[i]))->MBR, loc_v));
        }
    } else {
        //node_v is an inner-node.
        //Invoke the function recursively.
        p_list = is_relevant_node(node_v, q);
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                retrieve_sub_tree_sorted_dist((node_t*)(node_v->child[i]), obj_set_v, q, loc_v);
        }
    }
}

//===========================================================================================
/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 *
 *  1. objects are sorted by distance to loc_v, closest to farthest
 *  2. obj->dist storing the distance
 */
obj_set_t* range_query_sorted_dist(disk_t* disk_v1, disk_t* disk_v2, query_t* q, loc_t* loc_v)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    range_query_sorted_dist_sub(IRTree_v.root, disk_v1, disk_v2, obj_set_v, q, loc_v);

    return obj_set_v;
}

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_sorted_dist_sub(node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t*& obj_set_v, query_t* q, loc_t* loc_v)
{
    int i;
    BIT_TYPE p_list;
    range* MBR;

    if (node_v->parent == NULL)
        MBR = get_MBR_node(node_v, IRTree_v.dim);
    else
        MBR = node_v->parent->MBRs[node_v->loc];

    //No overlapping.
    if (!is_overlap(MBR, disk_v1) || !is_overlap(MBR, disk_v2))
        return;

    //Enclosed entrely.
    if (is_enclosed(MBR, disk_v1) && is_enclosed(MBR, disk_v2)) {
        retrieve_sub_tree_sorted_dist(node_v, obj_set_v, q, loc_v);
        if (node_v->parent == NULL) {
            free(MBR);

            /*s*/
            stat_v.memory_v -= IRTree_v.dim * sizeof(range);
            /*s*/
        }

        return;
    }

    //The remaining cases.
    if (node_v->level == 0) //node_v is a leaf-node.
    {
        ///for each object inside the leaf node
        for (i = 0; i < node_v->num; i++) {
            ///if inside the disk and relevant, then add into obj_set_v
            if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v1) && is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v2) && is_relevant_obj((obj_t*)(node_v->child[i]), q)) {
                add_obj_set_entry_sorted((obj_t*)(node_v->child[i]), obj_set_v, calc_minDist(((obj_t*)(node_v->child[i]))->MBR, loc_v));
            }
        }
    } else //node_v is an inner-node.
    {
        ///retrieve a list of relevant childean
        p_list = is_relevant_node(node_v, q);

        ///recursive call for each child in the list
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                range_query_sorted_dist_sub((node_t*)(node_v->child[i]), disk_v1, disk_v2, obj_set_v, q, loc_v);
        }
    }

    if (node_v->parent == NULL) {
        free(MBR);
        /*s*/
        stat_v.memory_v -= IRTree_v.dim * sizeof(range);
        /*s*/
    }
}

//===========================================================================================
/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query_sorted_dist(query_t* q, loc_t* loc_v)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    retrieve_sub_tree_sorted_dist(IRTree_v.root, obj_set_v, q, loc_v);

    return obj_set_v;
}

//===========================================================================================
// range query for intersection of two circle (disk) region

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_sub(node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t*& obj_set_v, query_t* q)
{
    int i;
    BIT_TYPE p_list;
    range* MBR;

    if (node_v->parent == NULL)
        MBR = get_MBR_node(node_v, IRTree_v.dim);
    else
        MBR = node_v->parent->MBRs[node_v->loc];

    //No overlapping.
    if (!is_overlap(MBR, disk_v1) || !is_overlap(MBR, disk_v2))
        return;

    //Enclosed entrely.
    if (is_enclosed(MBR, disk_v1) && is_enclosed(MBR, disk_v2)) {
        retrieve_sub_tree(node_v, obj_set_v, q);
        if (node_v->parent == NULL) {
            free(MBR);

            /*s*/
            stat_v.memory_v -= IRTree_v.dim * sizeof(range);
            /*s*/
        }

        return;
    }

    //The remaining cases.
    if (node_v->level == 0) //node_v is a leaf-node.
    {
        ///for each object inside the leaf node
        for (i = 0; i < node_v->num; i++) {
            ///if inside the disk and relevant, then add into obj_set_v
            if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v1) && is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v2) && is_relevant_obj((obj_t*)(node_v->child[i]), q))
                add_obj_set_entry((obj_t*)(node_v->child[i]), obj_set_v);
        }
    } else //node_v is an inner-node.
    {
        ///retrieve a list of relevant childean
        p_list = is_relevant_node(node_v, q);

        ///recursive call for each child in the list
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i))
                range_query_sub((node_t*)(node_v->child[i]), disk_v1, disk_v2, obj_set_v, q);
        }
    }

    if (node_v->parent == NULL) {
        free(MBR);
        /*s*/
        stat_v.memory_v -= IRTree_v.dim * sizeof(range);
        /*s*/
    }
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query(disk_t* disk_v1, disk_t* disk_v2, query_t* q)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    range_query_sub(IRTree_v.root, disk_v1, disk_v2, obj_set_v, q);

    return obj_set_v;
}

/*
 * Implementation of query generation in Cao TODS15
 *
 */
query_t** gen_query_set4(int query_n, int n, data_t* data_v)
{
    int i, j, k, rand_tmp, rand_tmp2;
    int* rand_v; ///array storing the randomed number
    loc_t* loc_v;
    psi_t* psi_v;
    query_t** q_set;
    k_node_t* k_head;

    /*t*/

    //Generate the queries.
    q_set = (query_t**)malloc(sizeof(query_t*) * query_n);
    memset(q_set, 0, sizeof(query_t*) * query_n);

    ///for each query
    for (i = 0; i < query_n; i++) {
        // printf( "------- Gen_query %d ------\n",i);

        ///allocate memory
        q_set[i] = alloc_query();
        psi_v = alloc_psi();
        rand_v = (int*)malloc(n * sizeof(int));
        memset(rand_v, 0, sizeof(int));

        /// find n objects
        for (k = 0; k < n; k++) {
            // printf( "k=%d\n",k);
            do {
                //-------- find an obj randomly ------
                rand_tmp = rand_i(0, data_v->obj_n - 1);
                //   printf("rand_tmp:%d id:%d\n", rand_tmp, data_v->obj_v[rand_tmp].id);

                while (is_old(rand_v, k, rand_tmp))
                    rand_tmp = rand_i(0, data_v->obj_n - 1);

                rand_v[k] = rand_tmp;
                k_head = data_v->obj_v[rand_tmp].k_head->next;

                //-------- find a keyword in this obj randomly ------
                int cnt = 0; ///count the number keywords in this object
                while (k_head != NULL) {
                    cnt++;
                    k_head = k_head->next;
                }
                rand_tmp2 = rand_i(1, cnt);
                //    printf("rand_tmp2:%d \t cnt:%d\n", rand_tmp2, cnt);

                k_head = data_v->obj_v[rand_tmp].k_head->next;
                int pos = 1;
                while (k_head != NULL && pos < rand_tmp2) {
                    pos++;
                    k_head = k_head->next;
                }
                //    printf("key:%f\n",k_head->key);
            } while (is_old2(psi_v, k_head->key)); //find another object if the keyword is already in psi_v

            add_psi_entry(psi_v, k_head->key);
        }

        q_set[i]->psi_v = psi_v;

        //print_k_list(psi_v->k_head, stdout);

        //---------query location------------------
        ///center point of the n objects
        loc_v = alloc_loc(data_v->dim);
        //for each dim
        for (j = 0; j < data_v->dim; j++) {
            double sum = 0.0;
            for (int x = 0; x < n; x++) {
                sum += data_v->obj_v[rand_v[x]].MBR[j].min;
            }
            loc_v->coord[j] = sum / n;
        }
        q_set[i]->loc_v = loc_v;

        // print_loc(loc_v, stdout);

        free(rand_v);
    }

    return q_set;
}

/******
 *	check whether the key @key is in the list @psi_v
 ******/
bool is_old2(psi_t* psi_v, B_KEY_TYPE key)
{
    k_node_t* k_node_v;

    k_node_v = psi_v->k_head->next;
    while (k_node_v != NULL) {
        if (key == k_node_v->key)
            return true;

        k_node_v = k_node_v->next;
    }

    return false;
}

/*
 *	Sort the objects in @obj_set_v by their distances to @loc_v.
 *
 *	Method: heap-sort.
 *	Alternative method: based on the binary search tree.
 */
b_heap_t* heap_sort_obj_set(obj_set_t* obj_set_v, loc_t* loc_v2)
{
    int cur;
    B_KEY_TYPE dist;
    b_heap_t* b_h;
    obj_node_t* obj_node_v;
    loc_t* loc_v;

    b_h = alloc_b_heap(obj_set_v->obj_n + 1);

    cur = 1;
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL) {
        loc_v = get_obj_loc(obj_node_v->obj_v);
        dist = calc_dist_loc(loc_v, loc_v2);
        release_loc(loc_v);

        b_h->obj_arr[cur].key = dist;
        b_h->obj_arr[cur].obj_v = obj_node_v->obj_v;

        b_h_insert(b_h, cur);

        cur++;

        obj_node_v = obj_node_v->next;
    }

    return b_h;
}

/*
 *	Sort the objects in @obj_set_v by their distances to @loc_v.
 * *** object with dist > dist_max will NOT be added into the heap
 *
 *	Method: heap-sort.
 *	Alternative method: based on the binary search tree.
 */
b_heap_t* heap_sort_obj_set2(obj_set_t* obj_set_v, loc_t* loc_v2, B_KEY_TYPE dist_max)
{
    int cur;
    B_KEY_TYPE dist;
    b_heap_t* b_h;
    obj_node_t* obj_node_v;
    loc_t* loc_v;

    b_h = alloc_b_heap(obj_set_v->obj_n + 1);

    cur = 1;
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL) {
        loc_v = get_obj_loc(obj_node_v->obj_v);
        dist = calc_dist_loc(loc_v, loc_v2);
        release_loc(loc_v);

        if (dist > dist_max) {
            obj_node_v = obj_node_v->next;
            continue;
        }

        b_h->obj_arr[cur].key = dist;
        b_h->obj_arr[cur].obj_v = obj_node_v->obj_v;

        b_h_insert(b_h, cur);

        cur++;

        obj_node_v = obj_node_v->next;
    }

    return b_h;
}

/*
 *	Sort the objects in @obj_set_v by their distances to @q.
 *
 *	Method: heap-sort.
 *	Alternative method: based on the binary search tree.
 *
 *  maxDist: max d(o,q) for the objects in the set
 *
 */
b_heap_t* heap_sort_obj_set3(obj_set_t* obj_set_v, query_t* q, B_KEY_TYPE& maxDist)
{
    int cur;
    B_KEY_TYPE dist;
    b_heap_t* b_h;
    obj_node_t* obj_node_v;
    loc_t* loc_v;

    b_h = alloc_b_heap(obj_set_v->obj_n + 1);

    cur = 1;
    obj_node_v = obj_set_v->head->next;
    while (obj_node_v != NULL) {
        loc_v = get_obj_loc(obj_node_v->obj_v);
        dist = calc_dist_loc(loc_v, q->loc_v);
        release_loc(loc_v);

        //--
        if (dist > maxDist)
            maxDist = dist;
        //--
        b_h->obj_arr[cur].key = dist;
        b_h->obj_arr[cur].obj_v = obj_node_v->obj_v;

        b_h_insert(b_h, cur);

        cur++;

        obj_node_v = obj_node_v->next;
    }

    return b_h;
}

/*
 *	Constrained NN search.
 *	Find the NN that contains the keyword @key and
 *	is located in the dist @disk_v.
 *
 *	disk_v == NULL indicates that the disk is the whole space,
 *	thus, in this case, this constraint is invalid (interface consideration).
 *
 *	Method: best-first search.
 *
 * Modified from const_NN_key
 * Finding NN that is located in intersection of 2 disks
 */
obj_t* const_NN_key2(loc_t* loc_v, KEY_TYPE key, disk_t* disk_v1, disk_t* disk_v2)
{
    int size, top, i, rear;
    BIT_TYPE p_list;

    B_KEY_TYPE min_dist, c_dist;
    range* MBR;
    obj_t* obj_v;
    node_t* node_v;
    b_heap_t* b_h;

    //Checking 1.
    if (IRTree_v.root->num == 0)
        return NULL;

    //Checking 2: involved in the search process implicitly.
    //if( !( bst_node_v = bst_search( IRTree.root->bst_v, key)))
    //return NULL;

    //Checking 3:
    MBR = get_MBR_node(IRTree_v.root, IRTree_v.dim);
    if (disk_v1 != NULL && !is_overlap(MBR, disk_v1) && disk_v2 != NULL && !is_overlap(MBR, disk_v2))
        return NULL;

    size = IRTree_v.obj_n / M + 1;
    b_h = alloc_b_heap(size);

    rear = 1;
    b_h->obj_arr[rear].node_v = IRTree_v.root;
    b_h->obj_arr[rear].key = calc_minDist(MBR, loc_v);
    free(MBR);

    /*s*/
    stat_v.memory_v -= sizeof(IRTree_v.dim * sizeof(range));
    /*s*/

    b_h_insert(b_h, rear++);

    min_dist = DBL_MAX;
    obj_v = NULL;

    while (!b_h_is_empty(b_h)) {
        top = b_h_get_top(b_h);
        if (b_h->obj_arr[top].key >= min_dist)
            break;

        node_v = b_h->obj_arr[top].node_v;

        //Keyword constraint (for all its entries).
        //if( !( bst_node_v = bst_search( node_v->bst_v, key)))
        //continue;
        if ((p_list = has_key_node(node_v, key)) == 0)
            continue;

        //bst_node_v != NULL.
        for (i = 0; i < node_v->num; i++) {
            //Keyword constraint (for a specific entry).
            if (!get_k_bit(p_list, i))
                continue;

            //the i'th entry contains the keyword.
            if (node_v->level == 0) {
                //node_v is a leaf-node.
                MBR = ((obj_t*)(node_v->child[i]))->MBR;

                //Disk constraint (if any, i.e., disk_v != NULL).
                if (disk_v1 != NULL && !is_overlap(MBR, disk_v1) && disk_v2 != NULL && !is_overlap(MBR, disk_v2))
                    continue;

                c_dist = calc_minDist(MBR, loc_v);
                if (c_dist < min_dist) {
                    //Update the current best solution.
                    min_dist = c_dist;
                    obj_v = (obj_t*)(node_v->child[i]);
                }
            } else {
                //node_v is an inner-node.
                MBR = node_v->MBRs[((node_t*)(node_v->child[i]))->loc];

                //Disk constraint (if any).
                if (disk_v1 != NULL && !is_overlap(MBR, disk_v1) && disk_v2 != NULL && !is_overlap(MBR, disk_v2))
                    continue;

                c_dist = calc_minDist(MBR, loc_v);
                if (c_dist < min_dist) {
                    //Enqueue the child.
                    b_h->obj_arr[rear].node_v = (node_t*)(node_v->child[i]);
                    b_h->obj_arr[rear].key = c_dist;

                    b_h_insert(b_h, rear++);
                }
            }
        } //for
    } //while

    release_b_heap(b_h);

    return obj_v;
}

/*
 *	Check the distance constraint.
 */
bool check_dist_constraint(obj_set_t* obj_set_v, obj_t* obj_v, B_KEY_TYPE d)
{
    obj_node_t* obj_node_iter;

    obj_node_iter = obj_set_v->head->next;
    while (obj_node_iter != NULL) {
        if (calc_dist_obj(obj_node_iter->obj_v, obj_v) > d)
            return false;

        obj_node_iter = obj_node_iter->next;
    }

    return true;
}
/*
 *	Update the heap structure @heap by removing the
 *  objects having number of intersection = 0
 */
void update_heap_obj(bst_t* heap, bst_node_t* bst_node_v, psi_t* psi_v, query_t* q, bst_node_list_t* bst_node_list_v, psi_t* psi_excluded)
{
    bst_node_list_t* tmp;

    if (bst_node_v == NULL)
        return;

    if (bst_node_v->left != NULL)
        update_heap_obj(heap, bst_node_v->left, psi_v, q, bst_node_list_v, psi_excluded);
    if (bst_node_v->right != NULL)
        update_heap_obj(heap, bst_node_v->right, psi_v, q, bst_node_list_v, psi_excluded);

    //--pre checking: updated is not needed for objects do not contain excluded keyword in the current iteration
    bool flag = false;
    k_node_t* k_head = psi_excluded->k_head->next;
    while (k_head != NULL) {
        if (has_key_obj(bst_node_v->obj_v1, k_head->key)) {
            flag = true;
            break;
        }
        k_head = k_head->next;
    }
    if (!flag)
        return;

    //--

    int numOfInt = is_relevant_obj(bst_node_v->obj_v1, psi_v);
    if (numOfInt == 0) {
        bst_delete(heap, bst_node_v);
    } else {
        B_KEY_TYPE key_new = calc_minDist(bst_node_v->obj_v1->MBR, q->loc_v) / ((double)numOfInt);
        //   printf("old:%0.4lf\t new:%0.4lf \n",bst_node_v->key,key_new );
        //if key is updated, del then re-insert
        if (bst_node_v->key != key_new) {
            bst_node_v->key = key_new;
            bst_delete(heap, bst_node_v);
            //bug.
            bst_node_v->p = NULL;
            bst_node_v->left = NULL;
            bst_node_v->right = NULL;

            tmp = (bst_node_list_t*)malloc(sizeof(bst_node_list_t));
            memset(tmp, 0, sizeof(bst_node_list_t));

            tmp->bst_node_v = bst_node_v;
            tmp->next = bst_node_list_v->next;
            bst_node_list_v->next = tmp;
        }
    }

    return;
}

//========================================================================
//range query + remove object dominated

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */
void retrieve_sub_tree_dominant(node_t* node_v, obj_set_t*& obj_set_v, query_t* q)
{
    int i;
    BIT_TYPE p_list;

    if (node_v->level == 0) {
        //node_v is a leaf-node.
        //Retrieve all its objects.
        for (i = 0; i < node_v->num; i++) {
            if (is_relevant_obj((obj_t*)(node_v->child[i]), q))
                add_obj_set_entry_domanint((obj_t*)(node_v->child[i]), obj_set_v, q);
        }
    } else {
        //node_v is an inner-node.
        //Invoke the function recursively.
        p_list = is_relevant_node(node_v, q);
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i) && !check_dominant((node_t*)(node_v->child[i]), obj_set_v, q))
                retrieve_sub_tree_dominant((node_t*)(node_v->child[i]), obj_set_v, q);
        }
    }
}

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_dominant_sub(node_t* node_v, disk_t* disk_v, obj_set_t*& obj_set_v, query_t* q)
{
    int i;
    BIT_TYPE p_list;
    range* MBR;

    if (node_v->parent == NULL)
        MBR = get_MBR_node(node_v, IRTree_v.dim);
    else
        MBR = node_v->parent->MBRs[node_v->loc];

    //No overlapping.
    if (!is_overlap(MBR, disk_v))
        return;

    //Enclosed entrely.
    if (is_enclosed(MBR, disk_v)) {
        retrieve_sub_tree_dominant(node_v, obj_set_v, q);
        if (node_v->parent == NULL) {
            free(MBR);

            /*s*/
            stat_v.memory_v -= IRTree_v.dim * sizeof(range);
            /*s*/
        }

        return;
    }

    //The remaining cases.
    if (node_v->level == 0) //node_v is a leaf-node.
    {
        ///for each object inside the leaf node
        for (i = 0; i < node_v->num; i++) {
            ///if inside the disk and relevant, then add into obj_set_v
            if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v) && is_relevant_obj((obj_t*)(node_v->child[i]), q))
                add_obj_set_entry_domanint((obj_t*)(node_v->child[i]), obj_set_v, q);
        }
    } else //node_v is an inner-node.
    {
        ///retrieve a list of relevant childean
        p_list = is_relevant_node(node_v, q);

        ///recursive call for each child in the list
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i) && !check_dominant((node_t*)(node_v->child[i]), obj_set_v, q))
                range_query_dominant_sub((node_t*)(node_v->child[i]), disk_v, obj_set_v, q);
        }
    }

    if (node_v->parent == NULL) {
        free(MBR);
        /*s*/
        stat_v.memory_v -= IRTree_v.dim * sizeof(range);
        /*s*/
    }
}

/*
 * *** only dominant objects in the obj_set_v
 * *** the object set is sorted by d(o,q) (from closest to farthest)
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query_dominant(disk_t* disk_v, query_t* q)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    range_query_dominant_sub(IRTree_v.root, disk_v, obj_set_v, q);

    return obj_set_v;
}

/*
 *	Add an object entry @obj_v to @obj_set_v.
 */
void add_obj_set_entry_domanint(obj_t* obj_v, obj_set_t* obj_set_v, query_t* q)
{
    obj_node_t *obj_node_v, *obj_node_temp, *obj_node_new;

    //we are actually checking obj_node_v->next in each iteration
    obj_node_v = obj_set_v->head;
    while (obj_node_v->next != NULL) {
        bool remove = false; //incdicate whether obj_node_v->next is removed or not

        obj_t* obj_v2 = obj_node_v->next->obj_v;

        loc_t* loc_v1 = get_obj_loc(obj_v);
        loc_t* loc_v2 = get_obj_loc(obj_v2);
        B_KEY_TYPE d_o1_q = calc_dist_loc(loc_v1, q->loc_v);
        B_KEY_TYPE d_o2_q = calc_dist_loc(loc_v2, q->loc_v);
        release_loc(loc_v1);
        release_loc(loc_v2);

        if (d_o1_q > d_o2_q) {
            if (isDominated(obj_v, obj_v2, q->psi_v)) {
                //obj_v is dominated by an object in obj_set_v
                //no need to add this obj_v into obj_set_v
                return;
            }
        } else {
            if (isDominated(obj_v2, obj_v, q->psi_v)) {
                //the object in obj_set_v is dominated by the new obj obj_v
                //remove the object in obj_set_v
                obj_node_temp = obj_node_v->next;
                obj_node_v->next = obj_node_v->next->next;
                free(obj_node_temp);
                obj_set_v->obj_n--;
                remove = true;
            }
        }
        if (!remove)
            obj_node_v = obj_node_v->next;
    }
    //--------------

    //obj_v is not dominated by any object in obj_set_v
    //can insert

    //-----insert at head (no ordering)----
    obj_node_v = (obj_node_t*)malloc(sizeof(obj_node_t));
    memset(obj_node_v, 0, sizeof(obj_node_t));

    obj_node_v->obj_v = obj_v;
    obj_node_v->next = obj_set_v->head->next;
    obj_set_v->head->next = obj_node_v;
    obj_set_v->obj_n++;

    /*s*/
    stat_v.memory_v += sizeof(obj_node_t);
    if (stat_v.memory_v > stat_v.memory_max)
        stat_v.memory_max = stat_v.memory_v;
    /*s*/

    ///-------insert to sorted order-----
    //
    //    B_KEY_TYPE dist = calc_minDist(obj_v->MBR, q->loc_v);
    //
    //    obj_node_v = obj_set_v->head;
    //    while (obj_node_v->next != NULL){
    //        if ( dist < obj_node_v->next->dist )
    //            break;
    //        obj_node_v = obj_node_v->next;
    //    }
    //    //obj_node_v pointing to the node before insert position
    //    //-----
    //    obj_node_new = ( obj_node_t*)malloc( sizeof( obj_node_t));
    //    memset( obj_node_new, 0, sizeof( obj_node_t));
    //
    //    /*s*/
    //    stat_v.memory_v += sizeof( obj_node_t);
    //    if( stat_v.memory_v > stat_v.memory_max)
    //        stat_v.memory_max = stat_v.memory_v;
    //    /*s*/
    //
    //    obj_node_new->obj_v = obj_v;
    //    obj_node_new->dist = dist;
    //
    //    obj_node_new->next = obj_node_v->next;
    //    obj_node_v->next = obj_node_new;
    //
    //    obj_set_v->obj_n++;
    //    //-----
}

/*
 * check whether node_v is dominated by any single object in obj_set_v
 *
 */
bool check_dominant(node_t* node_v, obj_set_t* obj_set_v, query_t* q)
{
    obj_node_t* obj_node_v;

    //we are actually checking obj_node_v->next in each iteration
    obj_node_v = obj_set_v->head;
    while (obj_node_v->next != NULL) {

        obj_t* obj_v2 = obj_node_v->next->obj_v;

        loc_t* loc_v2 = get_obj_loc(obj_v2);
        B_KEY_TYPE d_o1_q = calc_minDist_node(node_v, q->loc_v);
        B_KEY_TYPE d_o2_q = calc_dist_loc(loc_v2, q->loc_v);
        release_loc(loc_v2);

        if (d_o1_q > d_o2_q && isDominated(node_v, obj_v2, q->psi_v)) {
            //node_v is dominated by an object in obj_set_v
            return true;
        }

        obj_node_v = obj_node_v->next;
    }

    return false;
}

/* check whether obj_v is dominated by obj_v2 */
bool isDominated(node_t* node_v, obj_t* obj_v2, psi_t* psi_v)
{

    k_node_t* k_node_v;
    B_KEY_TYPE key;

    k_node_v = psi_v->k_head->next;

    //for each keyword in  psi_v
    while (k_node_v != NULL) {
        key = k_node_v->key;
        //if there exist a keyword that is covered by obj_v but not obj_v2
        //obj_v is not dominated by obj_v2
        if (has_key_node(node_v, key) && !has_key_obj(obj_v2, key))
            return false;
        k_node_v = k_node_v->next;
    }
    return true;
}

//========================================================================
//range query + remove object dominated + 2 disks intersection

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_dominant_sub(node_t* node_v, disk_t* disk_v1, disk_t* disk_v2, obj_set_t*& obj_set_v, query_t* q)
{
    int i;
    BIT_TYPE p_list;
    range* MBR;

    if (node_v->parent == NULL)
        MBR = get_MBR_node(node_v, IRTree_v.dim);
    else
        MBR = node_v->parent->MBRs[node_v->loc];

    //No overlapping.
    if (!is_overlap(MBR, disk_v1) || !is_overlap(MBR, disk_v2))
        return;

    //Enclosed entrely.
    if (is_enclosed(MBR, disk_v1) && is_enclosed(MBR, disk_v1)) {
        retrieve_sub_tree_dominant(node_v, obj_set_v, q);
        if (node_v->parent == NULL) {
            free(MBR);

            /*s*/
            stat_v.memory_v -= IRTree_v.dim * sizeof(range);
            /*s*/
        }

        return;
    }

    //The remaining cases.
    if (node_v->level == 0) //node_v is a leaf-node.
    {
        ///for each object inside the leaf node
        for (i = 0; i < node_v->num; i++) {
            ///if inside the disk and relevant, then add into obj_set_v
            if (is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v1) && is_enclosed(((obj_t*)(node_v->child[i]))->MBR, disk_v2) && is_relevant_obj((obj_t*)(node_v->child[i]), q))
                add_obj_set_entry_domanint((obj_t*)(node_v->child[i]), obj_set_v, q);
        }
    } else //node_v is an inner-node.
    {
        ///retrieve a list of relevant childean
        p_list = is_relevant_node(node_v, q);

        ///recursive call for each child in the list
        for (i = 0; i < node_v->num; i++) {
            if (get_k_bit(p_list, i) && !check_dominant((node_t*)(node_v->child[i]), obj_set_v, q))
                range_query_dominant_sub((node_t*)(node_v->child[i]), disk_v1, disk_v2, obj_set_v, q);
        }
    }

    if (node_v->parent == NULL) {
        free(MBR);
        /*s*/
        stat_v.memory_v -= IRTree_v.dim * sizeof(range);
        /*s*/
    }
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query_dominant(disk_t* disk_v1, disk_t* disk_v2, query_t* q)
{
    obj_set_t* obj_set_v;

    obj_set_v = alloc_obj_set();
    range_query_dominant_sub(IRTree_v.root, disk_v1, disk_v2, obj_set_v, q);

    return obj_set_v;
}

//=======================
/* check_cost
 * check whether adding the new object obj_v to S will satisfy the pdist_max constraint
 */

bool check_cost(obj_t* obj_v, obj_set_t* S, B_KEY_TYPE pdist_max)
{
    obj_node_t* obj_node_v = S->head->next;

    while (obj_node_v != NULL) {
        if (calc_dist_obj(obj_v, obj_node_v->obj_v) > pdist_max)
            return false;
        obj_node_v = obj_node_v->next;
    }
    return true;
}

//----

/* check whether obj_v is dominated by obj_v2 */
bool isDominated(obj_t* obj_v, obj_t* obj_v2, psi_t* psi_v)
{

    k_node_t* k_node_v;
    B_KEY_TYPE key;

    k_node_v = psi_v->k_head->next;

    //for each keyword in  psi_v
    while (k_node_v != NULL) {
        key = k_node_v->key;
        //if there exist a keyword that is covered by obj_v but not obj_v2
        //obj_v is not dominated by obj_v2
        if (has_key_obj(obj_v, key) && !has_key_obj(obj_v2, key))
            return false;
        k_node_v = k_node_v->next;
    }
    return true;
}

/*
 *	Remove the specified entry from the list of objects @obj_set_v.
 *  Modified from remove_obj_set_entry
 */
void remove_obj_set_entry(obj_set_t* obj_set_v, obj_t* obj_v)
{
    obj_node_t *obj_node_prev, *obj_node_temp;

    obj_node_prev = obj_set_v->head;
    // printf("finding object to remove\n");

    while (obj_node_prev->next != NULL) {
        if (obj_node_prev->next->obj_v == obj_v) {
            obj_node_temp = obj_node_prev->next;
            obj_node_prev->next = obj_node_prev->next->next;
            obj_set_v->obj_n--;
            free(obj_node_temp);
            //  printf("an object removed\n");
            /*s*/
            stat_v.memory_v -= sizeof(obj_node_t);
            /*s*/

            return;
        }

        obj_node_prev = obj_node_prev->next;
    }
}

/*
 *	Release the binary search tree T.
 * modified from bst_release
 */
void release_IF(bst_t* T)
{
    if (T != NULL) {
        if (T->root != NULL)
            release_IF_sub(T->root);

        free(T);

        /*s*/
        stat_v.memory_v -= sizeof(bst_t);
        /*s*/
    }
}

void release_IF_sub(bst_node_t* x)
{
    if (x->left != NULL)
        bst_release_sub(x->left);
    if (x->right != NULL)
        bst_release_sub(x->right);

    release_obj_set(x->p_list_obj);
    free(x);

    /*s*/
    stat_v.memory_v -= sizeof(bst_node_t);
    /*s*/
}
