//
//  sum_alg.cpp
//

#include "sum_alg.h"

/*
 * Sum_Exact_new
 * CCG+15
 */
obj_set_t* Sum_Exact_new( query_t* q)
{
    
	int  rear, top;
	b_heap_t* U;
	void* e;
	node_t* node_v;
	obj_t* obj_v;
    psi_t* psi_v;
    bit_node_t* ks,* bit_node_s,* bit_node_prev, * bit_node_ss;
    bit_set_t* markedSet,* valuedSet,* ssSet ;


    int qSize = q->psi_v->key_n;
    int powerSetSize = pow(2, qSize)-1;
    double cost[powerSetSize+1];
    obj_set_t* group[powerSetSize+1];
    
    for(int i=0;i<=powerSetSize;i++) {
        cost[i] =INFINITY;
        group[i] = NULL;
    }
    
    markedSet=alloc_bit_set();
    valuedSet=alloc_bit_set();
    //----creating an array storing the mapping between query keyword and index of bit
    
    k_node_t* k_node_v;
    B_KEY_TYPE A[qSize];
    k_node_v = q->psi_v->k_head->next;
    
    for( int i=0;i<qSize;i++){
        A[i] = k_node_v->key;
        k_node_v=k_node_v->next;
    }
    //----------------------------------
    
    
	U = alloc_b_heap( INI_HEAP_SIZE);
    
	rear = 1;
	U->obj_arr[ rear].element = ( void*)IRTree_v.root;
	U->obj_arr[ rear].e_tag = 1;
	U->obj_arr[ rear].key = calc_minDist_node( IRTree_v.root, q->loc_v);
    
	b_h_insert( U, rear++);
    
    while(!b_h_is_empty( U))
    {
        /// printf("enter while\n");
        top = b_h_get_top( U);
        e = U->obj_arr[ top].element;
        
        ///-------ks = q.\psi \cap e.\psi----------------------
       if( U->obj_arr[ top].e_tag == 2){
            //e is an obj_t*.
            obj_v = ( obj_t*)e;
            psi_v = key_intersection( q->psi_v->k_head, obj_v->k_head);
            
        }else{
            //e is an node_t*.
            node_v = (node_t*)e;
            psi_v = node_intersection(node_v, q->psi_v);
        }
        
        ks = psi_to_bit_node(A, qSize, psi_v);
        release_psi(psi_v);

        ///---------------------------------------------------------
        ///if(ks in markedSet)
        if(bit_set_find(markedSet, ks)){
            free(ks);
            continue;
        }
        ///---------------------------------------------------------
        
        if( U->obj_arr[ top].e_tag == 2){
            //e is an obj_t*.
            obj_v = ( obj_t*)e;
            
            ///---------------------------------------------------------
            ///for each s in valued set
            //   printf("Step 1\n");
            
            bit_node_s=valuedSet->p_head->next;
            bit_node_prev=valuedSet->p_head;
            
            while( bit_node_s!=NULL){
                unsigned int i=map_to_integer(bit_node_s);
                
                loc_t* loc_v = get_obj_loc(obj_v);
                
                if (cost[i]< calc_dist_loc(loc_v, q->loc_v)){
                    if(i == powerSetSize){
                        release_loc(loc_v);
                        goto E;//return group[i]
                    }
                    //remove from valuedSet
                    bit_node_prev->next = bit_node_s->next;
                    bit_node_s->next=NULL;
                    valuedSet->bit_n--;
                    
                    //add to markedSet
                    bit_set_insert(markedSet, bit_node_s);
                    bit_node_s=bit_node_prev->next;
                    
                }else{
                    bit_node_prev=bit_node_prev->next;
                    bit_node_s=bit_node_s->next;
                }
                release_loc(loc_v);
            }

            ///---------------------------------------------------------
            //find all subset in ks
            //    printf("Step 2\n");
            
            ssSet = find_all_bit_set(ks);
            
            bit_node_ss=ssSet->p_head->next;
            bit_node_prev=ssSet->p_head;

            ///for each ss in ks
            while(bit_node_ss != NULL){
                if(!bit_set_find(markedSet, bit_node_ss)){
                    ///insert into markedSet
                    unsigned int i = map_to_integer(bit_node_ss);
                    
                    bit_node_t* bit_node_temp2 = copy_bit_node(bit_node_ss);
                    bit_set_insert(markedSet, bit_node_temp2);
                    if(bit_set_find(valuedSet, bit_node_ss))
                        bit_set_remove(valuedSet,bit_node_ss);
                    
                    //cost[i] <- Dist(e,q)
                    cost[i] = calc_minDist(obj_v->MBR, q->loc_v);
                    //   printf("i:%d\tcost[i]%0.3lf\n",i,cost[i]);
                    
                    //group[i] <- {e}
                    obj_set_t* obj_set_v = alloc_obj_set();
                    add_obj_set_entry(obj_v, obj_set_v);
                    release_obj_set(group[i]);
                    group[i] = obj_set_v;
                    
                }
                
                bit_node_ss=bit_node_ss->next;
                bit_node_prev=bit_node_prev->next;
            }
            
            release_bit_set(ssSet);
            
            ///---------------------------------------------------------
            //  printf("Step 3\n");
            
            unsigned int j = map_to_integer(ks);
            if (j == powerSetSize)
            {
                goto E;
            }
            
            for(unsigned int i=1;i<=powerSetSize;i++){
                if( cost[i]==INFINITY) continue;
                unsigned int unionKey = i;
                union_bit(unionKey, j);
                //    printf("i:%u\tj:%u\tunionKey:%u\n",i,j,unionKey);
                if(unionKey == i || unionKey==j){
                    continue;
                }
                    
                //--
                if (cost[unionKey] == INFINITY){
                    //valuedSet += unionKey
                    bit_node_t* bit_node_v = alloc_bit_node();
                    bit_node_v->bit_string = unionKey;
                    bit_node_v->next = NULL;
                    bit_set_insert(valuedSet, bit_node_v);
                }
                //--
                    
                double D = cost[i] + calc_minDist(obj_v->MBR, q->loc_v);
                //     printf("cost[unionKey]:%0.3lf\tD:%0.3lf\n",cost[unionKey],D);
                    
                if(cost[unionKey]>D){
                    //      printf("cost[unionKey]>D\n");
                    cost[unionKey] = D;
                    //group[unionkey] = group[i] \cap {e}
                    release_obj_set(group[unionKey]);
                    group[unionKey] = copy_obj_set(group[i]);
                    add_obj_set_entry(obj_v, group[unionKey]);
                }
            }

            //------------------------------------
        }
        else
        {
            //e is an node_t*.
            node_v = ( node_t*)e;
            
            //for each child node of e
            if( node_v ->level > 0)
            {
                //   printf("inserting inner node child (total number:%d)\n",node_v->num);
                
                for(int i=0; i<node_v->num; i++)
                {
                    //printf("inserting inner node (i:%d)\n",i);
                    
                    //node_v is an inner-node.
                    node_t* child = (node_t*)node_v->child[i];
                    psi_t* psi_v_temp = node_intersection(child, q->psi_v);
                    
                    bit_node_t* bit_node_v_temp = psi_to_bit_node(A, q->psi_v->key_n, psi_v_temp);
                    
                    if( psi_v_temp->key_n > 0 && !bit_set_find(markedSet, bit_node_v_temp)){
                    
                        U->obj_arr[ rear].element = node_v->child[ i];
                        U->obj_arr[ rear].e_tag = 1;
                        U->obj_arr[ rear].key = calc_minDist_node(child, q->loc_v);
                        //Enqueue.
                        b_h_insert( U, rear++);
                    }
                    release_psi(psi_v_temp);
                    free(bit_node_v_temp);
                    
                }//end for
             }
            else
            {
                for(int i=0; i<node_v->num; i++)
                {
                    //node_v is a leaf-node.
                    obj_t* child = (obj_t*)node_v->child[i];
                    psi_t* psi_v_temp = key_intersection(q->psi_v->k_head, child->k_head);
             
                    bit_node_t* bit_node_v_temp = psi_to_bit_node(A, q->psi_v->key_n, psi_v_temp);
                    
                    if( psi_v_temp->key_n >=0 && !bit_set_find(markedSet, bit_node_v_temp)){
                    
                        loc_t* loc_v = get_obj_loc(child);
                    
                        U->obj_arr[ rear].element = node_v->child[ i];
                        U->obj_arr[ rear].e_tag = 2;
                        U->obj_arr[ rear].key =calc_dist_loc(loc_v, q->loc_v);
                        //Enqueue.
                        b_h_insert( U, rear++);
                        
                        release_loc(loc_v);
                    }
                    release_psi(psi_v_temp);
                    free(bit_node_v_temp);
                
                }//end for
            }
        }//else
        free(ks);
    }//while
E:
    //release resource.
    release_bit_set(markedSet);
    release_bit_set(valuedSet);
    release_b_heap(U);

    return group[powerSetSize];
}


//===========================================================================
//===========================================================================
//===========================================================================

/*
 * Sum_Exact
 */
obj_set_t* Sum_Exact( query_t* q)
{
    
	int  rear, top;
	b_heap_t* U;
	void* e;
	node_t* node_v;
	obj_t* obj_v;
    psi_t* psi_v;
    bit_node_t* ks,* bit_node_s,* bit_node_prev, * bit_node_ss,* bit_node_ts;
    bit_set_t* markedSet,* valuedSet,* tempSet,* ssSet ;
    
    int qSize = q->psi_v->key_n;
    int powerSetSize = pow(2, qSize)-1;
    double cost[powerSetSize+1];
    obj_set_t* group[powerSetSize+1];
    
    for(int i=0;i<=powerSetSize;i++) {
        cost[i] =INFINITY;
        group[i] = NULL;
    }
    
    
    /*t/
     printf("powerSetSize:%d\n",powerSetSize);
     /*t*/
    //---
    
    markedSet=alloc_bit_set();
    valuedSet=alloc_bit_set();
    //----creating an array storing the mapping between query keyword and index of bit
    
    k_node_t* k_node_v;
    B_KEY_TYPE A[qSize];
    k_node_v = q->psi_v->k_head->next;
    
    for( int i=0;i<qSize;i++){
        A[i] = k_node_v->key;
        k_node_v=k_node_v->next;
    }
    //----------------------------------
    
    
	U = alloc_b_heap( INI_HEAP_SIZE);
    
	rear = 1;
	U->obj_arr[ rear].element = ( void*)IRTree_v.root;
	U->obj_arr[ rear].e_tag = 1;
	U->obj_arr[ rear].key = calc_minDist_node( IRTree_v.root, q->loc_v);
    
	b_h_insert( U, rear++);
    
    
    
    while(!b_h_is_empty( U))
    {
        /// printf("enter while\n");
        
        top = b_h_get_top( U);
        e = U->obj_arr[ top].element;
        /*
         printf("===============================\n");
         printf("top picked:%d\trear:%d\n",top,U->rear);
         
         printf("markedSet:%d\n",markedSet->bit_n);
         print_bit_set(markedSet);
         printf("valuedSet:%d\n",valuedSet->bit_n);
         print_bit_set(valuedSet);
         */
        /*t/
         if( e == NULL)
         printf( "");
         /*t*/
        
        ///-------ks = q.\psi \cap e.\psi----------------------
        //  printf("query:\n");
        //  print_k_list(q->psi_v->k_head, stdout);
        if( U->obj_arr[ top].e_tag == 2){
            //e is an obj_t*.
            obj_v = ( obj_t*)e;
            psi_v=key_intersection(obj_v->k_head,q->psi_v->k_head);
            
        }else{
            //e is an node_t*.
            node_v = (node_t*)e;
            k_node_t* k_node_temp = collect_keywords_bst(node_v->bst_v);
            // printf("k_node_v:\n");
            // print_k_list(k_node_v, stdout);
            psi_v=key_intersection(k_node_temp, q->psi_v->k_head);
            release_k_list(k_node_temp);
        }
        
        ks =psi_to_bit_node(A,qSize,psi_v);
        release_psi(psi_v);
        /*t/
         printf("A:\n");
         for( int a=0;a<q->psi_v->key_n;a++)
         printf("%d: %.0lf\n",a, A[a]);
         printf("psi_v:\n");
         print_k_list(psi_v->k_head, stdout);
         printf("ks:%u\n",ks->bit_string);
         /*t*/
        
        ///---------------------------------------------------------
        ///if(ks in markedSet)
        if(bit_set_find(markedSet, ks)){
            continue;
        }
        ///---------------------------------------------------------
        
        if( U->obj_arr[ top].e_tag == 2){
            // printf("if....(obj)\n");
            
            //e is an obj_t*.
            obj_v = ( obj_t*)e;
            
            ///---------------------------------------------------------
            ///for each s in valued set
            //   printf("Step 1\n");
            
            bit_node_s=valuedSet->p_head->next;
            bit_node_prev=valuedSet->p_head;
            
            while( bit_node_s!=NULL){
                unsigned int i=map_to_integer(bit_node_s);
                
                loc_t* loc_v = get_obj_loc(obj_v);
                
                if (cost[i]< calc_dist_loc(loc_v, q->loc_v)){
                    if(i==powerSetSize){
                        release_loc(loc_v);
                        goto E;//return group[i]
                    }
                    //remove from valuedSet
                    bit_node_prev->next = bit_node_s->next;
                    bit_node_s->next=NULL;
                    valuedSet->bit_n--;
                    
                    //add to markedSet
                    bit_set_insert(markedSet, bit_node_s);
                    bit_node_s=bit_node_prev->next;
                    
                }else{
                    bit_node_prev=bit_node_prev->next;
                    bit_node_s=bit_node_s->next;
                }
                release_loc(loc_v);
            }
            /*t/
             printf("markedSet:%d\n",markedSet->bit_n);
             print_bit_set(markedSet);
             printf("valuedSet:%d\n",valuedSet->bit_n);
             print_bit_set(valuedSet);
             /*t*/
            
            ///---------------------------------------------------------
            //find all subset in ks
            //    printf("Step 2\n");
            
            tempSet=alloc_bit_set();
            ssSet = find_all_bit_set(ks);
            
            bit_node_ss=ssSet->p_head->next;
            bit_node_prev=ssSet->p_head;
            
            /*t/
             printf("ks:\n");
             print_bit_node(ks);
             printf("ssSet:%d\n",ssSet->bit_n);
             print_bit_set(ssSet);
             /*t*/
            
            ///for each ss in ks
            while(bit_node_ss!=NULL){
                if(!bit_set_find(markedSet, bit_node_ss)){
                    ///insert into markedSet
                    unsigned int i=map_to_integer(bit_node_ss);
                    
                    bit_node_t* bit_node_temp2 = copy_bit_node(bit_node_ss);
                    bit_set_insert(markedSet, bit_node_temp2);
                    bit_node_t* bit_node_temp3 = copy_bit_node(bit_node_ss);
                    bit_set_insert(tempSet, bit_node_temp3);
                    if(bit_set_find(valuedSet, bit_node_ss))
                        bit_set_remove(valuedSet,bit_node_ss);
                    
                    //cost[i] <- Dist(e,q)
                    cost[i] = calc_minDist(obj_v->MBR, q->loc_v);
                    //   printf("i:%d\tcost[i]%0.3lf\n",i,cost[i]);
                    
                    //group[i] <- {e}
                    obj_set_t* obj_set_v=alloc_obj_set();
                    add_obj_set_entry(obj_v, obj_set_v);
                    release_obj_set(group[i]);
                    group[i] = obj_set_v;
                    
                    if(i==powerSetSize)
                        goto E;
                }
                
                bit_node_ss=bit_node_ss->next;
                bit_node_prev=bit_node_prev->next;
            }
            release_bit_set(ssSet);
            
            /*t/
             printf("tempSet:%d\n",tempSet->bit_n);
             print_bit_set(tempSet);
             printf("markedSet:%d\n",markedSet->bit_n);
             print_bit_set(markedSet);
             printf("valuedSet:%d\n",valuedSet->bit_n);
             print_bit_set(valuedSet);
             /*t*/
            
            ///---------------------------------------------------------
            
            ///for each ts in tempset
            //  printf("Step 3\n");
            
            bit_node_ts=tempSet->p_head->next;
            
            ///for each ss in ks
            while(bit_node_ts!=NULL){
                unsigned int j=map_to_integer(bit_node_ts);
                for(unsigned int i=1;i<=powerSetSize;i++){
                    if( cost[i]==INFINITY) continue;
                    unsigned int unionKey =i;
                    union_bit(unionKey, j);
                    //    printf("i:%u\tj:%u\tunionKey:%u\n",i,j,unionKey);
                    if(unionKey == i && unionKey==j){
                        continue;
                    }
                    
                    double D = cost[i] + calc_minDist(obj_v->MBR, q->loc_v);
                    //     printf("cost[unionKey]:%0.3lf\tD:%0.3lf\n",cost[unionKey],D);
                    
                    if(cost[unionKey]>D){
                        //      printf("cost[unionKey]>D\n");
                        cost[unionKey]=D;
                        //group[unionkey] = group[i] \cap {e}
                        release_obj_set(group[unionKey]);
                        group[unionKey]=copy_obj_set(group[i]);
                        add_obj_set_entry(obj_v, group[unionKey]);
                        
                        //valuedSet += unionKey
                        bit_node_t* bit_node_v = alloc_bit_node();
                        bit_node_v->bit_string=unionKey;
                        bit_node_v->next=NULL;
                        bit_set_insert(valuedSet, bit_node_v);
                    }
                }
                bit_node_ts=bit_node_ts->next;
            }
            release_bit_set(tempSet);
            /*t/
             printf("markedSet:%d\n",markedSet->bit_n);
             print_bit_set(markedSet);
             printf("valuedSet:%d\n",valuedSet->bit_n);
             print_bit_set(valuedSet);
             /*t*/
            //------------------------------------
        }
        else
        {
            //printf("else....(node)\n");
            //e is an node_t*.
            //enqueue something into priority queue
            node_v = ( node_t*)e;
            
            //for each child node of e
            
            if( node_v ->level > 0)
            {
                //   printf("inserting inner node child (total number:%d)\n",node_v->num);
                
                for(int i=0; i<node_v->num; i++)
                {
                    
                    //printf("inserting inner node (i:%d)\n",i);
                    
                    //node_v is an inner-node.
                    node_t* child = (node_t*)node_v->child[i];
                    k_node_t* k_node_temp = collect_keywords_bst(child->bst_v);
                    
                    psi_t* psi_v_temp = key_intersection(q->psi_v->k_head, k_node_temp);
                    //    printf("psi_v_temp:\n");
                    //   print_k_list(psi_v_temp->k_head, stdout);
                    release_k_list(k_node_temp);
                    
                    bit_node_t* bit_node_v_temp = psi_to_bit_node(A, q->psi_v->key_n, psi_v_temp);
                    if( psi_v_temp->key_n ==0 ||bit_set_find(markedSet, bit_node_v_temp))continue;
                    
                    U->obj_arr[ rear].element = node_v->child[ i];
                    U->obj_arr[ rear].e_tag = 1;
                    U->obj_arr[ rear].key = calc_minDist_node(child, q->loc_v);
                    //Enqueue.
                    b_h_insert( U, rear++);
                    
                }//end for
                
            }
            else
            {
                // printf("inserting leaf node child (total number:%d)\n",node_v->num);
                
                for(int i=0; i<node_v->num; i++)
                {
                    
                    //   printf("inserting leaf node (i:%d)\n",i);
                    
                    //node_v is a leaf-node.
                    //the child are objs
                    obj_t* child = (obj_t*)node_v->child[i];
                    psi_t* psi_v_temp = key_intersection(q->psi_v->k_head, child->k_head);
                    //  printf("psi_v_temp:\n");
                    //  print_k_list(psi_v_temp->k_head, stdout);
                    
                    bit_node_t* bit_node_v_temp = psi_to_bit_node(A, q->psi_v->key_n, psi_v_temp);
                    
                    if( psi_v_temp->key_n ==0 ||bit_set_find(markedSet, bit_node_v_temp))continue;
                    
                    U->obj_arr[ rear].element = node_v->child[ i];
                    loc_t* loc_v = get_obj_loc(child);
                    U->obj_arr[ rear].e_tag = 2;
                    U->obj_arr[ rear].key =calc_dist_loc(loc_v, q->loc_v);
                    release_loc(loc_v);
                    
                    //Enqueue.
                    b_h_insert( U, rear++);
                    
                }//end for
            }
            //   printf("else end\n");
        }//else
        /*t/
         printf("U:\n");
         print_b_heap(stdout, U);
         printf("\n");
         /*t*/
    }//while
E:
    
    //release resource.
    release_bit_set(markedSet);
    release_bit_set(valuedSet);
    release_b_heap(U);
    /*t/
     printf("E:\n");
     print_obj_set(group[powerSetSize], stdout);
     /*t*/
    return group[powerSetSize];
}

/*
 *	The implementation of the "Cao-Appro" algorithm.
 *  CCG+15
 *
 */
obj_set_t* Sum_Appro( query_t* q)
{
	int i, rear, top;
	B_KEY_TYPE costV;
	b_heap_t* U;
	obj_set_t* V;
	void* e,* e2;
	node_t* node_v,* node_v2;
	obj_t* obj_v,* obj_v2;
    psi_t* mSet, * pSet;
 
	U = alloc_b_heap( INI_HEAP_SIZE);
    
	rear = 1;
	U->obj_arr[ rear].element = ( void*)IRTree_v.root;
	U->obj_arr[ rear].e_tag = 1;
	U->obj_arr[ rear].key = 0;
    
	b_h_insert( U, rear++);
    
    V = alloc_obj_set();
	costV = 0;

    mSet = alloc_psi();
    copy_k_list(mSet->k_head, q->psi_v->k_head);
    mSet->key_n = q->psi_v->key_n;
    
    pSet = alloc_psi();
    copy_k_list(pSet->k_head, q->psi_v->k_head);
    pSet->key_n = q->psi_v->key_n;

    //Best-first search process.
    while(mSet->key_n!=0 && !b_h_is_empty( U))
    {
        top = b_h_get_top( U);
        e = U->obj_arr[ top].element;
     
        B_KEY_TYPE key =U->obj_arr[ top].key;
        U->obj_arr[top].key =-1; //**
        
        if( U->obj_arr[ top].e_tag == 2){
            costV = costV +key;
            
            
            //e is an obj_t*.
            obj_v = ( obj_t*)e;
            
            add_obj_set_entry(obj_v, V);
            
            //pSet <- mSet
            release_psi(pSet);
            pSet = alloc_psi();
            copy_k_list(pSet->k_head, mSet->k_head);
            pSet->key_n = mSet->key_n;
            
            //mSet <- mSet - o.\psi
            psi_exclusion(mSet, obj_v);
   
            //for each element in the heap
            //  U->rear != rear
            for(int i=1;i<=U->rear;i++){
                int k = U->h_arr[i];
               // printf("============i:%d\tk:%d\t\tobj_arr[k]:%0.3lf=================\n",i,k,U->obj_arr[k].key);
                
                e2 =  U->obj_arr[k].element;
                
                if( U->obj_arr[k].e_tag == 2){//obj
                    obj_v2 = ( obj_t*)e2;
                    
                        int cnt = is_relevant_obj(obj_v2, mSet);//cnt= number of intersection of keyword in obj_v2 and mSet
                         if(cnt==0)
                            U->obj_arr[k].key = -1;
                        else
                            U->obj_arr[k].key = U->obj_arr[k].key * ((double)is_relevant_obj(obj_v2, pSet)) / ((double)cnt);
                    
                    
                }else{//node
                    node_v2 = (node_t*)e2;
               
                        int cnt = is_relevant_node(node_v2, mSet);

                        if(cnt==0)
                            U->obj_arr[k].key = -1;
                        else
                            U->obj_arr[k].key = U->obj_arr[k].key *    ((double)is_relevant_node(node_v2, pSet)) / ((double)cnt);
                    
                }
            }//end for
            b_h_restruct_heap(U, rear);
        }
        else
        {
             //e is an node_t*.
            node_v = ( node_t*)e;
            for( i=0; i<node_v->num; i++)
            {
                if( node_v->level == 0)
                {
                    //node_v is a leaf-node.
                    obj_t* child = (obj_t*)node_v->child[i];
            
                    int cnt = is_relevant_obj(child, mSet);
                   
                    if (cnt > 0){
                        loc_t* loc_v = get_obj_loc(child);
                    
                        U->obj_arr[ rear].element = node_v->child[ i];
                        U->obj_arr[ rear].e_tag = 2;
                        U->obj_arr[ rear].key = calc_dist_loc(loc_v, q->loc_v)/((double)cnt);
                    
                        //Enqueue.
                        b_h_insert( U, rear++);
                        release_loc(loc_v);
                    }
                }
                else
                {
                    //node_v is an inner-node.
                    node_t* child = (node_t*)node_v->child[i];
                    
                    int cnt = is_relevant_node(child, mSet);
                    
                    if (cnt > 0){
                        U->obj_arr[ rear].element = node_v->child[ i];
                        U->obj_arr[ rear].e_tag = 1;
                        U->obj_arr[ rear].key = calc_minDist_node(child, q->loc_v)/((double)cnt);
                        //Enqueue.
                        b_h_insert( U, rear++);
                    
                    }
                }
            }//for
        }//else
    }//while
    
	release_b_heap( U);
    release_psi(mSet);
    release_psi(pSet);
    
    //printf("costV:%f \t comp_cost:%f\n",costV, comp_cost(cost_tag, V, q));
    
	return V;
}


//===========================================================
//===========================================================
//===========================================================


bit_node_t* alloc_bit_node(){
    bit_node_t* bit_node_v;
    
    bit_node_v = ( bit_node_t*)malloc( sizeof(  bit_node_t));
	memset( bit_node_v, 0, sizeof(  bit_node_t));
    
    bit_node_v->bit_string=0;
    bit_node_v->next = NULL;
    
    /*s*/
	stat_v.memory_v += sizeof( bit_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
    
    return bit_node_v;
    
}

/*
 * return a new bit node with the same bit_string of bit_node_v
 */
bit_node_t* copy_bit_node(bit_node_t* bit_node_v){
    bit_node* bit_node_v2;
    bit_node_v2=alloc_bit_node();
    bit_node_v2->bit_string=bit_node_v->bit_string;
    bit_node_v2->next=NULL;
    return bit_node_v2;
}

bit_set_t* alloc_bit_set(){
    bit_set_t* bit_set_v;
    
    
    bit_set_v = ( bit_set_t*)malloc( sizeof(  bit_set_t));
	memset( bit_set_v, 0, sizeof(  bit_set_t));
    
    bit_set_v->bit_n=0;
    bit_set_v->p_head = alloc_bit_node();
    
    /*s*/
	stat_v.memory_v += sizeof( bit_set_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
    
    return bit_set_v;
    
}

///return true if bit_v \in bit_set_v, false otherwise
bool bit_set_find(bit_set_t* bit_set_v, bit_node_t* bit_v){
    bit_node_t* bit_v_temp;
    bit_v_temp = bit_set_v->p_head->next;
    while(bit_v_temp!=NULL){
        if(bit_v_temp->bit_string==bit_v->bit_string)
            return true;
        bit_v_temp=bit_v_temp->next;
    }
    
    return false;
}

///insert a bit_node into a bit set (after the dummy head)
void bit_set_insert(bit_set_t* bit_set_v, bit_node_t* bit_v){
    bit_node_t* bit_v_temp;
    bit_v_temp = bit_set_v->p_head->next;
    
    bit_set_v->p_head->next = bit_v;
    bit_v->next=bit_v_temp;
    
    bit_set_v->bit_n++;
    
    return;
    
}

///return true if bit_v \in bit_set_v and removed successfully, false otherwise
bool bit_set_remove(bit_set_t* bit_set_v, bit_node_t* bit_v){
    bit_node_t* bit_v1,* bit_v2;
    bit_v1 = bit_set_v->p_head->next;
    bit_v2 = bit_set_v->p_head;
    
    while(bit_v1!=NULL){
        if(bit_v1->bit_string==bit_v->bit_string){
            bit_v2->next=bit_v1->next;
            bit_v1->next=NULL;
            
            bit_set_v->bit_n--;
            return true;
        }
        bit_v1=bit_v1->next;
        bit_v2=bit_v2->next;
    }
    return false;
}


// append v2 at the tail of v1
void bit_set_append(bit_set_t* bit_set_v1, bit_set_t* bit_set_v2){
    
    bit_node_t* bit_node_v;
    bit_node_v=bit_set_v1->p_head->next;
    while(bit_node_v!=NULL){
        bit_node_v=bit_node_v->next;
    }
    
    bit_node_v->next = bit_set_v2->p_head->next;
    bit_set_v1->bit_n +=bit_set_v2->bit_n;
}


// k =[31,0]
void find_all_bit_set_sub(bit_set_t* bit_set_v, int k, BIT_TYPE bit){
    /*t/
     printf("find_all_bit_set_sub\n");
     printf("k: %d  bit: %u\n",k,bit);
     print_bit_set(bit_set_v);
     /*t*/
    
    if(k==32){
        bit_node_t* bit_node_v;
        bit_node_v=alloc_bit_node();
        bit_node_v->bit_string=bit;
        bit_node_v->next=NULL;
        bit_set_insert(bit_set_v, bit_node_v);
        return;
    }
    
    if(get_k_bit(bit, k)==0){
        //skip the bit originally is 0
        //   printf("k-th bit =0\n");
        find_all_bit_set_sub(bit_set_v,k+1,bit);
    }else{
        //k-th bit ==1
        // printf("k-th bit =  1\n");
        
        find_all_bit_set_sub(bit_set_v,k+1,bit);
        //printf("--\n");
        delete_k_bit(bit, k);///change the k-th bit to 0
        //printf("k: %d  bit: %u\n",k,bit);
        find_all_bit_set_sub(bit_set_v,k+1,bit);
    }
    
    return ;
}

bit_set_t* find_all_bit_set(bit_node_t* bit_node_v){
    /*t/
     printf("=====find_all_bit_set=====\n");
     printf("bit:\n");
     print_bit_node(bit_node_v);
     /*t*/
    
    bit_set_t* bit_set_v;
    bit_set_v = alloc_bit_set();
    
    find_all_bit_set_sub(bit_set_v,0,bit_node_v->bit_string);
    
    /*t/
     printf("bit_set_v:%d\n",bit_set_v->bit_n);
     print_bit_set(bit_set_v);
     /*t*/
    return bit_set_v;
}

/*
 * a is the mapping array
 * size is the array size
 
 */
bit_node_t* psi_to_bit_node(B_KEY_TYPE a[],int size, psi_t* psi_v){
    
    bit_node_t* bit_node_v;
    k_node_t* k_node_v;
    int i=0;
    
    
    bit_node_v=alloc_bit_node();
    k_node_v=psi_v->k_head->next;
    
    while(k_node_v!=NULL){
        i=0;
        while( i<size){
            if( a[i] ==k_node_v->key){
                //   printf("i:%d\t%.0lf\n",i,k_node_v->key);
                insert_k_bit(bit_node_v->bit_string, 31-i);
                break;
            }
            i++;
        }
        k_node_v=k_node_v->next;
    }
    
    return bit_node_v;
}




/*
 *	Release a bit_set structure.
 */
void release_bit_set( bit_set_t* bit_set_v)
{
    if (!bit_set_v) {
        return;
    }
    
	bit_node_t* bit_node_v1, *bit_node_v2;
    
	bit_node_v1 = bit_set_v->p_head;
	while( bit_node_v1->next != NULL)
	{
		bit_node_v2 = bit_node_v1->next;
		free( bit_node_v1);
		bit_node_v1 = bit_node_v2;
	}
	free( bit_node_v1);
    free(bit_set_v);
}

/*
 *	print a bit_node
 */
void print_bit_node( bit_node_t* bit_node_v)
{
    if (!bit_node_v) {
        return;
    }
    
    printf("%u",bit_node_v->bit_string);
}



/*
 *	print a bit_set one by one.
 */
void print_bit_set( bit_set_t* bit_set_v)
{
    if (!bit_set_v) {
        return;
    }
    
	bit_node_t* bit_node_v1;
    
	bit_node_v1 = bit_set_v->p_head->next;
	while( bit_node_v1 != NULL)
	{
        print_bit_node(bit_node_v1);
        printf("\t");
        bit_node_v1=bit_node_v1->next;
	}
    printf("\n");
    return;
}


/*
 map a bit string to integer
 */

unsigned int map_to_integer(bit_node_t* bit_node_v){
    
    return bit_node_v->bit_string;
    //
    //note that the bit start from left
    
    unsigned int x=0;
    
    for( int i=0;i<31;i++){
        if(get_k_bit(bit_node_v->bit_string, i)==1)
            x+=pow(2, i);
    }
    
    return x;
}




