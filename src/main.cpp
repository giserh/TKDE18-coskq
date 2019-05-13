/*
 *	Author: Harry Kai-Ho Chan
 *	Email: khchanak@cse.ust.hk
 */

#include "sum_alg.h"
#include "unified.h"
#include "cao_alg_new.h"
#include <iostream>
#include <fstream>

using namespace std;

IRTree_t IRTree_v;
coskq_stat_t stat_v;

int cost_tag;

///
data_t* data_v;
unordered_map<KEY_TYPE, KEY_TYPE>* keyfreq_hashmap;

void coskq();

int main(int argc, char* argv[])
{
    coskq();
    return 0;
}



/*
 *	The interface for the CoSKQ query.
 *
 *  @cost_opt:
 =   1: MaxMax
 =   2: MaxMax2
 =   3: Sum
 =   4: max
 =   5: MinMax
 =   6: MinMax2
 =   7: SumMax
 =   8: SumMax2
 *
 *	@alg_opt:
 = 1: Unified-E
 = 2: Unified-A
 = 3: Cao-E
 = 4: Cao-A1
 = 5: Cao-A2
 = 6: Long-E
 = 7: Long-A
 = 8: Sum-E
 = 9: Sum-A
 = 32: Cao-E-new (in TODS15)
 = 52: Cao-A-new (in TODS15)
 = 82: Sum-E-new (in TODS15)
 *
 *  @ratio_tag:
 =1: compare the optimal/approximte ratio
 // extract the cost in result.txt to compute approximation ratio
 */
void coskq()
{
    int cntTimeOut = 0;
    int i, j, size_sum;
    B_KEY_TYPE cost, cost_sum;
    coskq_config_t* cfg;
    //coskq_stat_t* sta_v;
    query_t** q_set;
    range* MBR;
    obj_set_t *S;
    FILE* r_fp;
    //--

    memset(&stat_v, 0, sizeof(coskq_stat_t));

    size_sum = 0;
    //Read the cofig.
    ///read the config.txt, return coskq_config_t pointer
    printf("Reading configuration ...\n");
    cfg = read_config_coskq();

    //??
    cfg->prune_opt = 0;

    cost_tag = cfg->cost_measure;
    if (cost_tag == Cost_MaxMax)
        printf("The MaxMax measurement:\n");
    else if (cost_tag == Cost_MaxMax2)
        printf("The MaxMax2 measurement:\n");
    else if (cost_tag == Cost_Sum)
        printf("The Sum measurement:\n");
    else if (cost_tag == Cost_Max)
        printf("The max measurement:\n");
    else if (cost_tag == Cost_MinMax)
        printf("The MinMax measurement:\n");
	else if (cost_tag == Cost_MinMax2)
		printf("The MinMax2 measurement:\n");
	else if (cost_tag == Cost_SumMax)
        printf("The SumMax measurement:\n");
    else if (cost_tag == Cost_SumMax2)
        printf("The SumMax2 measurement:\n");
    else {
        printf("Cost measurement undefined\n");
        printf("Exiting now...\n");
        return;
    }

    //Read the data.
    ///read the "loc" and "doc" files, return data_t pointer

    printf("Reading data ...\n");
    keyfreq_hashmap = new unordered_map<KEY_TYPE, KEY_TYPE>();
    data_v = read_data_coskq(cfg, keyfreq_hashmap);
    printf("#obj:%d\n", data_v->obj_n);
    printf("#key (in hashmap):%lu\n", keyfreq_hashmap->size());

#ifndef WIN32
    float sys_t, usr_t, usr_t_sum = 0;
    struct rusage IR_tree_sta, IR_tree_end;

    GetCurTime(&IR_tree_sta);
#endif
    ///    printf("tree_tag:%d\n",cfg->tree_tag);
    //Option 1: Build the tree from scratch.
    //Build the IR-tree.
    if (cfg->tree_tag == 0) {

        printf("Building IR-tree ...\n");
        build_IRTree(data_v);
        //IR-tree is pointed by global var IRTree_v now
        print_and_check_tree(1, cfg->tree_file);
        //check_IF( );
    } else if (cfg->tree_tag == 1) {
        //Option 2: Read an existing tree.
        ///read tree file
        printf("Reading IR-Tree ...\n");
        read_tree(cfg->tree_file);
    }

#ifndef WIN32
    GetCurTime(&IR_tree_end);
    GetTime(&IR_tree_sta, &IR_tree_end, &stat_v.irtree_build_time, &sys_t);
#endif
    //Get the whole range.
    MBR = get_MBR_node(IRTree_v.root, IRTree_v.dim);

    //Generate the set of querys.
    ///within the MBR
    printf("Generating queries ...\n");

    q_set = gen_query_set2(cfg->q_set_size, cfg->q_key_n, MBR, data_v, cfg->low, cfg->high);
    if (q_set == NULL) {
        printf("Query generation failed!\n");
        exit(0);
    }

    //Query.
    printf("Performing Queries ...\n");
	if (cfg->alg_opt == 1)
		printf("Unified-Exact:\n");
	else if (cfg->alg_opt == 2)
		printf("Unified-Appro:\n");
    else if (cfg->alg_opt == 3)
        printf("Cao-Exact:\n");
    else if (cfg->alg_opt == 4)
        printf("Cao-Appro1:\n");
    else if (cfg->alg_opt == 5)
        printf("Cao-Appro2:\n");
	else  if (cfg->alg_opt == 6)
		printf("Long-Exact:\n");
	else if (cfg->alg_opt == 7)
		printf("Long-Appro:\n");
    else if (cfg->alg_opt == 8)
        printf("Cao-Sum-Exact:\n");
    else if (cfg->alg_opt == 9)
        printf("Cao-Sum-Appro:\n");
    else if (cfg->alg_opt == 32)
        printf("Cao-Exact-new:\n");
    else if (cfg->alg_opt == 52)
        printf("Cao-Appro2-new:\n");
    else if (cfg->alg_opt == 82)
        printf("Sum-Exact-new:\n");

    cost_sum = 0;
    ///for each query
    for (i = 0; i < cfg->q_set_size; i++) {
        if (cntTimeOut > 5) {
            printf("Total number of time out exceeds 5! Exit instance.\n");
            break;
        }
        printf("Query #%i ...\n", i + 1);

#ifndef WIN32
        struct rusage query_sta, query_end;

        GetCurTime(&query_sta);
#endif
	
		if (cfg->alg_opt == 1) ///Unified-Exact
		{
			S = Unified_Approach(q_set[i], 1, cost_tag, keyfreq_hashmap);//1 = Exact
		} else if (cfg->alg_opt == 2) ///Unified-Appro
		{
			S = Unified_Approach(q_set[i], 2, cost_tag, keyfreq_hashmap); //2 = Appro
		} else  if (cfg->alg_opt == 3) ///Cao-Exact
        {
            S = Cao_Exact(q_set[i], query_sta, cntTimeOut);			
        } else if (cfg->alg_opt == 4) ///Cao-Appro1
        {
            S = Cao_Appro1(q_set[i]);
        } else if (cfg->alg_opt == 5) ///Cao-Appro2
        {
            S = Cao_Appro2(q_set[i]);
        } else if (cfg->alg_opt == 6) ///Long-Exact
		{
			S = CostEnum_Exact(q_set[i], cfg->prune_opt);
		} else if (cfg->alg_opt == 7) ///Long-Appro
		{
			S = CostEnum_Appro(q_set[i]);
		} else if (cfg->alg_opt == 8) ///Cao-Sum-Exact
        {
            S = Sum_Exact(q_set[i]);
        } else if (cfg->alg_opt == 9) ///Cao-Sum-Appro
        {
            S = Sum_Appro(q_set[i]);
        } else if (cfg->alg_opt == 32) //Cao-Exact-new
        {
            S = Cao_Exact_new(cost_tag, q_set[i], keyfreq_hashmap);
        } else if (cfg->alg_opt == 52) //Cao-Appro-new
        {
            S = Cao_Appro2_new(q_set[i], keyfreq_hashmap);
        } else if (cfg->alg_opt == 82) //Sum-Exact-new
        {
            S = Sum_Exact_new(q_set[i]);
        }

        if (S == NULL) {
            printf("---- S is null ----\n");
            return;
        }

#ifndef WIN32
        GetCurTime(&query_end);

        GetTime(&query_sta, &query_end, &usr_t, &sys_t);
        usr_t_sum += usr_t;
#endif
        cost = comp_cost(cost_tag, S, q_set[i]);
        cost_sum += cost;
        size_sum += S->obj_n;
     
        //Print the accumulated query results.
        if (i == 0) {
            if ((r_fp = fopen(COSKQ_RESULT_FILE, "w")) == NULL) {
                fprintf(stderr, "Cannot open the coskq_result file.\n");
                exit(0);
            }

        } else {
            if ((r_fp = fopen(COSKQ_RESULT_FILE, "a")) == NULL) {
                fprintf(stderr, "Cannot open the coskq_result file.\n");
                exit(0);
            }
        }

        fprintf(r_fp, "Query #%i:\n", i + 1);
        fprintf(r_fp, "Keywords: ");
        print_k_list(q_set[i]->psi_v->k_head, r_fp);
        ///
        fprintf(r_fp, "Locations: ");
        print_loc(q_set[i]->loc_v, r_fp);
        ///

        //Print the query result.
        print_obj_set(S, r_fp);
        fprintf(r_fp, "Cost: %0.4lf\n", cost);
        printf("Cost: %0.4lf\n", cost);
#ifndef WIN32
        fprintf(r_fp, "Time: %f\n\n", usr_t);
        printf("Time: %f\n\n", usr_t);
#endif

        fclose(r_fp);

        //Print the statistics.
        ///it is updated after each quey performed
        stat_v.aver_cost = cost_sum / (i + 1);
        stat_v.aver_size = (float)size_sum / (float)(i + 1);

#ifndef WIN32
        stat_v.q_time = usr_t_sum / (i + 1);
#endif
        print_coskq_stat(cfg, i + 1);

        release_obj_set(S);

    } ///end for each query

    //Release the resources.
    for (i = 0; i < cfg->q_set_size; i++)
        release_query(q_set[i]);
    free(q_set);
    free(MBR);
    //printf("Memory balance: %f\n", stat_v.memory_v / cfg->q_set_size);
    free(cfg);
}
