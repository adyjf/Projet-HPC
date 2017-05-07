#include "projet.h"

/* 2017-05-04 : version 1.0 */

unsigned long long int node_searched = 0;

//Chronomètre
double my_gettimeofday()
{
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

void evaluate(tree_t * T, result_t *result){
  node_searched++;

  move_t moves[MAX_MOVES];
  int n_moves;

  result->score = -MAX_SCORE - 1;
  result->pv_length = 0;

  //evaluate_beginning(T, result, &n_moves, moves);
  if (test_draw_or_victory(T, result))
    return;

  if (TRANSPOSITION_TABLE && tt_lookup(T, result))     /* la réponse est-elle déjà connue ? */
    return;

  compute_attack_squares(T);

  /* profondeur max atteinte ? si oui, évaluation heuristique */
  if (T->depth == 0) {
    result->score = (2 * T->side - 1) * heuristic_evaluation(T);
    return;
  }

  n_moves = generate_legal_moves(T, &moves[0]);

  /* absence de coups légaux : pat ou mat */
  if (n_moves == 0) {
    result->score = check(T) ? -MAX_SCORE : CERTAIN_DRAW;
    return;
  }

  if (ALPHA_BETA_PRUNING)
  sort_moves(T, n_moves, moves);

  /* évalue récursivement les positions accessibles à partir d'ici */
  for (int i = 0; i < n_moves; i++) {
    tree_t child;
    result_t child_result;

    play_move(T, moves[i], &child);

    evaluate(&child, &child_result);

    int child_score = -child_result.score;

    if (child_score > result->score) {
      result->score = child_score;
      result->best_move = moves[i];
      result->pv_length = child_result.pv_length + 1;
      for(int j = 0; j < child_result.pv_length; j++)
        result->PV[j+1] = child_result.PV[j];
      result->PV[0] = moves[i];
    }

    if (ALPHA_BETA_PRUNING && child_score >= T->beta)
      break;    

    T->alpha = MAX(T->alpha, child_score);
  }

  if (TRANSPOSITION_TABLE)
    tt_store(T, result);
}

void deep_evaluate(tree_t *T, result_t *result, tree_t nodes[], result_t results[], int *i_nodes, int allowed_to_dig, int my_rank)
{
  node_searched++;

  move_t moves[MAX_MOVES];
  int n_moves;

  result->score = -MAX_SCORE - 1;
  result->pv_length = 0;

	//evaluate_beginning(T, result, &n_moves, moves);
  if (test_draw_or_victory(T, result))
    return;

  if (TRANSPOSITION_TABLE && tt_lookup(T, result))     /* la réponse est-elle déjà connue ? */
    return;

  compute_attack_squares(T);

  /* profondeur max atteinte ? si oui, évaluation heuristique */
  if (T->depth == 0) {
    result->score = (2 * T->side - 1) * heuristic_evaluation(T);
    return;
  }

  n_moves = generate_legal_moves(T, &moves[0]);

  /* absence de coups légaux : pat ou mat */
  if (n_moves == 0) {
    result->score = check(T) ? -MAX_SCORE : CERTAIN_DRAW;
    return;
  }

  if (ALPHA_BETA_PRUNING)
    sort_moves(T, n_moves, moves);

  /* évalue récursivement les positions accessibles à partir d'ici */
  for (int i = 0; i < n_moves; i++) {
    tree_t child;
    result_t child_result;

    play_move(T, moves[i], &child);

    if (allowed_to_dig > 0){
      fprintf(stderr, "processeur #%d deep_evaluate\tallowed_to_dig : %d\n", my_rank, allowed_to_dig);
    	deep_evaluate(&child, &child_result, nodes, results, i_nodes, allowed_to_dig-1, my_rank);
      int child_score = -child_result.score;

	    if (child_score > result->score) {
	      result->score = child_score;
	      result->best_move = moves[i];
	      result->pv_length = child_result.pv_length + 1;
	      for(int j = 0; j < child_result.pv_length; j++)
	        result->PV[j+1] = child_result.PV[j];
	      result->PV[0] = moves[i];
	    }
	    if (ALPHA_BETA_PRUNING && child_score >= T->beta)
      	break;  
      T->alpha = MAX(T->alpha, child_score); 
    }
    else{
    	nodes[(*i_nodes)] = child;
      /*for(int j=0; j<128; j++){
        (nodes[(*i_nodes)]).pieces[j] = child.pieces[j];
        (nodes[(*i_nodes)]).colors[j] = child.colors[j];
        (nodes[(*i_nodes)]).attack[j] = child.attack[j];
      }
      nodes[(*i_nodes)].side = child.side;
      nodes[(*i_nodes)].depth = child.depth;
      nodes[(*i_nodes)].height = child.height;
      nodes[(*i_nodes)].alpha = child.alpha;
      nodes[(*i_nodes)].beta = child.beta;
      nodes[(*i_nodes)].alpha_start = child.alpha_start;
      for(int j=0; j<2; j++){
        (nodes[(*i_nodes)]).king[j] = child.king[j];
        (nodes[(*i_nodes)]).pawns[j] = child.pieces[j];
      }
      nodes[(*i_nodes)].suggested_move = child.suggested_move;
      nodes[(*i_nodes)].hash = child.hash;
      for (int j=0; j<MAX_DEPTH; j++)
      {
        (nodes[(*i_nodes)]).history[j] = child.history[j];
      }*/

    	//results[(*i_nodes)] = child_result;
    	results[(*i_nodes)] = *result;
    	*i_nodes = *i_nodes+1;
      //fprintf(stderr, "deep_evaluate >> i_nodes = %d\n", *i_nodes);
    	//fprintf(stderr, "processeur #%d deep_evaluate\tallowed_to_dig : %d\n", my_rank, allowed_to_dig);
    }
  }

  if (TRANSPOSITION_TABLE)
    tt_store(T, result);

  //fprintf(stderr, "processeur #%d deep_evaluate 1\n", my_rank);
}

void evaluate_first(tree_t * T, result_t *result, int my_rank, int p, MPI_Status *status, MPI_Request *request, MPI_Datatype mpi_result_t, MPI_Datatype mpi_tree_t)
{
  //fprintf(stderr, "Ok evaluate_first processus #%d\n", my_rank);
  /*-----maitre-----*/
  if(my_rank==0){
    node_searched++; //faire un reduce à la fin dans une prochaine version pour rassembler les noeuds explores
    tree_t node_tmp;
    result_t result_tmp;
    int alpha_tmp;

    //move_t moves[MAX_MOVES];
    int n_moves;
    int n_nodes;
    int i_nodes;
    int allowed_to_dig = 0;

		tree_t nodes[MAX_NODES];
		result_t results[MAX_NODES];

    /*result->score = -MAX_SCORE - 1;
    result->pv_length = 0;

    evaluate_beginning(T, result, &n_moves, moves);*/
    //n_nodes = n_moves; //nb de move possible à la prodondeur 1, souvent insuffisant pour tous les processeurs
    
    //fprintf(stderr, "Ok evaluate_first2 processus #%d, p=%d, n_nodes=%d, T->depth=%d\n", my_rank, p, n_nodes, T->depth);
  	/* granulométrie adaptative */
    /* si nb de processeurs superieur au nb de moves on lance une evaluation de la profondeur adaptee */
		//while ((p > n_nodes) && (T->depth > 1) && (T->depth > allowed_to_dig)){
    while (T->depth > allowed_to_dig){
      //fprintf(stderr, "processeur #%d deep_evaluate 0\tallowed_to_dig : %d\n", my_rank, allowed_to_dig);
			/* remise à zeros des compteurs de nodes */
			i_nodes = 0;
			n_nodes = 0;
      deep_evaluate(T, result, nodes, results, &i_nodes, allowed_to_dig, my_rank);
      allowed_to_dig++;
			n_nodes = i_nodes; 
			//fprintf(stderr, "processeur #%d evaluate_first\tallowed_to_dig=%d, n_nodes=%d\n", my_rank, allowed_to_dig, n_nodes);
			if (p-1 <= n_nodes) 
        break;
		}

    fprintf(stderr, "N_NODES : %d\n", n_nodes);
    int iproc, imoves=0;
    for(iproc=1; iproc<p; iproc++){
      i_nodes = iproc-1;
      fprintf(stderr, "i_nodes envoye : %d\n", i_nodes);
      MPI_Send(&nodes[i_nodes], 1, mpi_tree_t, iproc, TAG_DATA, MPI_COMM_WORLD);
      MPI_Send(&results[i_nodes], 1, mpi_result_t, iproc, TAG_DATA, MPI_COMM_WORLD);
      //MPI_Send(&moves[imoves], 1, MPI_INT, iproc, TAG_DATA, MPI_COMM_WORLD);
      //MPI_Send(&moves[imoves], 1, MPI_INT, iproc, TAG_DATA, MPI_COMM_WORLD);
    }

    //while (imoves < n_moves-1){
    while(i_nodes < n_nodes-1){
      //imoves++;
      //fprintf(stderr, "processeur #%d evaluate_first while\n", my_rank);
      i_nodes++;
      fprintf(stderr, "i_nodes envoye : %d\n", i_nodes);
      MPI_Probe(MPI_ANY_SOURCE, TAG_DATA, MPI_COMM_WORLD, status);
      MPI_Recv(&node_tmp, 1, mpi_tree_t, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Recv(&result_tmp, 1, mpi_result_t, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Recv(&alpha_tmp, 1, MPI_INT, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      //MPI_Send(&moves[imoves], 1, MPI_INT, status->MPI_SOURCE, TAG_DATA, MPI_COMM_WORLD);
      MPI_Send(&nodes[i_nodes], 1, mpi_tree_t, status->MPI_SOURCE, TAG_DATA, MPI_COMM_WORLD);
      MPI_Send(&results[i_nodes], 1, mpi_result_t, status->MPI_SOURCE, TAG_DATA, MPI_COMM_WORLD);

      if(result_tmp.score > result->score){
        //result = &result_tmp;
        result->score = result_tmp.score;
        result->best_move = result_tmp.best_move;
        result->pv_length = result_tmp.pv_length;
        for(int j = 0; j < result_tmp.pv_length; j++)
          result->PV[j+1] = result_tmp.PV[j];
        printf("reactualisation tableau best_move\n");
        //T->side = node_tmp.side;
      }
      T->alpha = MAX(T->alpha, alpha_tmp);
    }
    
    for(iproc=1; iproc<p; iproc++){
      //fprintf(stderr, "processeur #%d evaluate_first dernier paquet\n", my_rank);
      MPI_Probe(MPI_ANY_SOURCE, TAG_DATA, MPI_COMM_WORLD, status);
      MPI_Recv(&node_tmp, 1, mpi_tree_t, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Recv(&result_tmp, 1, mpi_result_t, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Recv(&alpha_tmp, 1, MPI_INT, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      
      if(result_tmp.score > result->score){
        //result = &result_tmp;
        result->score = result_tmp.score;
        result->best_move = result_tmp.best_move;
        result->pv_length = result_tmp.pv_length;
        for(int j = 0; j < result_tmp.pv_length; j++)
          result->PV[j+1] = result_tmp.PV[j];
        printf("reactualisation tableau best_move\n");
        //T->side = node_tmp.side;
      }
      T->alpha = MAX(T->alpha, alpha_tmp);
    }

    //imoves=-1;
    i_nodes = 0;
    for(iproc=1; iproc<p; iproc++){
      //MPI_Send(&imoves, 1, MPI_INT, iproc, TAG_END, MPI_COMM_WORLD);
      MPI_Send(&nodes[i_nodes], 1, mpi_tree_t, iproc, TAG_END, MPI_COMM_WORLD);
      MPI_Send(&results[i_nodes], 1, mpi_result_t, iproc, TAG_END, MPI_COMM_WORLD);
      //fprintf(stderr, "processus #%d a envoye msg fin à %d\n", my_rank, iproc);
    }   
  }
  /*-----esclave-----*/
  else{
  	/*tree_t child;
    result_t child_result;*/
    tree_t node;
    result_t node_result;
    move_t move;

    node_result.score = -MAX_SCORE - 1;
    node_result.pv_length = 0;

    //MPI_Recv(&move, 1, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD, status);
    MPI_Recv(&node, 1, mpi_tree_t, 0, TAG_DATA, MPI_COMM_WORLD, status);
    MPI_Recv(result, 1, mpi_result_t, 0, TAG_DATA, MPI_COMM_WORLD, status);

    while(1){
      //MPI_Wait(request, status);
      //fprintf(stderr, "processus #%d a recu paquet %d\n", my_rank, moves);
      if(status->MPI_TAG == TAG_DATA){
        /*if (T->depth == 0) {
          result->score = (2 * T->side - 1) * heuristic_evaluation(T);
        }*/
        if(node.depth == 0){
          //node_result.score = (2 * node.side - 1) * heuristic_evaluation(&node);
          result->score = (2 * node.side - 1) * heuristic_evaluation(&node);
          fprintf(stderr, "processus #%d heuristiqueesclave\n", my_rank);
        }
        else{
          //play_move(T, move, &child);
          //evaluate(&child, &child_result);
          evaluate(&node, &node_result);

          /*int child_score = -child_result.score;
          if (child_score > result->score) {
            result->score = child_score;
            result->best_move = move;
            result->pv_length = child_result.pv_length + 1;
            for(int j = 0; j < child_result.pv_length; j++)
              result->PV[j+1] = child_result.PV[j];
            result->PV[0] = move;
          }
          if (ALPHA_BETA_PRUNING && child_score >= T->beta)
            break;
          T->alpha = MAX(T->alpha, child_score);*/

          int node_result_score = -node_result.score;
          if(node_result_score > result->score){
            result->score = node_result_score;
            result->best_move = node_result.best_move;
            result->pv_length = node_result.pv_length + 1;
            for(int j = 0; j < node_result.pv_length; j++)
              result->PV[j+1] = node_result.PV[j];
            result->PV[0] = node_result.best_move;
          }
          if (ALPHA_BETA_PRUNING && node_result_score >= T->beta)
            break;
          T->alpha = MAX(T->alpha, node_result_score);
        }
        fprintf(stderr, "processus #%d, score : %d\n", my_rank, node_result.score);
        //MPI_Send(&node_result, 1, mpi_result_t, 0, TAG_DATA, MPI_COMM_WORLD);
        MPI_Send(&node, 1, mpi_tree_t, 0, TAG_DATA, MPI_COMM_WORLD);
        MPI_Send(result, 1, mpi_result_t, 0, TAG_DATA, MPI_COMM_WORLD);
        MPI_Send(&T->alpha, 1, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD);
        MPI_Recv(&node, 1, mpi_tree_t, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status);
        MPI_Recv(result, 1, mpi_result_t, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status);
      }
      else {
        //fprintf(stderr, "processus #%d a recu msg fin\n", my_rank);
        break;
      }
    }
  }

  if (TRANSPOSITION_TABLE)
    tt_store(T, result);

  MPI_Barrier(MPI_COMM_WORLD);
  //fprintf(stderr, "processus #%d est sorti\n", my_rank);
}

void decide(tree_t * T, result_t *result, int my_rank, int p, MPI_Status *status, MPI_Request *request, MPI_Datatype mpi_result_t, MPI_Datatype mpi_tree_t, double *temps_calcul)
{
  for (int depth = 1;; depth++) {
    T->depth = depth;
    T->height = 0;
    T->alpha_start = T->alpha = -MAX_SCORE - 1;
    T->beta = MAX_SCORE + 1;

    if (my_rank==0)
      printf("=====================================\n");

  	evaluate_first(T, result, my_rank, p, status, request, mpi_result_t, mpi_tree_t);

    if (my_rank==0){
      printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
      print_pv(T, result);
    }

    MPI_Bcast(result,1,mpi_result_t,0,MPI_COMM_WORLD);

    if (DEFINITIVE(result->score))
      break;
  }
}

int main(int argc, char **argv){
  double temps_calcul=0;

  /* variables MPI */
  int my_rank;
  int p;
  MPI_Status status;
  MPI_Request request;

  /* initialisation MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  if (argc < 2) {
    if(my_rank==0)
      printf("usage: %s \"4k//4K/4P w\" (or any position in FEN)\n", argv[0]);
    exit(1);
  }

  if (my_rank==0 && ALPHA_BETA_PRUNING)
    printf("Alpha-beta pruning ENABLED\n");

  if (my_rank==0 && TRANSPOSITION_TABLE) {
    printf("Transposition table ENABLED\n");
    init_tt();
  }

  /* structure MPI_tree_t */
  int count_tree_t = 14;
  const int blocklengths_tree_t[] = {128,128,1,1,1,1,1,1,2,2,128,1,1,MAX_DEPTH};
  MPI_Aint displacements_tree_t[14];
    displacements_tree_t[0] = offsetof(tree_t, pieces);
    displacements_tree_t[1] = offsetof(tree_t, colors);
    displacements_tree_t[2] = offsetof(tree_t, side);
    displacements_tree_t[3] = offsetof(tree_t, depth);
    displacements_tree_t[4] = offsetof(tree_t, height);
    displacements_tree_t[5] = offsetof(tree_t, alpha);
    displacements_tree_t[6] = offsetof(tree_t, beta);
    displacements_tree_t[7] = offsetof(tree_t, alpha_start);
    displacements_tree_t[8] = offsetof(tree_t, king);
    displacements_tree_t[9] = offsetof(tree_t, pawns);
    displacements_tree_t[10] = offsetof(tree_t, attack);
    displacements_tree_t[11] = offsetof(tree_t, suggested_move);
    displacements_tree_t[12] = offsetof(tree_t, hash);
    displacements_tree_t[13] = offsetof(tree_t, history);
  const MPI_Datatype types_tree_t[] = {MPI_CHAR,MPI_CHAR,MPI_INT,
    MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,
    MPI_INT,MPI_INT,MPI_CHAR,
    MPI_INT,MPI_UNSIGNED_LONG_LONG,MPI_UNSIGNED_LONG_LONG};
  MPI_Datatype mpi_tree_t;
  MPI_Type_create_struct(count_tree_t, blocklengths_tree_t, displacements_tree_t, types_tree_t, &mpi_tree_t);
  MPI_Type_commit(&mpi_tree_t);

  /* structure MPI_result_t */
  int count_result_t = 4;
  const int blocklengths_result_t[4] = {1,1,1, MAX_DEPTH};
  MPI_Aint displacements_result_t[4];
    displacements_result_t[0] = offsetof(result_t, score);
    displacements_result_t[1] = offsetof(result_t, best_move);
    displacements_result_t[2] = offsetof(result_t, pv_length);
    displacements_result_t[3] = offsetof(result_t, PV);
  const MPI_Datatype types_result_t[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Datatype mpi_result_t;
  MPI_Type_create_struct(count_result_t, blocklengths_result_t, displacements_result_t, types_result_t, &mpi_result_t);
  MPI_Type_commit(&mpi_result_t);
  
  /*--------------------------------------------------*/

  tree_t root;
  result_t result;

  parse_FEN(argv[1], &root);
  if(my_rank==0)
    print_position(&root);

  //fprintf( stdout, "Processus #%d\tDEBUT\n", my_rank);
  decide(&root, &result, my_rank, p, &status, &request, mpi_result_t, mpi_tree_t, &temps_calcul);
  //fprintf(stderr, "ok1\n" );
  if (my_rank==0){
    printf("\nDécision de la position: ");
    switch(result.score * (2*root.side - 1)) {
      case MAX_SCORE: printf("blanc gagne\n"); break;
      case CERTAIN_DRAW: printf("partie nulle\n"); break;
      case -MAX_SCORE: printf("noir gagne\n"); break;
      default: printf("BUG\n");
    }
    printf("Node searched: %llu\n", node_searched);
  }

  if (TRANSPOSITION_TABLE)
    free_tt();

  fprintf( stdout, "Processus #%d\tTemps effectif de calcul : %g sec\n", my_rank, temps_calcul);
  
  /*Désactivation MPI*/
  MPI_Type_free(&mpi_result_t);
  MPI_Type_free(&mpi_tree_t);
  MPI_Finalize();

  return 0;
}