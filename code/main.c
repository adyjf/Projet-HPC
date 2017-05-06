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

void evaluate_beginning(tree_t * T, result_t *result, int * n_moves, move_t moves[]){
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

  *n_moves = generate_legal_moves(T, &moves[0]);
  
  /* absence de coups légaux : pat ou mat */
  if (n_moves == 0) {
    result->score = check(T) ? -MAX_SCORE : CERTAIN_DRAW;
    return;
  }

  if (ALPHA_BETA_PRUNING)
    sort_moves(T, *n_moves, moves);
}

void evaluate(tree_t * T, result_t *result){
  node_searched++;

  move_t moves[MAX_MOVES];
  int n_moves;

  result->score = -MAX_SCORE - 1;
  result->pv_length = 0;

  evaluate_beginning(T, result, &n_moves, moves);

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

	evaluate_beginning(T, result, &n_moves, moves);

  /* évalue récursivement les positions accessibles à partir d'ici */
  for (int i = 0; i < n_moves; i++) {
    tree_t child;
    result_t child_result;

    play_move(T, moves[i], &child);

    if (allowed_to_dig > 0){
      //fprintf(stderr, "processeur #%d deep_evaluate 0\tallowed_to_dig : %d\n", my_rank, allowed_to_dig);
    	deep_evaluate(&child, &child_result, nodes, results, i_nodes, allowed_to_dig-1, my_rank);
      //fprintf(stderr, "processeur #%d deep_evaluate 0\tallowed_to_dig : %d\n", my_rank, allowed_to_dig);
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
    	results[(*i_nodes)] = child_result;
    	*i_nodes++;
    }
  }

  if (TRANSPOSITION_TABLE)
    tt_store(T, result);

  //fprintf(stderr, "processeur #%d deep_evaluate 1\n", my_rank);
}

void evaluate_first(tree_t * T, result_t *result, int my_rank, int p, MPI_Status *status, MPI_Request *request, MPI_Datatype datatype)
{
  /*-----maitre-----*/
  if(my_rank==0){
    node_searched++; //faire un reduce à la fin dans une prochaine version pour rassembler les noeuds explores
    result_t result_tmp;
    int alpha_tmp;

    move_t moves[MAX_MOVES];
    int n_moves;
    int n_nodes;
    int i_nodes;
    int allowed_to_dig = 0;

		tree_t nodes[MAX_NODES];
		result_t results[MAX_NODES];

    result->score = -MAX_SCORE - 1;
    result->pv_length = 0;

    evaluate_beginning(T, result, &n_moves, moves);
    n_nodes = n_moves; //nb de move possible à la prodondeur 1, souvent insuffisant pour tous les processeurs
    
    /* granulométrie adaptative */
    /* si nb de processeurs superieur au nb de moves on lance une evaluation de la profondeur adaptee */
		while ((p > n_nodes) && (T->depth > 1) && (T->depth > allowed_to_dig)){
      fprintf(stderr, "processeur #%d deep_evaluate 0\tallowed_to_dig : %d\n", my_rank, allowed_to_dig);
			/* remise à zeros des compteurs de nodes */
			allowed_to_dig++;
			i_nodes = 0;
			n_nodes = 0;
      //fprintf(stderr, "processeur #%d evaluate first 1\n", my_rank);
			deep_evaluate(T, result, nodes, results, &i_nodes, allowed_to_dig, my_rank);
      //fprintf(stderr, "processeur #%d evaluate first 2\n", my_rank);
			n_nodes = i_nodes + 1; 
		}
  	
    int iproc, imoves=0;
    for(iproc=1; iproc<p; iproc++){
      imoves = iproc-1;
      MPI_Send(&moves[imoves], 1, MPI_INT, iproc, TAG_DATA, MPI_COMM_WORLD);
      //fprintf(stderr, "node %d envoye a %d\nimoves = %d\n", iproc-1, iproc, imoves);
    }

    while (imoves < n_moves-1){
      imoves++;
      MPI_Probe(MPI_ANY_SOURCE, TAG_DATA, MPI_COMM_WORLD, status);
      MPI_Recv(&result_tmp, 1, datatype, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      //fprintf(stderr, "processus #%d a bien recu de la part de #%d\n", my_rank, status->MPI_SOURCE);
      MPI_Recv(&alpha_tmp, 1, MPI_INT, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Send(&moves[imoves], 1, MPI_INT, status->MPI_SOURCE, TAG_DATA, MPI_COMM_WORLD);
      
      if(result_tmp.score > result->score){
        //result = &result_tmp;
        result->score = result_tmp.score;
        result->best_move = result_tmp.best_move;
        result->pv_length = result_tmp.pv_length;
        for(int j = 0; j < result_tmp.pv_length; j++)
          result->PV[j+1] = result_tmp.PV[j];
        //printf("result score : %d\n", result_tmp.score);
      }
      T->alpha = MAX(T->alpha, alpha_tmp);
    }

    for(iproc=1; iproc<p; iproc++){
      MPI_Probe(MPI_ANY_SOURCE, TAG_DATA, MPI_COMM_WORLD, status);
      MPI_Recv(&result_tmp, 1, datatype, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Recv(&alpha_tmp, 1, MPI_INT, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      
      if(result_tmp.score > result->score){
        //result = &result_tmp;
        result->score = result_tmp.score;
        result->best_move = result_tmp.best_move;
        result->pv_length = result_tmp.pv_length;
        for(int j = 0; j < result_tmp.pv_length; j++)
          result->PV[j+1] = result_tmp.PV[j];
        //printf("result score : %d\n", result_tmp.score);
      }
      T->alpha = MAX(T->alpha, alpha_tmp);
    }

    imoves=-1;
    for(iproc=1; iproc<p; iproc++){
      MPI_Send(&imoves, 1, MPI_INT, iproc, TAG_END, MPI_COMM_WORLD);
      //fprintf(stderr, "processus #%d a envoye msg fin à %d\n", my_rank, iproc);
    }   
  }
  /*-----esclave-----*/
  else{
    tree_t child;
    result_t child_result;
    move_t move;

    result->score = -MAX_SCORE - 1;
    result->pv_length = 0;

    MPI_Recv(&move, 1, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD, status);
    //fprintf(stderr, "processus #%d a bien recu paquet %d\n", my_rank, moves);
    while(1){
      //MPI_Wait(request, status);
      //fprintf(stderr, "processus #%d a recu paquet %d\n", my_rank, moves);
      if(status->MPI_TAG == TAG_DATA){
        //fprintf(stderr, "processus #%d a recu paquet %d\n", my_rank, moves);
        if (T->depth == 0) {
          result->score = (2 * T->side - 1) * heuristic_evaluation(T);
        }
        else{
          play_move(T, move, &child);
          evaluate(&child, &child_result);

          int child_score = -child_result.score;

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
          T->alpha = MAX(T->alpha, child_score);
        }

        MPI_Send(result, 1, datatype, 0, TAG_DATA, MPI_COMM_WORLD);
        MPI_Send(&T->alpha, 1, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD);
        MPI_Recv(&move, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status);
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

void decide(tree_t * T, result_t *result, int my_rank, int p, MPI_Status *status, MPI_Request *request, MPI_Datatype datatype, double *temps_calcul)
{
  for (int depth = 1;; depth++) {
    T->depth = depth;
    T->height = 0;
    T->alpha_start = T->alpha = -MAX_SCORE - 1;
    T->beta = MAX_SCORE + 1;

    if (my_rank==0)
      printf("=====================================\n");

    evaluate_first(T, result, my_rank, p, status, request, datatype);
    //fprintf(stderr, "ok decide processus #%d\n", my_rank);

    if (my_rank==0){
      printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
      print_pv(T, result);
    }

    MPI_Bcast(result,1,datatype,0,MPI_COMM_WORLD);

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

  decide(&root, &result, my_rank, p, &status, &request, mpi_result_t, &temps_calcul);
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