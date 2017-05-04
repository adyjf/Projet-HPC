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


void evaluate(tree_t * T, result_t *result)
{
  node_searched++;

  move_t moves[MAX_MOVES];
  int n_moves;

  result->score = -MAX_SCORE - 1;
  result->pv_length = 0;

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

void evaluate_first(tree_t * T, result_t *result, int my_rank, int p, MPI_Status *status, MPI_Request *request, MPI_Datatype datatype)
{
  /*-----maitre-----*/
  if(my_rank==0){
    node_searched++;
    result_t result_tmp;
    int alpha_tmp;

    move_t moves[MAX_MOVES];
    int n_moves;

    result->score = -MAX_SCORE - 1;
    result->pv_length = 0;

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

    int iproc, imoves=0;
    for(iproc=1; iproc<p; iproc++){
        printf("%d\n", iproc);
        fprintf(stderr, "!!!!!! ok first 0.1 !!!!!! \n" );
        MPI_Send(&moves[iproc-1], 1, MPI_INT, iproc, TAG_DATA, MPI_COMM_WORLD);
        fprintf(stderr, "!!!!!! ok first 0.2 !!!!!! \n" );
        imoves++;
      }
    
    while (imoves!=n_moves){
      MPI_Probe(MPI_ANY_SOURCE, TAG_DATA, MPI_COMM_WORLD, status);
      MPI_Recv(&result_tmp, 1, datatype, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Recv(&alpha_tmp, 1, MPI_INT, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Send(&moves[imoves], 1, MPI_INT, status->MPI_SOURCE, TAG_DATA, MPI_COMM_WORLD);
      imoves++;
    }
    for(iproc=1; iproc<p; iproc++){
      MPI_Probe(MPI_ANY_SOURCE, TAG_DATA, MPI_COMM_WORLD, status);
      MPI_Recv(&result_tmp, 1, datatype, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
      MPI_Recv(&alpha_tmp, 1, MPI_INT, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status);
    }

    imoves=-1;
    for(iproc=1; iproc<p; iproc++){
      MPI_Send(&imoves, 1, MPI_INT, iproc, TAG_END, MPI_COMM_WORLD);
    }   
  }
  /*-----esclave-----*/
  else{
    tree_t child;
    result_t child_result;
    move_t moves;
    
    MPI_Recv(&moves, 1, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD, status);
    while(1){
      //MPI_Wait(request, status);
      if(status->MPI_TAG == TAG_DATA){
        MPI_Recv(&moves, 1, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD, status);

        play_move(T, moves, &child);
        evaluate(&child, &child_result);

        int child_score = -child_result.score;

        if (child_score > result->score) {
          result->score = child_score;
          result->best_move = moves;
          result->pv_length = child_result.pv_length + 1;
          for(int j = 0; j < child_result.pv_length; j++)
            result->PV[j+1] = child_result.PV[j];
          result->PV[0] = moves;
        }

        if (ALPHA_BETA_PRUNING && child_score >= T->beta)
          break;    

        MPI_Send(result, 1, datatype, 0, TAG_DATA, MPI_COMM_WORLD);
        T->alpha = MAX(T->alpha, child_score);
        MPI_Send(&T->alpha, 1, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD);
      }
      else if(status->MPI_TAG == TAG_END)
        break;
    }
  }

  if (TRANSPOSITION_TABLE)
    tt_store(T, result);
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

    fprintf(stderr, "ok decide 2\n" );
    if (my_rank==0){
      printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
      print_pv(T, result);
    }

    if (DEFINITIVE(result->score))
      break;

  }
}

int main(int argc, char **argv){
  double temps_calcul=0;

  /*Variables MPI*/
  int my_rank;
  int p;
  MPI_Status status;
  MPI_Request request;

  /*Initialisation MPI*/
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

  /*structure MPI_tree_t*/
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

  /*structure MPI_result_t*/
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
  fprintf(stderr, "ok1\n" );
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

  MPI_Barrier(MPI_COMM_WORLD);
  fprintf( stdout, "Processus #%d\tTemps effectif de calcul : %g sec\n", my_rank, temps_calcul);
  
  /*Désactivation MPI*/
  MPI_Type_free(&mpi_result_t);
  MPI_Type_free(&mpi_tree_t);
  MPI_Finalize();

  return 0;
}