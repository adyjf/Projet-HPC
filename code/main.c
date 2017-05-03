#include "projet.h"

/* 2017-02-23 : version 1.0 */

unsigned long long int node_searched = 0;

//Chronomètre
double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

//Fonction pour gérer les cas de doublons d'affichage SANS communication supplémentaire:
// (en cas de score identique sur deux branches)
// Solution : on decide arbitrairement d'afficher le processeur au rang le plus petit
int Conflit_score_id(int score_init, int score_temp, int score_res, int iter, int nb_iter, int p, int my_rank)
{
	int rank_conflit = (iter+my_rank+1)%p; 	//rang du processeur en conflit
	
	if( (score_init == score_res) && 		//le processeur a un meilleur score identique à celui initial (avant les communications)
		(score_init == score_temp) &&			//le processeur a un meilleur score identique à celui qui arrive
		(iter != nb_iter-1) && 						//l'iteration n'est pas encore l'iteration finale (seule iteration où l'on retrouve sa propre valeur)
		(rank_conflit < my_rank) ) {			//le rang du processeur en conflit est plus petit que celui actuel
		return -100 ;		//renvoie -100, elimine l'affichage definitivement pour cette profondeur et pour ce processeur
	}			
	else {
		return 0 ;			//renvoie 0, aucune influence
	}	
	return 0;
}

void evaluate(tree_t * T, result_t *result, int my_rank, int p, MPI_Status status, int * boss, double * temps_calcul)
//Fonction qui doit repartir le travail avant de travailler.
{
	/* structure MPI pour result_t */
	const int count_items = 4;
	const int blocklengths_array[4] = {1,1,1, MAX_DEPTH};
	MPI_Aint offsets_array[4];
	const MPI_Datatype types_array[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype MPI_result_t;
	
	offsets_array[0] = offsetof(result_t, score);
	offsets_array[1] = offsetof(result_t, best_move);
	offsets_array[2] = offsetof(result_t, pv_length);
	offsets_array[3] = offsetof(result_t, PV);

	MPI_Type_create_struct(count_items, blocklengths_array, offsets_array, types_array, &MPI_result_t);
	MPI_Type_commit(&MPI_result_t);

	/* suite de la fonction */
	node_searched++;

	move_t moves[MAX_MOVES];
	int n_moves;

	result->score = -MAX_SCORE - 1;
	result->pv_length = 0;

	if (test_draw_or_victory(T, result))
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

	/* division des taches */
	int reste_n_child = n_moves%p; 					//calcul du reste
	int n_child = (n_moves - reste_n_child) / p; 	//calcul du quotient

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

		T->alpha = MAX(T->alpha, child_score);
	}
}

//Fonction evaluate PARALLELE qui sera appelée uniquement à la profondeur 1, pour dispatcher le travail des processeurs
void evaluate_first(tree_t * T, result_t *result, int my_rank, int p, MPI_Status status, int * boss, double * temps_calcul)
{
	/* Variables chronometre */
	double debut, fin;
	debut = my_gettimeofday();

	int i;
	node_searched++;

	move_t moves[MAX_MOVES];
	int n_moves;

	result->score = -MAX_SCORE - 1;
	result->pv_length = 0;

	if (test_draw_or_victory(T, result))
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

	/* division des taches */
	int reste_n_child = n_moves%p; 								//calcul du reste
	int n_child = (n_moves - reste_n_child) / p; 	//calcul du quotient
	int n_child_min, n_child_max;									//borne min incluse, borne max non incluse
	
	/* repartition automatique du travail selon le rang du processeur SANS communication */
	if ((reste_n_child!=0) && (my_rank<reste_n_child)) { //traite aussi le cas n_child = 0
		n_child_min = (my_rank)*n_child + my_rank; 
		n_child_max = (my_rank+1)*n_child + my_rank + 1;
	} else if (n_child==0){
		n_child_min = 0; 
		n_child_max = -1;
	} else {
		n_child_min = (my_rank)*n_child + reste_n_child;
		n_child_max = (my_rank+1)*n_child + reste_n_child;
	}

	/* évalue récursivement les positions accessibles à partir d'ici */
	for (int i = n_child_min; i < n_child_max; i++) {
		tree_t child;
		result_t child_result;

		play_move(T, moves[i], &child);

		evaluate(&child, &child_result); //on appelle bien la fonction evaluate initiale

		int child_score = -child_result.score;

		if (child_score > result->score) {
			result->score = child_score; //seul rappatrié dans le final
			result->best_move = moves[i]; 
			result->pv_length = child_result.pv_length + 1;
			for(int j = 0; j < child_result.pv_length; j++)
				result->PV[j+1] = child_result.PV[j];
			result->PV[0] = moves[i];
		}

		T->alpha = MAX(T->alpha, child_score);
	}

	/* transmission des résultats (communications en anneau) */
	/* nb de processeurs peut être pair ou impair */
	int nb_iter = (p%2 == 1) ? (p+1)/2 : p/2; //si nombre proc impair on le met à sup
	int tag = 0;
	int score_temp, score_init = result->score;
	
	fin = my_gettimeofday();
	*temps_calcul += fin-debut;

	/* Attente de tous les processeurs */
  MPI_Barrier(MPI_COMM_WORLD); //AJOUT
	
	debut = my_gettimeofday();

	for(int iter=0; iter<nb_iter; iter++) {
		// processeur 0
		if(my_rank == 0) {
			MPI_Send(&result->score, 1, MPI_INT, p-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2+1, nb_iter*2, p, my_rank);
			result->score = MAX(result->score, score_temp);
		}
		// dernier processeur soit pair soit impair
		else if((my_rank == p-1) && (my_rank%2 == 0)) {
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2+1, nb_iter*2, p, my_rank);
			result->score = MAX(result->score, score_temp);
		}
		else if((my_rank == p-1) && (my_rank%2 == 1)) {
			MPI_Recv(&score_temp, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2, nb_iter*2, p, my_rank);
			result->score = MAX(result->score, score_temp);
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
		}
		// autres processeurs pairs
		else if(my_rank%2 == 0) {
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, my_rank+1, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2+1, nb_iter*2, p, my_rank);
			result->score = MAX(result->score, score_temp);
		}
		// autres processeurs impairs
		else {
			MPI_Recv(&score_temp, 1, MPI_INT, my_rank+1, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2, nb_iter*2, p, my_rank);
			result->score = MAX(result->score, score_temp);
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
		}
	}
	/* le processeur devient un boss (boss > 0) si il a un meilleur score, conflit de score id déjà géré */
	*boss = (result->score == score_init) ?  *boss + 1 : *boss - 1;

	fin = my_gettimeofday();
	*temps_calcul += fin-debut;

}

//Fonction uniquement appelée 1fois
void decide(tree_t * T, result_t *result, int my_rank, int p, MPI_Status status, int * boss, double * temps_calcul)
{
	/* Variables chronometre */
	double debut, fin;

	for (int depth = 1;; depth++) {
		//Chacun des processeurs possède sa propre copie de l'arbre
		T->depth = depth;
		T->height = 0;
		T->alpha_start = T->alpha = -MAX_SCORE - 1;
		T->beta = MAX_SCORE + 1;

		*boss = 0; //variable remise à zero sur chaque processeur, à chaque profondeur
		
		if ((depth <= 5) && (my_rank == 0)){ //prof 1 : pas de parallelisme, rang arbitraire
			debut = my_gettimeofday();

			evaluate(T, result);

			//L'unique processeur affiche son score
			printf("=====================================\n");
			printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
			print_pv(T, result);

			fin = my_gettimeofday();
			*temps_calcul += fin-debut;

			if (DEFINITIVE(result->score)) //si score final
				break;
		} else if (depth > 5){ //prof > 1 : parallelisme 
		
			evaluate_first(T, result, my_rank, p, status, boss, temps_calcul);

			//Seul le meilleur processeur au rang le plus petit affiche son score (conflit de score id déjà géré)
			if(*boss > 0) {
				printf("=====================================\n");
				printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
				print_pv(T, result);
			}

			if (DEFINITIVE(result->score)) //si score final
				break;
		}

		/* Attente de tous les processeurs */
  	MPI_Barrier(MPI_COMM_WORLD); //AJOUT

	}
}

int main(int argc, char **argv)
{
	/* Variable chronometre */
	double temps_calcul;

	/* Variables MPI */
	int my_rank;
	int p;
	MPI_Status status;
	int boss; //variable qui montre si le processeur possede ou non le meilleur coup
	
	/* Initialisation MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	/* Gestion erreur */
	if (argc < 2) {
		if(my_rank == 0) printf("usage: %s \"4k//4K/4P w\" (or any position in FEN)\n", argv[0]);
		exit(1);
	}
	
	tree_t root;
	result_t result;
	parse_FEN(argv[1], &root);

	if (my_rank == 0) print_position(&root);
	
	/* Attente de tous les processeurs */
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);
        
	/* Travail */
	decide(&root, &result, my_rank, p, status, &boss, &temps_calcul);
	
	/* Attente de tous les processeurs */
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);
	
	/* Affichage du meilleur score uniquement par le processeur concerné */
	if(boss > 0) {
		switch(result.score * (2*root.side - 1)) {
			case MAX_SCORE: printf("blanc gagne\n"); break;
			case CERTAIN_DRAW: printf("partie nulle\n"); break;
			case -MAX_SCORE: printf("noir gagne\n"); break;
			default: printf("BUG\n");
		}
		printf("Node searched: %llu\n", node_searched);
	}

	/* Attente de tous les processeurs */
  MPI_Barrier(MPI_COMM_WORLD);
  sleep(1);

	/* Affichage temps de travail */
  fprintf( stdout, "Processus #%d\tTemps effectif de calcul : %g sec\n", my_rank, temps_calcul);

	/* Désactivation MPI */
	MPI_Finalize();

	return 0;
}
