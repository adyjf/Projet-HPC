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

//Fonction pour gérer les cas de doublons d'affichage :
// (en cas de score identique sur deux branches)
// Solution : on decide arbitrairement d'afficher le processeur au rang le plus petit
int Conflit_score_id(int score_init, int score_temp, int score_res, int iter, int nb_iter, int p, int my_rank)
{
	int rank_conflit = (iter+my_rank+1)%p; 	//rang du processeur en conflit
	
	if( (score_init == score_res) && 		//le processeur a un meilleur score identique à celui initial (avant les communications)
		(score_init == score_temp) &&		//le processeur a un meilleur score identique à celui qui arrive
		(iter != nb_iter-1) && 				//l'iteration n'est pas encore l'iteration finale (seule iteration où l'on retrouve sa propre valeur)
		(rank_conflit < my_rank) )			//le rang du processeur en conflit est plus petit que celui actuel
		return -1 ;			//renvoie -1, elimine l'affichage pour ce processeur
	else
		return 0 ;			//renvoie 0, aucune influence
	return 0;
}

//Fonction evaluate sequentielle, normale, qui sera appelée uniquement au sein d'un processeur
void evaluate(tree_t * T, result_t *result)
{
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
void evaluate_first(tree_t * T, result_t *result, int my_rank, int p, MPI_Status status, int * boss)
{
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
	int reste_n_child = n_moves%p; //calcul du reste
	int n_child = (n_moves - reste_n_child) / p; //calcul du quotient
	int n_child_min, n_child_max;
	
	//repartition automatique du travail selon le rang du processeur, SANS communication
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

	/* Procedure de debuggage */
	/*
	if(my_rank == 0) {
		printf("\nn_moves= %d.\n", n_moves);
	}
	printf("\nn_child_min = %d, n_child_max = %d depuis le proc %d.\n", n_child_min, n_child_max, my_rank);
	*/
	
	/* évalue récursivement les positions accessibles à partir d'ici */
	for (int i = n_child_min; i < n_child_max; i++) {
		tree_t child;
		result_t child_result;

		play_move(T, moves[i], &child);

		evaluate(&child, &child_result);

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
	/* tester nb impair : peut fonctionner */
	/* NB DE PROCESSEURS DOIT ETRE PAIR POUR LE MOMENT, sinon des processeurs seront bloqués */
	/* à ameliorer peut etre en reduces à terme */
	int nb_iter = (p%2 == 1) ? (p+1)/2 : p/2; //si nombre proc impair on le met à sup
	int tag = 0;
	int score_temp, score_init = result->score;
	for(int iter=0; iter<nb_iter; iter++) {
		// processeur 0
		if(my_rank == 0) {
			MPI_Send(&result->score, 1, MPI_INT, p-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2+1, nb_iter, p, my_rank);
			result->score = MAX(result->score, score_temp);
		}
		// dernier processeur soit pair soit impair
		else if((my_rank == p-1) && (my_rank%2 == 0)) {
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2+1, nb_iter, p, my_rank);
			result->score = MAX(result->score, score_temp);
		}
		else if((my_rank == p-1) && (my_rank%2 == 1)) {
			MPI_Recv(&score_temp, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2, nb_iter, p, my_rank);
			result->score = MAX(result->score, score_temp);
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
		}
		// autres processeurs pairs
		else if(my_rank%2 == 0) {
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, my_rank+1, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2+1, nb_iter, p, my_rank);
			result->score = MAX(result->score, score_temp);
		}
		// autres processeurs impairs
		else {
			MPI_Recv(&score_temp, 1, MPI_INT, my_rank+1, tag, MPI_COMM_WORLD, &status);
			*boss += Conflit_score_id(score_init, score_temp, result->score, iter*2, nb_iter, p, my_rank);
			result->score = MAX(result->score, score_temp);
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
		}
	}
	/* le processeur devient un boss (boss > 0) si il a un meilleur score, conflit de score id déjà géré */
	*boss = (result->score == score_init) ?  *boss + 1 : *boss - 1;
}

//Fonction uniquement appelée 1fois
void decide(tree_t * T, result_t *result, int my_rank, int p, MPI_Status status, int * boss)
{
	for (int depth = 1;; depth++) {
		//Chacun des processeurs possède sa propre copie de l'arbre
		T->depth = depth;
		T->height = 0;
		T->alpha_start = T->alpha = -MAX_SCORE - 1;
		T->beta = MAX_SCORE + 1;

		*boss = 0;
		
		if ((depth <= 1) && (my_rank == 0)){ //prof 1 : pas de parallelisme, rang arbitraire

			evaluate(T, result);

			//L'unique processeur affiche son score
			printf("=====================================\n");
			printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
			print_pv(T, result);

			if (DEFINITIVE(result->score)) //si score final
				break;
		} else if (depth > 1){ //prof > 1 : parallelisme 
		
			evaluate_first(T, result, my_rank, p, status, boss);

			//Seul le meilleur processeur au rang le plus petit affiche son score (conflit de score id déjà géré)
			if(*boss > 0) {
				printf("=%d====================================\n", my_rank);
				printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
				print_pv(T, result);
			}

			if (DEFINITIVE(result->score)) //si score final
				break;
		}
	}
}

int main(int argc, char **argv)
{
	/* Variables chronometre */
	double debut, fin;

	/* Variables MPI */
	int my_rank;
	int p;
	MPI_Status status;
	int boss;
	
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
	
	/* Travail et chronometrage */
	debut = my_gettimeofday();
	decide(&root, &result, my_rank, p, status, &boss);
	fin = my_gettimeofday();
	
	/* Attente de tous les processeurs */
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* Affichage temps de travail */
	fprintf( stderr, "Processus #%d\tTemps total de calcul : %g sec\n", my_rank, fin - debut);
	
	/* Attente de tous les processeurs */
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(boss > 0) {
		switch(result.score * (2*root.side - 1)) {
			case MAX_SCORE: printf("blanc gagne\n"); break;
			case CERTAIN_DRAW: printf("partie nulle\n"); break;
			case -MAX_SCORE: printf("noir gagne\n"); break;
			default: printf("BUG\n");
		}
		printf("Node searched: %llu\n", node_searched);
	}

	/* Désactivation MPI */
	MPI_Finalize();

	return 0;
}
