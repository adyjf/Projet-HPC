#include "projet.h"

/* 2017-02-23 : version 1.0 */

unsigned long long int node_searched = 0;

//Chronomètre
double my_gettimeofday(){
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

//Fonction qui sera appelée uniquement au sein d'un processeur
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
	n_child = (reste_n_child == 0) ? n_child : n_child + 1; //si le reste est non-nul, on ajoute +1 au quotient

	/* évalue récursivement les positions accessibles à partir d'ici */
	int n_child_min = (my_rank)*n_child; //borne inf incluse
	int n_child_max = MIN(n_moves, (my_rank+1)*n_child); //borne sup non incluse
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

	/* transmission des résultats */

	/* en anneau, NB DE PROCESSEURS DOIT ETRE PAIR sinon des processeurs seront bloqués */
	/* à ameliorer en reduce à terme */
	/* fonctionne actuellement uniquement pour le score du coup, à terme tout remplir pour le reste des données */
	int nb_iter = (p%2 == 1) ? (p+1)/2 : p/2; //si nombre proc impair on le met à sup
	int tag = 0;
	int score_temp;
	int score_initial = result->score;
	for(int i=0; i<nb_iter; i++) {
		// processeur 0
		if(my_rank == 0) {
			MPI_Send(&result->score, 1, MPI_INT, p-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &status);
			result->score = MAX(result->score, score_temp);
		}
		// dernier processeur soit pair soit impair
		else if((my_rank == p-1) && (my_rank%2 == 0)) {
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			result->score = MAX(result->score, score_temp);
		}
		else if((my_rank == p-1) && (my_rank%2 == 1)) {
			MPI_Recv(&score_temp, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
			result->score = MAX(result->score, score_temp);
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
		}
		// autres processeurs pairs
		else if(my_rank%2 == 0) {
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
			MPI_Recv(&score_temp, 1, MPI_INT, my_rank+1, tag, MPI_COMM_WORLD, &status);
			result->score = MAX(result->score, score_temp);
		}
		// autres processeurs impairs
		else {
			MPI_Recv(&score_temp, 1, MPI_INT, my_rank+1, tag, MPI_COMM_WORLD, &status);
			result->score = MAX(result->score, score_temp);
			MPI_Send(&result->score, 1, MPI_INT, my_rank-1, tag, MPI_COMM_WORLD);
		}
	}
	/* le processeur devient le boss (boss = 1) si il est le meilleur */
	*boss = (result->score == score_initial) ? 1 : 0;
}

//Fonction uniquement appelée 1fois
void decide(tree_t * T, result_t *result, int my_rank, int p, MPI_Status status)
{
	//NEW
	int boss;
	
	for (int depth = 1;; depth++) {
		//Chacun des processeurs possède sa propre copie de l'arbre
		T->depth = depth;
		T->height = 0;
		T->alpha_start = T->alpha = -MAX_SCORE - 1;
		T->beta = MAX_SCORE + 1;

		//NEW
		boss = 0;

		if(my_rank == 0) { //maitre
		printf("=====================================\n");
		}

		//Tous les processeurs lancent cette fonction
		//evaluate_first(T, result, my_rank, p, status);
		evaluate_first(T, result, my_rank, p, status, &boss);

		//if(my_rank == 0) { //maitre
		if(boss == 1) {
			printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
			print_pv(T, result);
		}

		//Communication peut etre à voir
		if (DEFINITIVE(result->score))
			break;
	}
}

int main(int argc, char **argv)
{
	/*Variables Chrono*/
	double debut, fin;

	/* Variables MPI*/
	int my_rank;
	int p;
	MPI_Status status;
	
	/* Initialisation MPI*/
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	if(my_rank == 0) { //maitre
		if (argc < 2) {
		printf("usage: %s \"4k//4K/4P w\" (or any position in FEN)\n", argv[0]);
		exit(1);
		}
	}

	tree_t root;
	result_t result;
	parse_FEN(argv[1], &root);
		
	if(my_rank == 0) { //maitre
		print_position(&root);
	}

	/* debut du chronometrage */
	debut = my_gettimeofday();

	decide(&root, &result, my_rank, p, status);

	/* fin du chronometrage */
	fin = my_gettimeofday();
	fprintf( stderr, "Processus #%d\tTemps total de calcul : %g sec\n", my_rank, fin - debut);

	if(my_rank == 0) { //maitre
		printf("\nDécision de la position: ");
		switch(result.score * (2*root.side - 1)) {
			case MAX_SCORE: printf("blanc gagne\n"); break;
			case CERTAIN_DRAW: printf("partie nulle\n"); break;
			case -MAX_SCORE: printf("noir gagne\n"); break;
			default: printf("BUG\n");
		}

		printf("Node searched: %llu\n", node_searched);
	}

	/*Désactivation MPI*/
	MPI_Finalize();

	return 0;
}
