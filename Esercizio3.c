/**
 * @brief Problem 3, Laboratory of Algorithms and Data Structures.
 * @author SCIAVICCO Guido (guido.sciavicco@unife.it)
 * @author STAN Ionel Eduard (ioneleduard.stan@unife.it)
 * @version Student
 */

/*
APPUNTI LEZIONE

	IMPLEMENTAZIONI:
		1) MINHEAP
		2) ARRAY SENZA STRUTTURA 
	
	Misurare il tempo medio al crescere della dimensione del grafo e confrontale le due implementazioni.
	Sorgente fissata.
	Densit? variabile. Per ogni densit? si crea un grafico diverso che confronta i tempi delle due implementazioni.
	
	COSTRUZIONE DEL GRAFO:
	Creo struttura di dimensione massima a 1000 e poi riempo in parte (es fino a 100)
	Per ogni coppia i-j da 1 a 100 (dimensione stabilita) decido probabilisticamente se esiste l'arco o no.
	Se arco esiste allora:
		1) Assegno all'arco un peso casuale
		2) Faccio insert casella sulla struttura. Nella lista puntata da i faccio Insert nuovo elemento di valore j e peso casuale
	I nodi sono inizializzati a infinito implementato come un numero negativo (da controllre ad ogni confronto altrimenti vincerebbe) o con un numero molto grande che non pu? essere superato.
	
	CODE DI PRIORIT? COME ARRAY:
	Array, ad ogni posizione corrisponde il nome del vertice e contiene il peso del vertice
	Decrease key = cambiare priorit? del vertice x = accedi in posizione x e assegna vuovo peso.
	Extract min =  scorri tutto array ed estrai numero posizione con peso minore.
	Su ogni casella flag che dice se esiste elemento o no (se l'ho estratto).
	
	CODE DI PRIORIT? COME MINHEAP:
	Min heap quindi ? semi ordinata su priorit? vertice (non sul nome). Ogni nodo della minheap ? una coppia.
	Decrease key = cambiare priorit? del vertice x = per evitare di scorrere tutto l`array,
	si usa un array di appoggio (P) che costruisco all'inizio e mantengo ad ogni operazione sulla heap,
	Posizioni di P corrispondono ai nomi dei vertici e ogni elemento contiene la posizione del veritne N nella heap.
	Quando chiamo ad esempio decreasekey(x), vado in posizione x di P, ottengo posizione di x nell'array della heap. H[P[x]].
	Per cambiare priorit? devo modificare il suo valore, ribilanciare MinHeap e aggiornale la sua posizione in P
	Extract min = estrai il primo elemento, ribilancia con minheapify e aggiorna array P
	Swap = prevede sia scambio su array heap che scambio su array P
	
	DIJKSTRA:
	Dijkstra lavora da una sorgente fissata (Vertice 0)
	Risultato = albero cammini minimi + pesi sui vertici per raggiungerli dalla sorgente
	Array dove la posizione corrisponde al numero del vertice. Ogni casella ? la testa della lista di adiacenza. La lista di adiacenza contiene il numero del vertice e il suo peso per raggiungerlo.
	Crea una coda per ogni implementazione.
	
	FUNZIONI ANTAGONISTE:
	1) Le due implementazioni devono dare gli stessi risultati in termini di pesi e di albero di cammini minimi
	2) Creazione manaule di un grafo con pesi e archi prestabiliti e controllo che il cammino minimo risultate sia corretto
	
	*/

// ##### LIBRARIES ##### //

#include <stdio.h>      // input-output library.
#include <stdlib.h>     // standard library.
#include <stdbool.h>    // standard boolean library.
#include <limits.h>     // limits library.
#include <time.h>       // time library.
#include <string.h>     // string library.

// ##### End of LIBRARIES ##### //

// ##### DATA STRUCTURES ##### //

// ----- MIN-HEAP ----- //

/**
 * @brief Min-heap node data structure.
 */
typedef struct min_heap_node_t {
    // Vertex number.
    int vertex_number;
    // Estimated distance from the source; needed for, e.g., extracting the minimum.
    int distance;
} min_heap_node_t;

/**
 * @brief Min-heap data structure.
 */
typedef struct min_heap_t {
    // Heap size.
    int heap_size;
    // Length of the array.
    int array_length;
    // Array of positions; needed for, e.g., decreasing the key.
    int* P;
    // Array of pointers of min-heap nodes.
    min_heap_node_t** A;
} min_heap_t;

// ----- End of MIN-HEAP ----- //

// ----- QUEUE ----- //
/**
 * @brief Queue node data structure.
 */
typedef struct queue_node_t {
    // Vertex number.
    int vertex_number;
    // Estimated distance from the source; needed for, e.g., extracting the minimum.
    int distance;
    // Flag needed for extracting the minimum.
    bool present;
} queue_node_t;

/**
 * @brief Queue data structure.
 */
typedef struct queue_t {
    // Queue size.
    int queue_size;
    // Array length.
    int array_length;
    // Array of pointers of queue nodes.
    queue_node_t** A;
} queue_t;

// ----- QUEUE ----- //

// ----- GRAPH ----- //

/**
 * @brief Adjacency list node data structure.
 */
typedef struct adj_list_node_t {
    // Target vertex number. (Note that this is called v in the text book.)
    int target;
    // Weight of the edge (u,v).
    int weight;
    // Next node.
    struct adj_list_node_t* next; 
} adj_list_node_t;

/**
 * @brief Adjacency list data structure.
 */
typedef struct adj_list_t {
    // Head of the list.
    adj_list_node_t* head;
} adj_list_t;

/**
 * @brief Graph data structure.
 */
typedef struct graph_t {
    // Number of vertices.
    int number_vertices;
    // Number of edges.
    int number_edges;
    // Adjacency list.
    adj_list_t* adj;
} graph_t;

// ----- End of GRAPH ----- //

// ----- AUXILIARY DATA STRUCTURES ----- //

/**
 * @brief Enumeration data type for the output.
 */
typedef enum output_enum_t {
    ONCONSOLE,  // On console.
    ONFILE      // On file.
} output_enum_t;

// ----- End of AUXILIARY DATA STRUCTURES ----- //

// ##### End of DATA STRUCTURES ##### //

// ##### GLOBAL VARIABLES ##### //

// Random seed (important for reproducibility).
const time_t RANDOM_SEED = 20;
// Minimum number of vertices.
const unsigned int MIN_NUM_VERTICES = 10; //10
// Maximum number of vertices.
const unsigned int MAX_NUM_VERTICES = 250; //1000
// Step from one experiment to another.
const unsigned int STEP_EXPERIMENTS = 10;
// How many experiments for a fixed number of vertices?
const unsigned int NUM_EXPERIMENTS = 50; //50 
// Source vertex number.
const unsigned int SOURCE_VERTEX_NUMBER = 0;
// Edge probability.
const unsigned int EDGE_PROBABILITY = 50; //1.0
// Maximum weight.
const unsigned int MAX_WEIGHT = 1000;
// Output type.
const output_enum_t output_type = ONCONSOLE;
// Output pointer (for printing).
FILE* output_pointer;

// ##### End of GLOBAL VARIABLES ##### //

// ##### IMPLEMENTATION OF THE FUNCTIONS ##### //

// ----- MIN-HEAP ----- //

/**
 * @brief Parent index for heap.
 * @param i Index for computing the parent index.
 * @return Parent index.
 */
int heap_parent_index(const unsigned int i) {
	 //OK - CORRETTO
    return (i-1)>>1; // Subtract 1 from i, then shift right the bit-representation by 1 position; English version of the CLRS text-book, pg. 152.
}

/**
 * @brief Left index for heap.
 * @param i Index for computing the left index.
 * @return Left index.
 */
int heap_left_index(const unsigned int i) {
	 //OK - CORRETTO
    return (i<<1)+1; // Shift left the bit-representation of i, then add 1 as lower bit; English version of the CLRS text-book, pg. 152.
}

/**
 * @brief Right index for heap.
 * @param i Index for computing the right index.
 * @return Right index.
 */
int heap_right_index(const unsigned int i) {
	 //OK - CORRETTO
    return (i<<1)+2; // Shift left the bit-representation of i, then add 2 as lower bit; English version of the CLRS text-book, pg. 152.
}

/**
 * @brief Create min-heap node.
 * @param vertex_number Vertex number.
 * @param distance Estimated distance from source.
 * @return Newly created node.
 */
min_heap_node_t* min_heap_create_node(const unsigned int vertex_number, const unsigned int distance) {
	min_heap_node_t *nodo = (min_heap_node_t *) malloc(sizeof(min_heap_node_t));
    nodo->vertex_number = vertex_number;
    nodo->distance = distance;
	return nodo;
}

/**
 * @brief Swap two min-heap nodes.
 * @param x First node.
 * @param y Second node.
 */
void min_heap_swap(min_heap_node_t** x, min_heap_node_t** y) {
    min_heap_node_t* temp = *x;
    *x = *y;
    *y = temp;
    return;
}

/**
 * @brief Min-heapify procedure.
 * @param H Heap to be heapified.
 * @param i Index needed by the procedure.
 */
void min_heap_heapify(min_heap_t* H, const unsigned int i) {
	int l = heap_left_index(i);
	int r = heap_right_index(i);
	int smallest = i;

	if(l<=H->heap_size && (H->A[l])<(H->A[i]))
	{
		smallest = l;
	}
	if(r<=H->heap_size && (H->A[r])<(H->A[smallest]))
	{
		smallest = r;
	}  
	if(smallest != i)
	{
        //aggionare p
		min_heap_swap(&H->A[i],&H->A[smallest]);
		min_heap_heapify(H,smallest);
	}
	 
    return;
}

/**
 * @brief Create empty min-heap.
 * @param array_length Array length.
 * @return Newly create min-heap.
 */
min_heap_t min_heap_create(const unsigned int array_length) {
    min_heap_t H;
    H.P = malloc(array_length * sizeof(int*));
	H.A = malloc(array_length * sizeof(min_heap_node_t**));
    H.array_length = array_length;
    H.heap_size=0;   
    return H;
}

/**
 * @brief Is min-heap empty?
 * @param H Min-heap.
 * @return true if it is.
 */
bool min_heap_is_empty(min_heap_t* H) {
	if(H->heap_size == 0)
	    return true;
	else
		return false;
}

/**
 * @brief Extract minimum from min-heap.
 * @param H Min-heap.
 * @return Min-heap node containing the minimum (estimated distance from the source).
 */
min_heap_node_t* min_heap_extract_min(min_heap_t* H) {
	if(!min_heap_is_empty(H))
	{
        int i;
		min_heap_node_t* min = H->A[0];
        min_heap_node_t* i_node;
		H->A[0] = H->A[H->heap_size];
		H->heap_size--;
        H->P[min->vertex_number] = i;
        H->P[i_node->vertex_number] = min->vertex_number;
		min_heap_heapify(H, 0);
    	return min;
	}
	exit(1);
}

/**
 * @brief Min-heap decrease-key.
 * @param H Min-heap.
 * @param vertex_number Vertex number.
 * @param distance New distance/key.
 */
void min_heap_decrease_key(min_heap_t* H, /*const*/ unsigned int vertex_number, const unsigned int distance) {
	int i;
    if(distance>(H->A[vertex_number]->distance))
	{
		exit(1);
	}
	else
	{
        i=H->P[vertex_number];
        H->A[i]->distance = distance;
		while(vertex_number > 0 && H->A[heap_parent_index(vertex_number)]->distance > H->A[vertex_number]->distance)
		{
			min_heap_swap(&H->A[i], &H->A[heap_parent_index(i)]);
			vertex_number = heap_parent_index(vertex_number);
            H->P[H->A[i]->vertex_number] = (i-1)/2;
            H->P[H->A[(i-1)/2]->vertex_number] = i;
		}
	}
    return;
}

/**
 * @brief Min-heap free.
 * @param H Min-heap.
 */
void min_heap_free(min_heap_t* H) {
	for(int i = 0; i<H->array_length; i++)
	{
		free(H->A[i]);
	}
    return;
}

// ----- End of MIN-HEAP ----- //

// ----- QUEUE ----- //

/**
 * @brief Create queue node.
 * @param vertex_number Vertex number.
 * @param distance Estimated distance from the source.
 * @return Newly create node.
 */
queue_node_t* queue_create_node(const unsigned int vertex_number, const unsigned int distance) {
	queue_node_t *nodo = (queue_node_t *) malloc(sizeof(queue_node_t));
	nodo->distance = distance;
	nodo->vertex_number = vertex_number;
	nodo->present = true;
    return nodo;
}

/**
 * @brief Create empty queue.
 * @param array_length Array length.
 * @return Newly create queue.
 */
queue_t queue_create(const unsigned int array_length) {
    queue_t Q;
    Q.array_length = array_length;
    Q.queue_size = 0;
    Q.A = (queue_node_t **) malloc(array_length * sizeof(queue_node_t**));
    return Q;
}

/**
 * @brief Is queue empty?
 * @param Q Queue.
 * @return true if it is.
 */
bool queue_is_empty(queue_t* Q) {
	if (Q->queue_size == 0)
		return true;
	else
		return false;
}

/**
 * @brief Extract minimum from queue.
 * @param Q Queue.
 * @return Queue node containing the minimum (estimated distance from the source).
 */
queue_node_t* queue_extract_min(queue_t* Q) {
	queue_node_t* min;
	int index = 0;
	int min_dis = INT_MAX;
	for(int i=0; i<Q->array_length; i++)
	{
		if(Q->A[i]->present && Q->A[i]->distance<min_dis)
		{
			min = Q->A[i];
			min_dis =  min->distance;
			index = i;
		}
	}
	if(Q->A[index]->present == true)
	{
        min = Q->A[index];
		Q->A[index]->present = false;
	}
    return min;
}

/**
 * @brief Queue decrease-key.
 * @param Q Queue.
 * @param vertex_number Vertex number.
 * @param distance New distance/key.
 */
void queue_decrease_key(queue_t* Q, const unsigned int vertex_number, const unsigned int distance) {
    if(Q->A[vertex_number]->distance > distance && (Q->A[vertex_number]->present))
    {
    	Q->A[vertex_number]->distance = distance;
	}
	return;
}

/**
 * @brief Queue free.
 * @param Q Queue.
 */
void queue_free(queue_t* Q) {
	for(int i = 0; i<Q->array_length;i++)
	{
        free(Q->A[i]);
	}
    free(Q->A);
    return;
}

// ----- End of QUEUE ----- //

// ----- GRAPH ----- //

/**
 * @brief Insert node in adjacency list.
 * @param L Adjacency list.
 * @param x Adjacency list node to be inserted.
 */
void adj_list_insert_node(adj_list_t* L, adj_list_node_t* x) {
    x->next = L->head;
    L->head = x;
   	return;

}

/**
 * @brief Create adjacency list node.
 * @param target Target vertex number.
 * @param weight Weight of the edge.
 * @return Newly created node.
 */
adj_list_node_t* adj_list_create_node(const unsigned int target, const unsigned int weight) {
	adj_list_node_t *nodo = (adj_list_node_t *) malloc(sizeof(adj_list_node_t));
    nodo->target = target;
	nodo->weight = weight;
    return nodo;
}

/**
 * @brief Add weighted edge.
 * @param G Graph.
 * @param source Source vertex number.
 * @param target Target vertex number.
 * @param weight Weight of the edge.
 */
void graph_add_edge(graph_t* G, const unsigned int source, const unsigned int target, const unsigned int weight) {
    adj_list_node_t* nodo = adj_list_create_node(target,weight);
    adj_list_insert_node(&G->adj[source], nodo);
	return;
}

/**
 * @brief Create (non-empty) graph with random edge weights.
 * @param number_vertices Number of vertices.
 * @param edge_prob Edge probability
 * @return Newly create graph.
 */
graph_t graph_create(unsigned const int number_vertices, const double edge_prob) {
    graph_t G;
    G.number_vertices = number_vertices;
    G.adj = malloc(number_vertices * sizeof(adj_list_node_t));
    int prob = 0;
    int weight;
    for(int i = 0; i<number_vertices; i++)
    {
    	G.adj[i].head = NULL;
    	for(int j = 0; j<number_vertices; j++)
    	{
    		prob = rand()%100;
				if(i!=j && (prob>edge_prob))
				{
                    weight = rand()%MAX_WEIGHT;
                    graph_add_edge(&G,i,j,weight);
				}
		}
    	
	}
    return G;
}

/**
 * @brief Free graph.
 * @param G Graph.
 */
void graph_free(graph_t* G) {
    adj_list_t list;
    adj_list_node_t* x;
    for (int i=0; i<G->number_vertices; i++) 
    {
        list = G->adj[i];
        while (list.head!=NULL)
        {
            x= list.head;
            list.head=list.head->next;
            free(x);
        }
    }
    free(G->adj);
    return;
}

// ----- End of GRAPH ----- //

// ----- ANTAGONISTIC FUNCTIONS ----- //

/**
 * @brief Print graph.
 * @param G Graph.
 */
void graph_print(graph_t* G) {

    fprintf(stdout, "G->number_vertices=%d\n", G->number_vertices);
    for (int u=0; u<G->number_vertices; u++) {
        fprintf(stdout, "adj[u=%d] ==> ", u);
        adj_list_node_t* x = G->adj[u].head;
        if (!x) printf("NULL\n");
        while (x) {
            if (x->next)
                fprintf(stdout, "(v=%d, w=%d), ", x->target, x->weight);
            else
                fprintf(stdout, "(v=%d, w=%d)\n", x->target, x->weight);
            x = x->next;
        }
    }
    return;
}

/**
 * @brief Print distances of vertices from the source.
 * @param distance Array of distances.
 * @param n Length of the array (i.e., number of vertices of the graph).
 */
void print_distances(int* distance, unsigned const int num_vertices) {
    printf("Vertex \t\t Distance\n");
    for (int u=0; u<num_vertices; u++)
        printf("%d \t\t %d\n", u, distance[u]);
}

// ----- End of ANTAGONISTIC FUNCTIONS ----- //

// ----- CORE FUNCTIONS ----- //

/**
 * @brief Dijkstra's single-source shortest-path algorithm with min-heap.
 * @param G Graph.
 * @param source Source vertex number.
 */
void dijkstra(graph_t* G, unsigned const int source) {
    return;
}

/**
 * @brief Dijkstra's single-source shortest-path algorithm with queue.
 * @param G Graph.
 * @param source Source vertex number.
 */

/*
	Array, ad ogni posizione corrisponde il nome del vertice e contiene il peso del vertice
	Decrease key = cambiare priorit? del vertice x = accedi in posizione x e assegna vuovo peso.
	Extract min =  scorri tutto array ed estrai numero posizione con peso minore.
	Su ogni casella flag che dice se esiste elemento o no (se l'ho estratto).
	
*/
void dijkstra_with_queue(graph_t* G, const unsigned int source) {
 	queue_t Q = queue_create(G->number_vertices);
    queue_node_t* q_node;
    adj_list_node_t* v_node;
    int u,v;
    int distances[G->number_vertices];
    for(int i=0;i<G->number_vertices; i++)
    {
        distances[i]=INT_MAX;
        Q.A[i]=queue_create_node(i,distances[i]);
        Q.queue_size ++;
        Q.A[i]->present = true;
    }
    queue_decrease_key(&Q,source,0); //Sorgente ha peso zero
    distances[source] = 0;
    while(!queue_is_empty(&Q))
    {
        q_node = queue_extract_min(&Q);
        Q.queue_size--;
        u = q_node->vertex_number;
        v_node = G->adj[u].head;
        while(v_node)
        {
            v = v_node->target;
            if(v_node->weight+distances[u]<distances[v])
            {
                distances[v] = distances[u]+v_node->weight;
                queue_decrease_key(&Q,v,distances[v]);
            }
            v_node=v_node->next;
        }
    }
    queue_free(&Q);
    return;
}

/**
 * @brief Polymorphic function that calls different versions of Dijkstra's algorithm.
 * @param G Graph.
 * @param priority_type Priority type.
 * @return Elapsed time in clocks.
 */
time_t do_experiment(graph_t* G, char* priority_type) {
    clock_t start_time, end_time = 0;
    
    start_time = clock();
    if (strcmp(priority_type, "min-heap") == 0) dijkstra(G, SOURCE_VERTEX_NUMBER);
    else if (strcmp(priority_type, "queue") == 0) dijkstra_with_queue(G, SOURCE_VERTEX_NUMBER);
    else {
        fprintf(stderr, "ERROR: The type of the priority can be either min-heap or queue: %s is not allowed\n", priority_type);
        exit(-1);
    }
    end_time = clock();
    
    return end_time - start_time;
}

// ----- End of CORE FUNCTIONS ----- //

// ##### End of IMPLEMENTATION OF THE FUNCTIONS ##### //

/**
 * @brief Main function.
 * @return 0 if all ok.
 */
int main() {
    // Random seed initialization.
    srand(RANDOM_SEED);
    // Elapsed time using min heaps.
    clock_t time_min_heap = 0;
    // Elapsed time using queues.
    clock_t time_queue = 0;
    
    // What is the outputPointer?
    if (output_type == ONCONSOLE || output_type == ONFILE) {
        // On console.
        if (output_type== ONCONSOLE) output_pointer = stdout;
        // On file.
        else {
            // Open file.
            output_pointer = fopen("results.txt", "w");
            // Have we opened the file?
            if (output_pointer == NULL) {
                fprintf(stderr, "ERROR: The output_pointer has not been created\n");
                exit(-1);
            }
        }
    }
    // Error.
    else {
        fprintf(stderr, "ERROR: The output_type can be only ONCONSOLE or ONFILE\n");
        exit(-1);
    }
    
    // Print the header, only if it is on console.
    if (output_type == ONCONSOLE) {
        fprintf(output_pointer, "+--------------------+---------------------+---------------------+\n");
        fprintf(output_pointer, "| Number of vertices | Min heap            | Queue               |\n");
        fprintf(output_pointer, "+--------------------+---------------------+---------------------+\n");
    }
    
    for (int num_vertices=MIN_NUM_VERTICES; num_vertices<=MAX_NUM_VERTICES; num_vertices+=STEP_EXPERIMENTS) {
        // Reset the elapsed times.
        time_min_heap = time_queue = 0;
        
        for (int experiment=0; experiment<NUM_EXPERIMENTS; experiment++) {
            // Create the graph.
            graph_t G = graph_create(num_vertices, EDGE_PROBABILITY);
            //graph_print(&G);
            // Time with min heap.
            time_min_heap += do_experiment(&G, "min-heap");
            // Time with queue.
            time_queue += do_experiment(&G, "queue");
            graph_free(&G);
        }
        
        // Printing the (sample mean as) result. Use TAB (\t) on file.
        if (output_type == ONCONSOLE)
            fprintf(output_pointer, "| %17d  | %19f | %19f |\n", 
                num_vertices, 
                (float) time_min_heap/NUM_EXPERIMENTS, 
                (float) time_queue/NUM_EXPERIMENTS);
        
        else
            fprintf(output_pointer, "%d \t%f \t%f \n",
                num_vertices,
                (float) time_min_heap/NUM_EXPERIMENTS, 
                (float) time_queue/NUM_EXPERIMENTS);
                 
    }
    
    return 0;
}


