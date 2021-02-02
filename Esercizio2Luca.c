/**
 * @brief Problem 2, Laboratory of Algorithms and Data Structures.
 * @author SCIAVICCO Guido (guido.sciavicco@unife.it)
 * @author STAN Ionel Eduard (ioneleduard.stan@unife.it)
 * @author MARCHIORI Luca (luca01.marchiori@edu.unife.it)
 * @version Student
 */

// ##### LIBRARIES ##### //

// Standard input-output library (e.g., fprintf).
#include <stdio.h>
// Time library (e.g., time, clock()).
#include <time.h>
// Standard library (e.g., rand, srand).
#include <stdlib.h>
// Boolean library (e.g., bool).
#include <stdbool.h>
// String library (e.g., strcmp)
#include <string.h>

// ##### End of LIBRARIES ##### //

// ##### DATA STRUCTURES ##### //

// ----- LINKED LIST ----- //

/**
 * @brief Linked list node data type.
 */
typedef struct linkedListNode_t {
    // Value contained in the node.
    int value;
    // Pointer to next node.
    struct linkedListNode_t* next;
    // Pointer to previous node.
    struct linkedListNode_t* prev;
} linkedListNode_t;

/**
 * @brief Linked list data type.
 */
typedef struct linkedList_t {
    // Size in number of nodes of the list.
    unsigned int size;
    // Pointer to the head node of the list.
    struct linkedListNode_t* head;
} linkedList_t;

// ----- End of LINKED LIST ----- //

// ----- HASHTABLE ----- //

/**
 * @brief Hashtable entry data type.
 */
typedef struct hashtableEntry_t {
    // Pointer to the list.
    struct linkedList_t* list;
} hashtableEntry_t;

/**
 * @brief Hashtable data type.
 */
typedef struct hashtable_t {
    // Size in number of entries of the hashtable.
    unsigned int size;
    // Array of pointers to entries.
    struct hashtableEntry_t** entry;
} hashtable_t;

// ----- End of HASHTABLE ----- //

// ----- RED BLACK TREE (RBT) ----- //

/**
 * @brief RBT node data type.
 */
typedef struct rbtNode_t {
    // Value contained in the node.
    int value;
    // Color of the node.
    char color;
    // Pointer to the parent node.
    struct rbtNode_t* parent;
    // Pointer to the left node.
    struct rbtNode_t* left;
    // Pointer to the right node.
    struct rbtNode_t* right;
} rbtNode_t;

/**
 * @brief RBT data type.
 */
typedef struct rbt_t {
    // Size in number of nodes of the RBT.
    unsigned int size;
    // Pointer to the root node.
    struct rbtNode_t* root;
    // Pointer to RBT NIL node.
    struct rbtNode_t* nil;
} rbt_t;

/**
 * @brief RBT test data structure.
 */
typedef struct rbtTestStructure_t {
    // Array that contains the in order visit values of the RBT.
    int* A;
    // Current index of the array.
    int index;
} rbtTestStructure_t;

// ----- End of RBT ----- //

// ----- AUXILIARY DATA STRUCTURES ----- //

/**
 * @brief Enumeration data type for the output.
 */
typedef enum outputEnum_t {
    ONCONSOLE,  // On console.
    ONFILE      // On file.
} outputEnum_t;

// ----- End of AUXILIARY DATA STRUCTURES ----- //

// ##### End of DATA STRUCTURES ##### //

// ##### GLOBAL VARIABLES ###### //

// Random seed (important for reproducibility).
time_t RANDOM_SEED = 20;
// Maximum random number allowed.
const unsigned int MAX_RANDOM_NUMBER = 100;
// Minimum number of operations.
const unsigned int MIN_OPERATIONS = 100;
// Maximum number of operations.
const unsigned int MAX_OPERATIONS = 2000;
// Step of the experiment.
const unsigned int STEP = 100;
// Number of experiments.
const unsigned int NUM_EXPERIMENTS = 50;
// Percentage of insert operations.
const unsigned int PERCENTAGE_INSERTIONS = 40;
// Size of the hashtable.
const unsigned int NUM_ENTRIES = 3;    
// Test data structures?
const bool TEST_DATA_STRUCTURES = true;
// Number of elements for testing.
const unsigned int NUM_ELEMENTS_FOR_TEST = 1000;
// Output type.
const outputEnum_t outputType = ONCONSOLE;
//conteggio
unsigned int CONT = 0;
//risultato
bool R = false;
// Output pointer (for printing).
FILE* outputPointer;

// ##### End of GLOBAL VARIABLES #####

// ##### PROTOTYPES OF THE FUNCTIONS ##### //

// ----- LINKED LIST ----- //

/**
 * @brief Create a new linked list node.
 * @param Value that the linked list node should contain.
 * @return Created linked list node.
 */
linkedListNode_t* createLinkedListNode(const int);

/**
 * @brief Create a new linked list.
 * @return Created linked list.
 */
linkedList_t* createLinkedList();

/**
 * @brief Insert linked list node in the head of the linked list.
 * @param The linked list.
 * @param Linked list node to be inserted.
 */
void linkedListInsert(linkedList_t* head, linkedListNode_t* nodo);

/**
 * @brief Search for a value in the linked list.
 * @param The linked list.
 * @param Value to be searched for.
 * @return First linked list node containing such value, if it exists; otherwise, NULL.
*/
linkedListNode_t* linkedListSearch(linkedList_t*, const int);

/**
 * @brief Delete a linked list node from linked list.
 * @param The linked list.
 * @param The linked list node to be deleted.
 */
void linkedListDelete(linkedList_t*, linkedListNode_t*);

/**
 * @brief Print the linked list.
 * @param Linked list to be printed.
 */
void linkedListPrint(linkedList_t*);

/**
 * @brief Free the linked list.
 * @param Linked list to be freed.
 */
void linkedListFree(linkedList_t*);

// ----- End of LINKED LIST ----- //

// ----- HASHTABLE ----- //

/**
 * @brief Create a new hashtable.
 * @param The size of the hashtable (i.e., the number of entries).
 * @return The created hashtable.
 */
hashtable_t* createHashtable(const unsigned int);

/**
 * @brief Hash function computing the key for a given integer.
 * @param The hashtable needed to access the size of it.
 * @param The integer for which the key must be computed.
 * @return The computed key.
 */
const unsigned int hashFunction(hashtable_t*, const int);

/**
 * @brief Insert value in the hashtable.
 * @param The hashtable.
 * @param Value to be inserted.
 */
void hashtableInsert(hashtable_t*, const int);

/**
 * @brief Search for a value in the hashtable.
 * @param The hashtable.
 * @param Value to be searched.
 * @return Linked list node containing such value, if it exists; otherwise, NULL.
 */
linkedListNode_t* hashtableSearch(hashtable_t*, const int);

/**
 * @brief Delete value from hashtable.
 * @param The hashtable.
 * @param Linked list node to be deleted.
 */
void hashtableDelete(hashtable_t*, linkedListNode_t*);

/**
 * @brief Print the hashtable.
 * @param Hashtable to be printed.
 */
void hashtablePrint(hashtable_t*);

/**
 * @brief Test hashtable implementation.
 * @return True if it is correct; otherwise, false.
 */
bool hashtableTest();

/**
 * @brief Free hashtable.
 * @param Hashtable to be freed.
 */
void hashtableFree(hashtable_t*);

// ----- End of HASHTABLE ----- //

// ----- RBT ----- //

/**
 * @brief Create new RBT node.
 * @param Value that the RBT node should contain.
 * @return Created RBT node.
 */
rbtNode_t* createRbtNode(const int);

/**
 * @brief Create new RBT.
 * @return Created RBT.
 */
rbt_t* createRbt();

/**
 * @brief Left rotate operation.
 * @param The RBT.
 * @param The RBT node to rotate on.
 */
void rbtLeftRotate(rbt_t*, rbtNode_t*);
/**
 * @brief Right rotate operation.
 * @param The RBT.
 * @param The RBT node to rotate on.
 */
void rbtRightRotate(rbt_t*, rbtNode_t*);

/**
 * @brief Insert RBT node in th RBT.
 * @param The RBT.
 * @param The RBT node to be inserted.
 */
void rbtInsert(rbt_t*, rbtNode_t*);

/**
 * @brief Fixup function for RBT insertion.
 * @param The RBT the be fixed.
 * @param The initial RBT node to be fixed.
 */
void rbtInsertFixup(rbt_t*, rbtNode_t*);

/**
 * @brief Search for a value in the RBT.
 * @param The RBT.
 * @param Value to be searched.
 * @return RBT node containing the value, if it exists; otherwise, NULL.
 */
rbtNode_t* rbtSearch(rbt_t*, const int);

/**
 * @brief Print RBT in order.
 * @param RBT to be printed.
 * @param RBT node to be printed.
 */
void rbtInOrder(rbt_t*, rbtNode_t*);

/**
 * @brief Test RBT implementation.
 * @return True if it is correct; otherwise, false.
 */
bool rbtTest();

/**
 * @brief Check if the tree is actually a RBT.
 * @param Tree to be checked.
 * @return True if it is; otherwise, false.
 */
bool isRbt(rbt_t*);

/**
 * @brief Function that checks if the tree has the BST property (i.e., x->left->value < x->value <= x->right->value, for all x).
 * @param Tree to be checked.
 * @return True if it is; otherwise, false.
 */
bool rbtHasBstProperty(rbt_t*);

/**
 * @brief Utility function for checking if the tree has the BST property.
 * @param Tree to be checked.
 * @param Current RBT node.
 * @param RBT test data structure.
 */
void rbtHasBstPropertyUtil(rbt_t*, rbtNode_t*, rbtTestStructure_t*);

/**
 * @brief Function that computes the black height of the RBT.
 * @param The RBT.
 * @param Current RBT node.
 * @return Black height if all paths have the same black height; otherwise, -1.
 */
int rbtComputeBlackHeight(rbt_t*, rbtNode_t*);

/**
 * @brief Free RBT nodes.
 * @param RBT whose nodes must be freed.
 * @param RBT node to be freed.
 */
void rbtFreeNodes(rbt_t*, rbtNode_t*);

/**
 * @brief Free RBT.
 * @param RBT to be freed.
 */
void rbtFree(rbt_t*);

// ----- End of RBT ----- //

// ----- AUXILIARY FUNCTIONS ----- //
/**
 * @brief Generate a collection of random numbers.
 * @param Array of random numbers.
 * @param Size of the array.
 */
void generateRandomArray(int*, const int);

/**
 * @brief Unit test: check if the input array is sorted.
 * @param Array to be checked if sorted.
 * @param Size of the array.
 * @return True if it is sorted; otherwise, false
 */
bool isSorted(const int*, const int);

// ----- End of AUXILIARY FUNCTIONS ----- //

// ----- CORE FUNCTIONS ----- //

/**
 * @brief Function that does the experiment.
 * @param Array of random numbers.
 * @param Number of insertion operations.
 * @param Number of search operations.
 * @param Data structure to be used. The possible values are:
 * @return Elapsed time for the experiment.
 */
clock_t doExperiment(int*, const unsigned int, const unsigned int, char*);

// ----- End of CORE FUNCTIONS ----- //

// ##### End of PROTOTYPES OF THE FUNCTIONS ##### //

int main() {

    // Random seed initialization.
    srand(RANDOM_SEED);
    // Elapsed time for hashtable.
    clock_t timeHashtable = 0;
    // Elapsed time for RBT.
    clock_t timeRbt = 0;
    // Number of insert operations.
    unsigned int numInsertions = 0;
    // Number of search operations.
    unsigned int numSearches = 0;

    // What is the outputPointer?
    if (outputType == ONCONSOLE || outputType == ONFILE) {
        // On console.
        if (outputType== ONCONSOLE) outputPointer = stdout;
        // On file.
        else {
            // Open file.
            outputPointer = fopen("results.txt", "w");
            // Have we opened the file?
            if (outputPointer == NULL) {
                fprintf(stderr, "ERROR: The outputPointer has not been created\n");
                exit(-1);
            }
        }
    }
    // Error
    else {
        fprintf(stderr, "ERROR: The outputType can be only ONCONSOLE or ONFILE\n");
        exit(-1);
    }

    // Print the header, only if itONFILE is on console.
    if (outputType == ONCONSOLE) {
        fprintf(outputPointer, "+-----------------------------+---------------------+---------------------+\n");
        fprintf(outputPointer, "| Operations - %%I & %%S        | Hashtable - %-5d   | Red Black Tree      |\n", NUM_ENTRIES);
        fprintf(outputPointer, "+-----------------------------+---------------------+---------------------+\n");
    }

    // For each number of operations in the interval [MIN_OPERATIONS, MAX_OPERATIONS] with STEP
    for (int numOps=MIN_OPERATIONS; numOps<=MAX_OPERATIONS; numOps+=STEP) {
        // Reset the times.
        timeHashtable = timeRbt = 0;
        // For each experiment
        for (int exper=1; exper<=NUM_EXPERIMENTS; exper++) {

            // Compute the number of insert operations.
            numInsertions = numOps*PERCENTAGE_INSERTIONS/100;
            // Compute the number of search operations.
            numSearches = numOps-numInsertions;
            // Allocate numInsertions memory cells for the array of random numbers.
            int* randomArray = malloc(numInsertions*sizeof(int));
            // Fill-in the array with random numbers.
            generateRandomArray(randomArray, numInsertions);
            // Hashtable experiment.
            timeHashtable += doExperiment(randomArray, numInsertions, numSearches, "hashtable");
            // RBT experiment.
            timeRbt += doExperiment(randomArray, numInsertions, numSearches, "rbt");
            // Free the array of random numbers.
            free(randomArray);
        }
        // Printing the (sample mean as) result. Use TAB (\t) on file.
        if (outputType == ONCONSOLE)
            fprintf(outputPointer, "| %15d - %-3d & %-3d | %19f | %19f |\n",
                numOps,
                PERCENTAGE_INSERTIONS,
                100-PERCENTAGE_INSERTIONS,
                (float) timeHashtable/NUM_EXPERIMENTS,
                (float) timeRbt/NUM_EXPERIMENTS);
        else
            fprintf(outputPointer, "%d \t%f \t%f \n",
                numOps,
                (float) timeHashtable/NUM_EXPERIMENTS,
                (float) timeRbt/NUM_EXPERIMENTS);
    }

    // Print the ending part, only if it is on console.
    if (outputType == ONCONSOLE) {
        fprintf(outputPointer, "+-----------------------------+---------------------+---------------------+\n");
        fprintf(outputPointer, "| Legend:                                                                 |\n");
        fprintf(outputPointer, "|                                                                         |\n");
        fprintf(outputPointer, "| %%I: Percentage of insertion operations                                  |\n");
        fprintf(outputPointer, "| %%S: Percentage of search operations                                     |\n");
        fprintf(outputPointer, "|                                                                         |\n");
        fprintf(outputPointer, "| The number near \"Hashtable\" is the number of entries in the hashtable   |\n");
        fprintf(outputPointer, "+-------------------------------------------------------------------------+\n");
    }

    if (TEST_DATA_STRUCTURES) {
        fprintf(outputPointer, "| Hashtable implementation: %-12s                                  |\n", hashtableTest() ? "correct" : "not correct");

        fprintf(outputPointer, "| Red black tree implementation: %-12s                             |\n", rbtTest() ? "correct" : "not correct");
        fprintf(outputPointer, "+-------------------------------------------------------------------------+\n");
    }

    // Return 0.
    return 0;
}

// ##### IMPLEMENTATION OF THE FUNCTIONS ##### //

// ----- LINKED LIST ----- //

/**
 * @brief Create a new linked list node.
 * @param v Value that the linked list node should contain.
 * @return Created linked list node.
 */
linkedListNode_t* createLinkedListNode(const int v) {
    linkedListNode_t *x = (linkedListNode_t *) malloc(sizeof(linkedListNode_t));
    x->value=v;
    x->next=NULL;
    x->prev=NULL;
    return x;
}

/**
 * @brief Create a new linked list.
 * @return Created linked list.
 */
linkedList_t* createLinkedList() {
    linkedList_t *l = (linkedList_t *) malloc(sizeof(linkedList_t));
    l->head=NULL;
    l->size=0;
    return l;
}

/**
 * @brief Insert linked list node in the head of the linked list.
 * @param list The linked list.
 * @param x Linked list node to be inserted.
 */
void linkedListInsert(linkedList_t* l, linkedListNode_t* x) {
    x->next=l->head;
    if(l->head != NULL)
        l->head->prev=x;
    l->head = x;
    x->prev=NULL;
    l->size++;
}

/**
 * @brief Search for a value in the linked list.
 * @param list The linked list.
 * @param v Value to be searched for.
 * @return First linked list node containing such value, if it exists; otherwise, NULL.
*/
linkedListNode_t* linkedListSearch(linkedList_t* list, const int v) {
    linkedListNode_t *x;
    x = list->head;
    while(x!=NULL && x->value != v)
        x = x->next;
    return x;
}

/**
 * @brief Delete a linked list node from linked list.
 * @param list The linked list.
 * @param x The linked list node to be deleted.
 */
void linkedListDelete(linkedList_t* list, linkedListNode_t* x) {
    if(x->prev != NULL)
    {
        x->prev->next=x->next;
    }else{
        list->head=x->next;
    }
    if(x->next != NULL)
        x->next->prev=x->prev;
    list->size = list->size -1;
    free(x);
}

/**
 * @brief Print the linked list.
 * @param list Linked list to be printed.
 */
void linkedListPrint(linkedList_t* list) {
    linkedListNode_t* x = list->head;
    while (x) {
        fprintf(stdout, "%d ", x->value);
        x = x->next;
    }
}

/**
 * @brief Free the linked list.
 * @param list Linked list to be freed.
 */
void linkedListFree(linkedList_t* list) {
    
    linkedListNode_t *del;
    while(list->head!=NULL){
        del = list->head;
        list->head = list->head->next;
        free(del);
    }
        free(list);

}

// ----- End of LINKED LIST ----- //

// ----- HASHTABLE ----- //

/**
 * @brief Create a new hashtable.
 * @param s The size of the hashtable (i.e., the number of entries).
 * @return The created hashtable.
 */
hashtable_t* createHashtable(const unsigned int s) {
  hashtable_t* hash = (hashtable_t*)malloc(sizeof(hashtable_t));
  hashtableEntry_t *head;
  if(!hash){
    return NULL;
  }
  hash->entry = malloc(sizeof(hashtableEntry_t)*s);
  if(!hash->entry){
    return NULL;
  }
  hash->size=s;
  for(int i = 0; i<hash->size;i++){
     head = malloc(sizeof(hashtableEntry_t));
     if(!head){
        return NULL;
  }
    head->list = createLinkedList();
    hash->entry[i] = head;
  }
    return hash;

}

/**
 * @brief Hash function computing the key for a given integer.
 * @param hashtbl The hashtable needed to access the size of it.
 * @param v The integer for which the key must be computed.
 * @return The computed key.
 */
const unsigned int hashFunction(hashtable_t* hashtbl, const int v) {
    return v % hashtbl->size;
}

/**
 * @brief Insert value in the hashtable.
 * @param hashtbl The hashtable.
 * @param v Value to be inserted.
 */
void hashtableInsert(hashtable_t* hashtbl, const int v) {
    int hash = hashFunction(hashtbl,v);
    linkedListNode_t *x = createLinkedListNode(v);
    linkedListInsert(hashtbl->entry[hash]->list,x);
}

/**
 * @brief Search for a value in the hashtable.
 * @param hashtbl The hashtable.
 * @param v Value to be searched.
 * @return Linked list node containing such value, if it exists; otherwise, NULL.
 */
linkedListNode_t* hashtableSearch(hashtable_t* hashtbl, const int v) {
    int i = hashFunction(hashtbl,v);
    return linkedListSearch(hashtbl->entry[i]->list,v);
}

/**
 * @brief Delete value from hashtable.
 * @param hashtbl The hashtable.
 * @param x Linked list node to be deleted.
 */
void hashtableDelete(hashtable_t* hashtbl, linkedListNode_t* x) {
    int i = hashFunction(hashtbl,x->value);
    linkedListDelete(hashtbl->entry[i]->list,x);
}

/**
 * @brief Print the hashtable.
 * @param hashtbl Hashtable to be printed.
 */
void hashtablePrint(hashtable_t* hashtbl) {
    for (int i=0; i<hashtbl->size; i++) {
        fprintf(stdout, "%d => ", i);
        linkedListPrint(hashtbl->entry[i]->list);
        fprintf(stdout, "\n");
    }
}

/**
 * @brief Test hashtable if it is correctly implemented.
 * @return True if it is correct; otherwise, false.
 */
bool hashtableTest() {
    bool test;
    int A[NUM_ELEMENTS_FOR_TEST];
    hashtable_t * hash = createHashtable(NUM_ELEMENTS_FOR_TEST);
    linkedListNode_t *nodo;
    for (int i = 0; i<NUM_ELEMENTS_FOR_TEST;i++){
        hashtableInsert(hash,A[i]);
        nodo = hashtableSearch(hash,A[i]);
        if(nodo->value==A[i])
            test=true;
        else
            test=false;
    }
    hashtableFree(hash);
    return test;
}

/**
 * @brief Free hashtable.
 * @param hashtbl Hashtable to be freed.
 */
void hashtableFree(hashtable_t* hashtbl) {
    for(int i=0;i<hashtbl->size;i++){
        linkedListFree(hashtbl->entry[i]->list);
        free(hashtbl->entry[i]);
    }
    free(hashtbl->entry);
    free(hashtbl);
}

// ----- End of HASHTABLE ----- //

// ----- RBT ----- //

/**
 * @brief Create new RBT node.
 * @param v Value that the RBT node should contain.
 * @return Created RBT node.
 */
rbtNode_t* createRbtNode(const int v) {
    rbtNode_t *nodo = (rbtNode_t*) malloc(sizeof(rbtNode_t));
    nodo->left = NULL;
    nodo->right = NULL;
    nodo->parent = NULL;
    nodo->value = v;
    nodo->color = 'R';
    return nodo;
}

/**
 * @brief Create new RBT.
 * @return Created RBT.
 */
rbt_t* createRbt() {
    rbt_t *rbt = (rbt_t*)malloc(sizeof(rbt_t));
    rbt->size = 0;
    rbt->nil = createRbtNode(0);
    rbt->nil->color = 'B';
    rbt->root=rbt->nil;
    return rbt;

}

/**
 * @brief Left rotate operation.
 * @param rbt The RBT.
 * @param x The RBT node to rotate on.
 */
void rbtLeftRotate(rbt_t* rbt, rbtNode_t* x) {
    rbtNode_t *y;
    y = x->right;
    x->right=y->left;
    if (!(y->left ==rbt->nil)){
        y->left->parent = x;
    }
    y->parent = x->parent;
    if(x->parent == rbt->nil){
        rbt->root= y;
    }
    if(!(x->parent == rbt->nil) && x == x->parent->left){
            x->parent->left= y;
    }
    if(!(x->parent == rbt->nil) && x == x->parent->right){
        x->parent->right = y;
    }
    y->left = x;
    x->parent  = y;

}

/**
 * @brief Right rotate operation.
 * @param rbt The RBT.
 * @param x The RBT node to rotate on.
 */
void rbtRightRotate(rbt_t* rbt, rbtNode_t* x) {
    rbtNode_t *y;
    y = x->left;
    x->left=y->right;
    if(y->right != rbt->nil){
        y->right->parent = x;
    }
    y->parent = x->parent;
    if(x->parent == rbt->nil){
        rbt->root= y;
    }
    if(!(x->parent == rbt->nil) && x == x->parent->right){
            x->parent->right= y;
    }
    if(!(x->parent == rbt->nil) && x == x->parent->left){
        x->parent->left = y;
    }
    y->right = x;
    x->parent  = y;

}

/**
 * @brief Insert RBT node in th RBT.
 * @param rbt The RBT.
 * @param z The RBT node to be inserted.
 */
void rbtInsert(rbt_t* rbt, rbtNode_t* z) {
    rbt->size++;
    rbtNode_t *y,*x;
    y = rbt->nil;
    x = rbt->root;
    while ( !(x == rbt->nil)){
        y = x;
        if(z->value < x->value){
            x = x->left;
        }else{
            x = x->right;
        }
    }
    z->parent = y;
    if( y == rbt->nil){
        rbt->root =z;
    }
    if(!(y==rbt->nil) && z->value < y->value){
        y->left = z;
    }
    if(!(y==rbt->nil) && z->value >= y->value){
        y->right = z;
    }
    z->left= rbt->nil;
    z->right=rbt->nil;
    z->color = 'R';
    rbtInsertFixup(rbt,z);

}

/**
 * @brief Fixup function for RBT insertion.
 * @param rbt The RBT the be fixed.
 * @param z The initial RBT node to be fixed.
 */
void rbtInsertFixup(rbt_t* rbt, rbtNode_t* z) {
    rbtNode_t *y;
    while(z->parent->color == 'R'){
        if(z->parent == z->parent->parent->left){
            y= z->parent->parent->right;
            if(y->color=='R'){
                z->parent->color = 'B';
                y->color = 'B';
                z->parent->parent->color = 'R';
                z = z->parent->parent;
            }else{
                if(z == z->parent->right){
                    z = z->parent;
                    rbtLeftRotate(rbt,z);
                }
                z->parent->color='B';
                z->parent->parent->color= 'R';
                rbtRightRotate(rbt,z->parent->parent);
            }
        }
        else
            {
            y = z->parent->parent->left;
            if(y->color == 'R'){
                z->parent->color = 'B';
                y->color = 'B';
                z->parent->parent->color = 'R';
                z = z->parent->parent;
            }
            else
            {
                if(z == z->parent->left){
                    z = z->parent;
                    rbtRightRotate(rbt,z);
                }
                z->parent->color='B';
                z->parent->parent->color= 'R';
                rbtLeftRotate(rbt,z->parent->parent);
            }
        }
    }
    rbt->root->color= 'B';
}


/**
 * @brief Search for a value in the RBT.
 * @param rbt The RBT.
 * @param v Value to be searched.
 * @return RBT node containing the value, if it exists; otherwise, NULL.
 */
rbtNode_t* rbtSearch(rbt_t* rbt, const int v) {
    rbtNode_t *x = rbt->root;
    while((x!=rbt->nil)&&(v!=x->value)){
        if(v<x->value)
            x=x->left;
        else
            x=x->right;     
    }
    return x;

}

/**
 * @brief Print RBT in order.
 * @param rbt RBT to be printed.
 * @param x RBT node to be printed.
 */
void rbtInOrder(rbt_t* rbt, rbtNode_t* x) {
    if (x != rbt->nil) {
        rbtInOrder(rbt, x->left);
        printf("%d ", x->value);
        rbtInOrder(rbt, x->right);
    }
}

/**
 * @brief Test RBT if it is correctly implemented.
 * @return True if it is correct; otherwise, false.
 */
bool rbtTest() {
    bool risultato,risultato1,risultato2;
    rbt_t * rbt = createRbt();
    rbtNode_t *nodo;
    int key = 0;
    for (int i = 0; i< NUM_ELEMENTS_FOR_TEST ;i++){
        key = i;
        nodo = createRbtNode(key);
        rbtInsert(rbt,nodo);
        nodo = rbtSearch(rbt,key);
        if(nodo->value == key){
            risultato = true;
        }else{
            risultato = false;
        }
    }
    risultato1 = isRbt(rbt);
    risultato2 = rbtHasBstProperty(rbt);
    rbtFree(rbt);
    return risultato && risultato1 && risultato2;
}

/**
 * @brief Check if the tree is actually a RBT.
 * @param rbt Tree to be checked.
 * @return True if it is; otherwise, false.
 */

bool isRbtp1(rbt_t * rbt,rbtNode_t *x){
    bool result;
    if (x != rbt->nil ){
        if(x->color == 'R' || x->color == 'B'){
        isRbtp1(rbt, x->left);
        isRbtp1(rbt, x->right);
        CONT++;
        result = true;
        }else{
            result = false;
        }
    }
    return result;
}
bool isRbtp4(rbt_t * rbt,rbtNode_t *x){
    if (x != rbt->nil){
    if(x->color == 'R')
    {
        if(x->left->color == 'B' && x->right->color == 'B')
        {
            R = true;
        }
        else
        {
            R = false;
        }
    }
    isRbtp4(rbt, x->left);
    isRbtp4(rbt, x->right);
    }
    return R;
}

bool isRbt(rbt_t* rbt) {
    bool risultato1, risultato2, risultato3, risultato4;
    rbtNode_t *node = rbt->root;
    //propriet? 1: ogni nodo ? rosso o nero
    risultato1 = isRbtp1(rbt,node);

    //propriet? 2: la radice ? nera
    if(rbt->root->color == 'B'){
        risultato2 = true;
    }else{
        risultato2 = false;
    }


    //propriet? 3: ogni foglia esterna = T.nil ? nera

    if(rbt->nil->color == 'B'){
        risultato3 = true;
    }else{
        risultato3 = false;
    }


    node = rbt->root;
    //propriet? 4: se un nodo ? rosso entrambi i figli sono neri
    risultato4 = isRbtp4(rbt,node);

    return risultato1 && risultato2 && risultato3 && risultato4;

}

/**
 * @brief Function that checks if the tree has the BST property (i.e., x->left->value < x->value <= x->right->value, for all x).
 * @param rbt Tree to be checked.
 * @return True if it is; otherwise, false.
 */
bool rbtHasBstProperty(rbt_t* rbt) {
    rbtNode_t *x=rbt->root;
    bool check1,check2;
    rbtTestStructure_t *testStructure = (rbtTestStructure_t*)malloc(sizeof(rbtTestStructure_t));
    testStructure->index = 0;
    testStructure->A = malloc(sizeof(rbtTestStructure_t)*rbt->size);
    rbtHasBstPropertyUtil(rbt,x,testStructure);
    check2 = isSorted(testStructure->A,testStructure->index);
    free(testStructure->A);
    free(testStructure);
    x = rbt->root;
    if(x!= rbt->nil && x->left != rbt->nil && x->right != rbt->nil){
        if(x->left->value < x->value && x->right->value >= x->value){
                check1 = true;
        }else{
            check1 = false;
        }
    }
    return check1 && check2;
}

/**
 * @brief Utility function for checking if the tree has the BST property.
 * @param rbt Tree to be checked.
 * @param x Current RBT node.
 * @param rbtTestStructure RBT test data structure.
 */
void rbtHasBstPropertyUtil(rbt_t* rbt, rbtNode_t* x, rbtTestStructure_t* rbtTestStructure) {
    if (x != rbt->nil) {
        rbtHasBstPropertyUtil(rbt, x->left,rbtTestStructure);
        rbtTestStructure->A[rbtTestStructure->index] = x->value;
        rbtTestStructure->index++;
        rbtHasBstPropertyUtil(rbt, x->right,rbtTestStructure);
    }

}

/**
 * @brief Function that computes the black height of the RBT.
 * @param rbt The RBT.
 * @param x Current RBT node.
 * @return Black height if all paths have the same black height; otherwise, -1.
 */
int rbtComputeBlackHeight(rbt_t* rbt, rbtNode_t* x) {
    int left,right;
    if(x==rbt->nil)
        return 1;
    else{
        right=rbtComputeBlackHeight(rbt,x->right);
        left=rbtComputeBlackHeight(rbt,x->left);
        if(left == -1 || right == -1 || left!=right)
            return -1;
        else
        {
            if(x->color == 'B')
               return left + 1;
            else
                return left;
       }
    }
}

/**
 * @brief Free RBT nodes.
 * @param T RBT whose nodes must be freed.
 * @param x RBT node to be freed.
 */
void rbtFreeNodes(rbt_t* T, rbtNode_t* x) {
    if(x!=T->nil){
        rbtFreeNodes(T,x->left);
        rbtFreeNodes(T,x->right);
        free(x);
    }
}

/**
 * @brief Free RBT.
 * @param T RBT to be freed.
 */
void rbtFree(rbt_t* T) {
    rbtNode_t* x = T->root;
    rbtFreeNodes(T,x);
    free(T->nil);
    free(T);
}

// ----- End of RBT ----- //

// ----- AUXILIARY FUNCTIONS ----- //

/**
 * @brief Generate a collection of random numbers.
 * @param A Array of random numbers.
 * @param n Size of the array.
 */
void generateRandomArray(int* A, const int n) {
    // For each i in 0..n-1, generate a random number.
    for (int i = 0; i < n; i++) A[i] = rand() % MAX_RANDOM_NUMBER;
}

/**
 * @brief Unit test: check if the input array is sorted.
 * @param A Array to be checked if sorted.
 * @param n Size of the array.
 * @return True if it is sorted; otherwise, false
 */
bool isSorted(const int* A, const int n) {
    // For each i in 0..n-2, if the current element is greater than the next one,
    // then it is unsorted.
    for (int i = 0; i < n-1; i++) if (A[i] > A[i+1]) return false;
    // Otherwise it is.
    return true;
}

// ----- End of AUXILIARY FUNCTIONS ----- //

// ----- CORE FUNCTIONS ----- //

/**
 * @brief Function that does the experiment.
 * @param randomArray Array of random numbers.
 * @param numInsertions Number of insertion operations.
 * @param numSearches Number of search operations.
 * @param dataStructure Data structure to be used. The possible values are: hashtable and rbt.
 * @return Elapsed time for the experiment.
 */
clock_t doExperiment(int* randomArray, const unsigned int numInsertions, const unsigned int numSearches, char* dataStructure) {
    int key;
    clock_t start,end=0;
    hashtable_t* hashTable = createHashtable(NUM_ENTRIES);
    rbt_t* rbt= createRbt();
    rbtNode_t* nodeRbt;
    linkedListNode_t* nodeHashTable;
    start=clock();
    if(strcmp(dataStructure, "hashtable")==0){
        for(int i=0; i<numInsertions; i++){
            key=rand() % MAX_RANDOM_NUMBER+1;
            hashtableInsert(hashTable,key);
        }
        for(int j=0; j<numSearches; j++){
            key=rand()% MAX_RANDOM_NUMBER +1;
            nodeHashTable=hashtableSearch(hashTable,key);
        }
    }
    else if (strcmp(dataStructure, "rbt")==0)
    {
        for(int i=0; i<numInsertions;i++){
            key=rand()%MAX_RANDOM_NUMBER +1;
            nodeRbt=createRbtNode(key);
            rbtInsert(rbt,nodeRbt);
        }
        for(int j=0;j<numSearches; j++){
            key=rand()%MAX_RANDOM_NUMBER +1;
            nodeRbt=rbtSearch(rbt,key);
        }
    }
    else{
        fprintf(stderr, "ERROR: There is no such sorting alghoritm called %s \n", dataStructure);
        exit(1);
    }
    end=clock();
    hashtableFree(hashTable);
    rbtFree(rbt);
    return end-start;
}

// ----- End of CORE FUNCTIONS ----- //

// ##### End of IMPLEMENTATION OF THE FUNCTIONS ##### //


