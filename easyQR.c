#include "easyQR.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////* private types *//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Matrix type in Compressesed Sparse Row (CSR) format
 * 
 * see https://en.wikipedia.org/wiki/Sparse_matrix for details of the data format
 * 
 */
typedef struct CSR{
	int rows;						/**< number of rows */ 
	int* row_start;					/**< indices that point to the beginning of each row in indices */
	int* indices;					/**< list of all indices */
	double* elements;				/**< list of all values */
}CSR_matrix;



/**
 * Matrix type in Compressesed Sparse Column (CSC) format
 * 
 * see https://en.wikipedia.org/wiki/Sparse_matrix for details of the data format
 * 
 */
typedef struct CSC{
	int cols;						/**< number of columns */ 
	int* col_start;					/**< indices that point to the beginning of each row in indices */
	int* indices;					/**< list of all indices */
	double* elements;				/**< list of all values */
}CSC_matrix;



/**
 * Matrix-element used for storage of the "Matrix-Market"-format
 * The element A_{ij} is represented by a row(i) and column(j) index and the corresponding value (A_{ij})
 * 
 * see http://math.nist.gov/MatrixMarket/formats.html for more info
*/
typedef struct MATRIX_ELEMENT{
	int row_index;					/**< row index */
	int col_index;					/**< column index */
	double value;					/**< real value */
}matrix_element;



/**
 * Element of a real eigensystem consisting of an eigenvalue, its corrsponding eigenvector, and its dimension
 */
typedef struct EIGEN{
	int dim;						/**< dimension of the eigenvector*/
	double value;					/**< (real) eigenvalue*/
	double* vector;					/**< (real) eigenvector as a list of <dim> doubles*/
}eigen;



/**
 * Color type for the red-black nodes
 */
typedef enum RBCOLOR{red,black}RBcolor;



/**
 * Relational-indicator type for red-black trees
 */
typedef enum RELATION{LEFTLEFT,LEFTRIGHT,RIGHTLEFT,RIGHTRIGHT}relation;



/**
 * Red-Black-Node type
 * 
 *see  https://en.wikipedia.org/wiki/Red%E2%80%93black_tree
 */
typedef struct RBNODE{
	struct RBNODE* parent; 	/**< pointer to parent node, if node is root then it contains NULL */
	struct RBNODE* left;	/**< pointer to left child */
	struct RBNODE* right;	/**< pointer to right child */
	RBcolor color;			/**< color of the node (either red or black) */
	void* data;				/**< pointer to the data represented by the node */
}rbNode;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////* private global variables *//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static int verbose_level = 0;						/**< Verbose level of plotting infromation during computation in range 0-3, 0 is quiet*/
static int (*compare)(void* X1,void* X2) = NULL;	/**< Global variable where a pointer to the compare function for red-black trees is stored*/
static void (*free_data)(void* X) = NULL;			/**< Global variable where a pointer to the free function for red-black trees is stored*/
static size_t data_size_in_bytes = 0;				/**< Global variable that holds the size of a data element of a red-black node in bytes*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////* private functions *////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////// red-black-tree functions /////////////////////////////////////////////////////////

/**
 * \brief Sets user-defined compare functione for red-black ordering
 * \param int (*Compare)(void* X1,void* X2): compare function
 */
static void RBTree_set_compare(int (*Compare)(void* X1,void* X2)){
	compare = Compare;	
}


/**
 * \brief Sets user-defined free function for the red-black tree.
 * \param void (*Free_data)(void* X): free function
 */
static void RBTree_set_free(void (*Free_data)(void* X)){
	free_data = Free_data;	
}


/**
 * \brief Sets the size of a user-defined data element in bytes
 * \param size_t data_size: size of a data element in bytes 
 */
static void RBTree_set_data_size(size_t data_size){
	data_size_in_bytes = data_size;
}


/**
 * \brief Creates on node of a red-black tree
 * \param rbNode* parent: parent node of the new node
 * \param void* data: pointer to the data that the node contains (note: that data is not copied, just the pointer is stored)
 * \return rbNode* : pointer to the newly created node
 */
static rbNode* RBTcreateNode(rbNode* parent,void* data){
	rbNode* node = (rbNode*)malloc(sizeof(rbNode));
	node->parent = parent;
	node->left = NULL;
	node->right = NULL;
	if (parent==NULL) node->color = black; else node->color = red;
	node->data = data;
	return node;
}


/**
 * \brief Frees a whole red-black tree
 * \param rbNode root: pointer to the root element of the red-black tree
 */
void RBTfree(rbNode* root){
	if (root!=NULL){
		if (free_data!=NULL && root->data!=NULL) (*free_data)(root->data);
		RBTfree(root->left);
		RBTfree(root->right);
		free(root);
	}
}


/**
 * \brief Gets the uncle of a given red-black node
 * \param rbNode* node: pointer to red-black node. 
 * \return rbNode* : the uncle node, if <node> is root then NULL
 */
static rbNode* getUncle(rbNode* node){
	rbNode* P = node->parent;
	if (P!=NULL){
		rbNode* G = P->parent;	
		if (G->left==P) return G->right;else return G->left;
	}
	else return NULL;
}


/** 
 * \brief Helper function for "RBfixTRee"
 * \param rbNode* child: child node 
 * \param rbNode* parent: parent node of <child>
 * \param rbNode* grand: parent node of <parent>
 * \return relation: relational case
 */
static relation getRelationCase(rbNode* child,rbNode* parent,rbNode* grand){
	if (parent==grand->left){
		if (child==parent->left) return LEFTLEFT; else return LEFTRIGHT;
	}
	else{
		if (child==parent->right) return RIGHTRIGHT; else return RIGHTLEFT;
	}
}


/**
 * \brief Rotates child and parent node in a red-black tree
 * \param rbNode* child: child node 
 * \param rbNode* parent: parent node of <child>
 */
static void RBTrotation(rbNode* child,rbNode* parent){
	rbNode* grand = parent->parent;
	if (grand!=NULL){			
		if (grand->left==parent) grand->left = child; else grand->right = child;
	}
	
	if (child==parent->left){				
		if (child->right!=NULL) child->right->parent = parent;			
		parent->left = child->right;
		child->right = parent;		
	}
	else{
		if (child->left!=NULL) child->left->parent = parent;
		parent->right = child->left;
		child->left = parent;
	}
	
	child->parent = grand;
	parent->parent = child;	
}


/**
 * \brief Helper function for "RBTinsert" that removes inconsistencies of red-black rules 
 * \param rbNode* node: node of a subtree that may be inconsistent 
 */
static void RBfixtree(rbNode* node){
	rbNode* parent = node->parent;
	if (parent==NULL){
		node->color = black;
		return;
	}
	
	if (node->color==red && parent->color==red){
		rbNode* grand = parent->parent;		
		rbNode* uncle = getUncle(node);
		
		if (uncle!=NULL && uncle->color==red){
			parent->color = black;
			uncle->color = black;
			if (grand->parent!=NULL) grand->color = red;
			RBfixtree(grand);			
		}
		else{
			relation rel = getRelationCase(node,parent,grand);
			switch(rel){
				case LEFTLEFT:
					parent->color = black;
					grand->color = red;
					RBTrotation(parent,grand);
					break;
				case LEFTRIGHT:
					node->color = black;
					grand->color = red;
					RBTrotation(node,parent);					
					RBTrotation(node,grand);				
					break;
				case RIGHTLEFT:
					node->color = black;
					grand->color = red;
					RBTrotation(node,parent);					
					RBTrotation(node,grand);
					break;
				case RIGHTRIGHT:
					parent->color = black;
					grand->color = red;
					RBTrotation(parent,grand);
					break;
			}					
		}		
	}
}


/**
 * \brief Finds the root of some given node in a red-black tree
 * \param rNode** root: pointer to address of some red-black node. The result is stored by the same pointer.
 */
static void RBTfindroot(rbNode** root){
	if (*root!=NULL){
		while((*root)->parent!=NULL){
			*root = (*root)->parent;
		}
	}
}


/**
 * \brief Helper function for "RBTinsertElement"
 * \param rbNode** root: pointer to address of the root of the red-black tree
 * \param rbNode* injector: node where the newly created node is inserted
 * \param void* data: pointer to a data structure that the new node contains
 * \return rbNode*: If a node with the same data is found, no new node is created and the pointer to this node is returned, else NULL is returned
 */
static rbNode* RBTinsert(rbNode** root,rbNode* injector,void* data){
	if (*root==NULL){
		*root = RBTcreateNode(NULL,data);
		return NULL;
	}	
	
	int cmp = (*compare)(data,injector->data);
	if (cmp!=0){
		rbNode* child = (cmp<0) ? injector->left : injector->right;
		if (child==NULL){	
			child = RBTcreateNode(injector,data);		
			if (cmp<0) injector->left = child; else injector->right = child;
			RBfixtree(child);
			RBTfindroot(root);
			return NULL;
		}
		else return RBTinsert(root,child,data);
	}
	else return injector;	
}


/**
 * \brief If no node with the same data exists in a red-black tree, a new node is created an inserted in the tree.
 * \param rbNode** root: pointer to address of the root node of the tree
 * \param void* data: pointer to a data structure that the new node contains
 * \return rbNode*: If a node with the same data is found, no new node is created and the pointer to this node is returned, else NULL is returned
 */
static rbNode* RBTinsertElement(rbNode** root,void* data){
	return RBTinsert(root,*root,data);
}

/**
 * \brief Gets the node with minimal data with respect to the defiend compare function
 * \param rbNode* root: pointer to root node of the red-black tree
 * \param rbNode*: pointer to minimal node
 */
rbNode* RBTminNode(rbNode* root){
	if (root==NULL) return NULL;
	while(root->left!=NULL){
		root = root->left;
	}
	return root;
}


/**
 * \brief Gets the node with maximal data with respect to the defiend compare function
 * \param rbNode* root: pointer to root node of the red-black tree
 * \param rbNode*: pointer to maximal node
 */
rbNode* RBTmaxNode(rbNode* root){
	if (root==NULL) return NULL;
	while(root->right!=NULL){
		root = root->right;	
	}
	return root;
}


/**
 * \brief Helper function for "RBTpredecessor"
 * \param rbNode* root: pointer to root node of the red-black tree
 * \return rbNode* : pointer to first upper node that has a right child, else NULL
 */
static rbNode* RBTupUntilRightChild(rbNode* node){
	if (node==NULL) return NULL;
	while(node->parent!=NULL && node!=node->parent->right){
		node = node->parent;
	}
	return node->parent;
}


/**
 * \brief Helper function for "RBTsuccessor"
 * \param rbNode* root: pointer to root node of the red-black tree
 * \return rbNode* : pointer to first upper node that has a left child, else NULL
 */
static rbNode* RBTupUntilLeftChild(rbNode* node){	
	if (node==NULL) return NULL;
	while(node->parent!=NULL && node!=node->parent->left){
		node = node->parent;
	}
	return node->parent;
}


/**
 * \brief Gets the previous node in order of the compare function
 * \param rbNode* root: pointer to root node of the red-black tree
 * \param rbNode* : ponter to the predecessor node, if none exists then NULL
 */
rbNode* RBTpredecessor(rbNode* node){
	if (node!=NULL){
		if (node->left!=NULL) return RBTmaxNode(node->left);
		else return RBTupUntilRightChild(node);				
	}
	else return NULL;
}


/**
 * \brief Gets the next node in order of the compare function
 * \param rbNode* root: pointer to root node of the red-black tree
 * \param rbNode* : ponter to the successor node, if none exists then NULL
 */
rbNode* RBTsuccessor(rbNode* node){
	if (node!=NULL){
		if (node->right!=NULL) return RBTminNode(node->right);
		else return RBTupUntilLeftChild(node);			
	}
	else return NULL;
}


/**
 * \brief Gets the number of nodes in a red-black tree
 * \param rbNode* root: pointer to root node of the red-black tree
 */
int RBTnodeCount(rbNode* root){
	if (root!=NULL){
		return RBTnodeCount(root->left)+RBTnodeCount(root->right)+1;
	}
	else return 0;
}


/**
 * \brief Helper function for "RBTtoIncArray"
 * \param void** Iterator: pointer to address of a data array
 * \param const size_t element_memsize: size of a data element in bytes
 */
static inline void inc_ptr(void** Iterator,const size_t element_memsize){
	const size_t CHARSIZE = sizeof(unsigned char);
	
	int d = element_memsize / CHARSIZE;
	unsigned char *C_Iter = (unsigned char*)*Iterator;
	C_Iter += d;
	*Iterator = C_Iter;
}


/**
 * \brief Extracts all the data of a red-black tree to an ascending-ordered array
 * \param rbNode* root: pointer to root node of the red-black tree
 * \param void** Array: pointer to address where the data array is stored
 * \param int* size: the size of the resulting array is stored under this pointer 
 */
void RBTtoIncArray(rbNode* root,void** Array,int* size){	
	*size = RBTnodeCount(root);
	*Array = (void*)realloc(*Array,(*size)*data_size_in_bytes);
	void* ptr = *Array;
	rbNode* it = RBTminNode(root);
	while(it!=NULL){
		memcpy(ptr,it->data,data_size_in_bytes);
		inc_ptr(&ptr,data_size_in_bytes);
		it = RBTsuccessor(it);
	}	
}

//////////////////////////////// Matrix functions /////////////////////////////////////////

static void print_vector(double* V,int n)__attribute__((unused));
static double CSRmatrixNorm(CSR_matrix* A)__attribute__((unused));


/**
 * \brief Scalar product 
 * \param double* x1: pointer to first vector
 * \param double* x2: pointer to second vector
 * \param int dim: dimension of the two vectors
 * \return double: the scalar product <x1,x2>
 * 
 * This function computes the euklidean scalar product p=<x1,x2> of 
 * two vectors x1 and x2.
 */
static double scalar(double* x1,double* x2,int dim){
	int i;
	double res = 0;
	for (i=0;i<dim;i++) res += x1[i]*x2[i];
	return res;
}



/**
 * \brief Compare function for index ordering
 * \param void* x1: pointer to the first matrix_element
 * \param void* x2: pointer to the second matrix_element
 * \return int: compare result
 * 
 * This function compares two matrix elements A_{ij} and A_{kl} 
 * according to their indices. The ersult is +1 if (i>k) or (i==k & j>l)
 * and -1 if (i<k) or (i==k & j<l). When (i==k & j==l) the function 
 * returns zero
 */
static int cmp_matr_ele(void* X1,void* X2){
	int i1 = ((matrix_element*)X1)->row_index;
	int j1 = ((matrix_element*)X1)->col_index;
	int i2 = ((matrix_element*)X2)->row_index;
	int j2 = ((matrix_element*)X2)->col_index;
	
	if (i1==i2){
		if (j1==j2) return 0; else return (j1>j2) ? 1 : -1;
	}
	else return (i1>i2) ? 1 : -1;	
}


/** 
 * \brief Frees memory of a matrix element
 * \param void* X: pointer to a matrix element
 */
static void free_matr_ele(void* X){
	matrix_element* a = (matrix_element*)X;
	free(a);
}


/**
 * \brief Compares eigenvalues of an eigensystem (for ordering)
 * \param void* X1: pointer to the first eigenvalue (cast from type eigen*)
 * \param void* X2: pointer to the second eigenvalue (cast from type eigen*)
 * \return int: compare result
 * 
 * This function compares two elements of an eigensystem for ordering. 
 * It returns +1 if eigenval(X1)>eigenval(X2) and -1 in the opposite case. 
 * If they equal within some tolerance it returns zero if the corresponding 
 * eigenvectors are linear dependent and else +1.
 */
static int cmp_eigen(void* X1,void* X2){
	const double tol = 1e-10;
	
	eigen* E1 = (eigen*)X1;
	eigen* E2 = (eigen*)X2;
	
	double d1 = E1->value;
	double d2 = E2->value;
	
	if (fabs(d1-d2)<tol){
		int dim = E1->dim;
		double* v1 = E1->vector;
		double* v2 = E2->vector;
		double s = scalar(v1,v2,dim);
		double c = scalar(v1,v1,dim)*scalar(v2,v2,dim)-s*s;
		return (fabs(c)<tol) ? 0 : 1;
	}
	else return (d1>d2) ? 1 : -1;	
}


/**
 * \brief Frees an element of an eigensystem
 * \param void* X: pointer to an eigenvalue (cast from type eigen*)
 */
static void free_eigen(void* X){
	eigen* E = (eigen*)X;
	free(E);
}


/**
 * \brief Create function for a real eigenystem element
 * \param double val: real eigenvalue
 * \param double* vec: pointer to the eigenvector
 * \param int dim: dimension of the eigenvector
 * 
 * This function creates an eigensystem element
 */
static eigen* create_eigen(double val,double* vec,int dim){
	eigen* res = (eigen*)malloc(sizeof(eigen));
	res->dim = dim;
	res->value = val;
	res->vector = vec;
	return res;
}


/**
 * \brief Allocates a list of integers initialized by zero
 * \param int n: number of entries of the list
 * \return int* : pointer to allocated list
 */
static int* zero_int_list(int n){													
	return (int*)calloc(n,sizeof(int));
}


/**
 * \brief Allocates a list of doubles initialized by zero
 * \param int n: dimension of vector
 * \return double* : pointer to allocated list
 */
static double* zero_vector(int n){												
	return (double*)calloc(n,sizeof(double));
}


/**
 * \brief Scales a vector (or dense matrix in linear meory) by some factor
 * \param double* x: pointer to a vector(/dense matrix)
 * \param double a: scaling factor
 * \param int n: dimension of vector (or total storage size of dense matrix)
 */
static void scale(double* x,double a,int n){
	int i;
	for (i=0;i<n;i++) x[i] *= a;
}


/**
 * \brief Computes a row-sum norm of a dense matrix
 * \param double* A: pointer to an array of matrix elements stored in linear row format
 * \param int rows: number of rows
 * \param int cols: number of cols
 * \return double: row-sum norm
 * 
 * The norm takes a row-ordered array (ptr(i,j) = i*cols+j) representing the matrix as first 
 * parameter. It computes all the row-sums  s_i=\sum_{j}|A_{ij}| and returns the largest. 
 */
static double denseMatrixNorm(double* A,int rows,int cols){
	int i,j;	
	double sum;
	
	double max = 0;
	for (i=0;i<rows;i++){
		sum = 0;
		for (j=0;j<cols;j++) sum += fabs(A[i*cols+j]);
		if (sum>max) max = sum;		
	}
	return max;
}


/**
 * \brief Normalizes a vector 
 * \param double* x: pointer to vector
 * \param int dim: dimension of vector
 * 
 * This function normalizes a vector with respect to the 
 * Euklidean norm: x-> x/sqrt(<x,x>)
 */
static void normalize(double* x,int dim){
	int i;
	double sum = 0;
	for (i=0;i<dim;i++) sum += x[i]*x[i];
	sum = sqrt(sum);
	if (sum!=0) for (i=0;i<dim;i++) x[i] /= sum;
}


/**
 * \brief Adds a multiple of a vector to another vector
 * \param double* x1: pointer to the vector that is modified
 * \param double* x2: pointer to the vector that is added
 * \param double a: scaling factor
 * \param int dim: dimenison of the vectors
 * 
 * The function computes x1 -> x1+a*x2
 */
static void addMultVec(double* x1,double* x2,double a,int dim){
	int i;
	for (i=0;i<dim;i++) x1[i] += a*x2[i];
}


/** 
 * \brief Dense identity matrix
 * \param int dim: dimension of vector space for the identety matrix
 * \return double: pointer to the linear row-indexed mememory
 * 
 * This function creates a dense version of a dim x dim identity matrix 
 * and returns a pointer to the newly created memory.
 */
static double* denseIDmatrix(int dim){
	int i;
	double* res = (double*)calloc(dim*dim,sizeof(double));
	for (i=0;i<dim;i++) res[i*(dim+1)] = 1.;
	return res;
}


/** 
 * \brief Adds a multiple of a dense matrix to another dense matrix
 * \param double* A: pointer to the dense matrix (row-ordered) that is modified 
 * \param double* B: pointer to the dense matrix (row-ordered) that is added
 * \param int rows: row number of the matrices
 * \param int cols: column number of the matrices
 * \param double factor: scaling factor
 * 
 * The function computes A -> A+factor*B. It does the same thing a the function 
 * "addMultVec" but is defined separately for clearity
 */
static void denseMatrixAdd(double* A,double* B,int rows,int cols,double factor){
	int k;
	for (k=0;k<rows*cols;k++) A[k] += factor*B[k];	
}


/**
 * \brief Frees all memory of a CSR-matrix struct and sets to NULL
 * \param CSR_matrix** A: pointer to the pointer where the matrix is stored 
 */
static void free_CSR_matrix(CSR_matrix** const A){
	if (*A!=NULL){
		if ((*A)->row_start!=NULL) free((*A)->row_start);
		if ((*A)->indices!=NULL) free((*A)->indices);
		if ((*A)->elements!=NULL) free((*A)->elements);
		free(*A);
		*A = NULL;
	}
}


/**
 * \brief Frees all memory of a CSC-matrix struct and sets to NULL
 * \param CSC_matrix** A: pointer to the pointer where the matrix is stored 
 */
static void free_CSC_matrix(CSC_matrix** const A){
	if (*A!=NULL){
		if ((*A)->col_start!=NULL) free((*A)->col_start);
		if ((*A)->indices!=NULL) free((*A)->indices);
		if ((*A)->elements!=NULL) free((*A)->elements);
		free(*A);
		*A = NULL;
	}
}


/**
 * \brief Initiates a CSR matrix
 * \param int row: number of rows
 * \param int init_size: initial number of elements that are allocated
 * \return CSR_matrix* : pointer to the allocated CSR matrix
 * 
 * The function allocates memory for storing <init_size> elements in a CSR matrix 
 * with <rows> rows. All row pointers are set to zero. 
 * 
 */
static CSR_matrix* CSRinit(int rows,int init_size){
	CSR_matrix* res = (CSR_matrix*)malloc(sizeof(CSR_matrix));
	res->rows = rows;
	res->row_start = (int*)calloc(rows+1,sizeof(int));
	res->indices = (int*)malloc(init_size*sizeof(int));
	res->elements = (double*)malloc(init_size*sizeof(double));
	//res->row_start[0] = 0;
	return res;
}


/**
 * \brief Initiates a CSC matrix
 * \param int col: number of columns
 * \param int init_size: initial number of elements that are allocated
 * \return CSC_matrix* : pointer to the allocated CSC matrix
 * 
 * The function allocates memory for storing <init_size> elements in a CSC matrix 
 * with <cols> columns. All column pointers are set to zero. 
 * 
 */
static CSC_matrix* CSCinit(int cols,int init_size){
	CSC_matrix* res = (CSC_matrix*)malloc(sizeof(CSC_matrix));
	res->cols = cols;
	res->col_start = (int*)calloc(cols+1,sizeof(int));
	res->indices = (int*)malloc(init_size*sizeof(int));
	res->elements = (double*)malloc(init_size*sizeof(double));
	//res->col_start[0] = 0;
	return res;
}


/**
 * \brief Resizes the element and index buffers of a CSR matrix
 * \param CSR_matrix* A: pointer to CSR matrix to be resized
 * \param int new_size: new buffer size
 */
static void CSRresizeBuffer(CSR_matrix* A,int new_size){
	A->indices = (int*)realloc(A->indices,new_size*sizeof(int));
	A->elements = (double*)realloc(A->elements,new_size*sizeof(double));	
}


/**
 * \brief Resizes the element and index buffers of a CSC matrix
 * \param CSC_matrix* A: pointer to CSC matrix to be resized
 * \param int new_size: new buffer size
 */
static void CSCresizeBuffer(CSC_matrix* A,int new_size){
	A->indices = (int*)realloc(A->indices,new_size*sizeof(int));
	A->elements = (double*)realloc(A->elements,new_size*sizeof(double));
}


/**
 * \brief Creates an identety matrix in CSR format
 * \param int dim: dimension of vector space of the matrix
 * \return CSR_matrix* : pointer to alloced CSR matrix
 * 
 * This function creates a dim x dim identity matrix in CSR format.
 */
static CSR_matrix* CSRid(int dim){
	int i;
	
	CSR_matrix* res = CSRinit(dim,dim);
	for (i=0;i<dim;i++){
		res->row_start[i] = i;
		res->indices[i] = i;
		res->elements[i] = 1.;
	}
	res->row_start[dim] = dim;
	return res;
}


/**
 * \brief Sets near-to-zero elements of a matrix to zero
 * \param CSR_matrix* A: pointer to a matrix in CSR format
 * \param double thres: threshold under which matrix elements are approsimated a zeros
 * 
 * Deletes all elaments with |A_{ij}|<thres frome the CSR matrix
 */
static void CSRchop(CSR_matrix* A,double thres){
	int i,j;
	double a;
	
	int rows = A->rows;
	int size = A->row_start[rows];
	int* start = (int*)malloc((rows+1)*sizeof(int));
	int* list = (int*)malloc(size*sizeof(int));
	double* eles = (double*)malloc(size*sizeof(double));
		
	int ptr = 0;
	for (i=0;i<rows;i++){
		start[i] = ptr;
		for (j=A->row_start[i];j<A->row_start[i+1];j++){
			a = A->elements[j];
			if (fabs(a)>thres){	
				list[ptr] = A->indices[j];
				eles[ptr] = a;
				ptr++;
			}					
		}		
	}
	start[rows] = ptr;
	
	free(A->row_start);
	free(A->indices);
	free(A->elements);
	A->row_start = start;
	A->indices = list;
	A->elements = eles;
	if (verbose_level>2) printf("reduction: %f\n",(double)ptr/size);	
}


/**
 * \brief Computes the row-sum norm of a CSR matrix
 * \param CSR_matrix* A: pointer a matrix in CSR format
 * \return double: row-sum norm
 * 
 * The function computes all the row-sums s_i=\sum_{j}|A_{ij}| and returns the largest. 
 */
static double CSRmatrixNorm(CSR_matrix* A){
	if (A!=NULL){
		int i,j,k;
		double sum;
		
		double max = 0;
		int n = A->rows;
		for (i=0;i<n;i++){
			sum = 0;
			for (j=A->row_start[i];j<A->row_start[i+1];j++){
				k = A->indices[j];
				sum += fabs(A->elements[k]);
			}
			if (sum>max) max = sum;
		}
		return max;
	}
	return 0;
}


/** 
 * \brief Prints diagonal elements of a CSR matrix to the output stream
 * \param CSR_matrix* A: pointer to matrix in CSR format
 * 
 * The function prints all the diagonal elements (comma separated) to the output stream. 
 */
static void print_CSR_diagnoals(CSR_matrix* A){
	if (A!=NULL){
		int i,j;
		int n = A->rows;
		for (i=0;i<n;i++){
			for (j=A->row_start[i];j<A->row_start[i+1];j++) if (A->indices[j]==i){
				if (i<n-1) printf("%f, ",A->elements[j]); else printf("%f\n",A->elements[j]);
			}
		}
		
	}
}


/**
 * \brief Checks if a CSR matrix is diagonal 
 * \param CSR_matrix* A: pointer to a matrix in CSR format
 * \param double tol: tolerance for considering elements as zeros
 * \return int: test result (0 or 1)
 * 
 * The function checks if the CSR matrix is diagonal within tolerance <tol> 
 * and return 1 if so, and 0 else.
 */ 
static int CSRisDiagonal(CSR_matrix* A,double tol){
	int i,j;
	double sum,diag;
	
	for (i=0;i<A->rows;i++){
		sum = 0;
		diag = 0;
		for (j=A->row_start[i];j<A->row_start[i+1];j++){
			if (A->indices[j]==i) diag = fabs(A->elements[j]); else sum += fabs(A->elements[j]);
		}
		if (diag!=0 && (sum/diag>tol)) return 0;
	}
	return 1;
}


/**
 * \brief Gets all the diagonal elements of a CSR matrix
 * \param CSR_matrix* A: pointer to CSR matrix
 * \return double* : pointer to vector where the diagonal elements are stored
 * 
 * The function returns a vector of length <A->rows> where all the diagonal elements are stored.
 */ 
static double* CSRgetDiagonal(CSR_matrix* A){
	int i,j;
	
	int n = A->rows;
	double* res = (double*)calloc(n,sizeof(double));
	for (i=0;i<n;i++){
		for (j=A->row_start[i];j<A->row_start[i+1];j++) if (A->indices[j]==i) res[i] = A->elements[j];
	}
	return res;
}


/**
 * \brief Gets the matrix element with spefic index of a CSR matrix
 * \param CSR_matrix* A: pointer to a CSR matrix
 * \param int row: row index
 * \param int col: column index
 * \return double: value of matrix element
 * 
 * This function simply returns the matrix element A_{row,col} of a CSR matrix.
 */
static double CSRgetElement(CSR_matrix* A,int row,int col){
	if (row<A->rows){
		int j;
		
		double res = 0;
		for (j=A->row_start[row];j<A->row_start[row+1];j++){
			if (A->indices[j]==col){
				res = A->elements[j];
				break;
			}			
		}
		return res;
		
	}
	else return 0;
}


/**
 * \brief Creates a 2x2 submatrix of a CSR matrix restricted to the given indices
 * \param CSR_matrix* A: pointer to a matrix in CSR format
 * \param int i: index of first base vector
 * \param int j: index of second base vector
 * \return double* : pointer to dense row-ordered 2x2 submatrix 
 * 
 * This function creates a submatrix acting on the subspace span{e_i,e_j} from 
 * the given CSR matrix A (e_i,e_j base vectors).
 */
static double* getSubMatrix2D(CSR_matrix* A,int i,int j){
	int k,l;
	
	double* res = (double*)calloc(4,sizeof(double));
	for (l=A->row_start[i];l<A->row_start[i+1];l++){
		k = A->indices[l];
		if (k==i) res[0] = A->elements[l];
		if (k==j) res[1] = A->elements[l];
	}
	for (l=A->row_start[j];l<A->row_start[j+1];l++){
		k = A->indices[l];
		if (k==i) res[2] = A->elements[l];
		if (k==j) res[3] = A->elements[l];
	}
	return res;
}


/**
 * \brief Converts a matrix from CSR to dense row-ordered format
 * \param CSR_matrix* A: pointer to CSR matrix to be converted
 * \param int cols: number of columns of the matrix <A>. 
 * \return double* : pointer to converted dense matrix
 * 
 * The function creates a copy of a CSR matrix in dense row-ordered format.
 */
static double* CSRtoDense(CSR_matrix* A,int cols){
	int i,j,k;
	
	double* res = (double*)calloc(A->rows*cols,sizeof(double));
	for (i=0;i<A->rows;i++){
		for (j=A->row_start[i];j<A->row_start[i+1];j++){
			k = A->indices[j];
			res[i*cols+k] = A->elements[j];
		}
	}
	return res;
}


/**
 * \brief Computes the transpose of a real dense matrix
 * \param double* M_row: pointer to a row-ordered dense matrix
 * \param int rows: row number 
 * \param int cols: column number
 * \return double* : pointer to the created dense transpose
 * 
 * The function return the transpose of the matrix <M_row> also 
 * in row-ordered dense format. 
 */
static double* denseTranspose(double* M_row,int rows,int cols){
	int i,j;
	double* T_row = (double*)malloc(rows*cols*sizeof(double));
	for (i=0;i<rows;i++){
		for (j=0;j<cols;j++) T_row[j*rows+i] = M_row[i*cols+j];
	}
	return T_row;
}


/**
 * \brief Replace a dense square matrix by its transpose
 * \param double** M_row: pointer to the adress of a row-ordered dense matrix
 * \param int n: row number of the nxn matrix
 */
static void DenseTranspose(double** M_row,int n){
	double* M = denseTranspose(*M_row,n,n);
	free(*M_row);
	*M_row = M;
}


/**
 * \brief Converts a CSR matrix to dense column-ordered format
 * \param CSR_matrix* A: pointer to a CSR matrix
 * \param int cols: number of columns of matrix <A>
 * \return double* : pointer to the converted matrix
 * 
 * This function creates a copy of the input CSR matrix converted to 
 * column-ordered dense format.
 */
static double* CSRtoDenseCol(CSR_matrix* A,int cols){
	double* AT = CSRtoDense(A,cols);
	double* res = denseTranspose(AT,A->rows,cols);
	free(AT);
	return res;
}


/**
 * \brief Converts a matrix in dense row-ordered format to CSR format
 * \param double* A: pointer to the row-ordered dense matrix
 * \param int col: number of columns of matrix <A>
 * \param int rows: number of rows of matrix <A>
 * \param double tol: tolerance for considering matrix elements as zeros
 * \return CSR_matrix* : pointer to converted matrix
 * 
 * The function creates a copy of a dense row-ordered matrix in CSR format.
 * Elements with absolute value smaller than <tol> are considered as zeros. 
 */
static CSR_matrix* denseToCSR(double* A,int cols,int rows,double tol){
	int i,j;
	double a;
	
	CSR_matrix* res = CSRinit(rows,rows*cols);
	int ptr = 0;
	int ptrA = 0;
	for (i=0;i<rows;i++){
		res->row_start[i] = ptr;
		for (j=0;j<cols;j++){
			a = A[ptrA++];
			if (fabs(a)>tol){
				res->indices[ptr] = j;
				res->elements[ptr] = a;
				ptr++;
			}
		}
	}
	res->row_start[rows] = ptr;
	CSRresizeBuffer(res,ptr);
	return res;
}


/**
 * \brief Converts a matrix in dense row-ordered format to CSC format
 * \param double* A: pointer to the row-ordered dense matrix
 * \param int col: number of columns of matrix <A>
 * \param int rows: number of rows of matrix <A>
 * \param double tol: tolerance for considering matrix elements as zeros
 * \return CSC_matrix* : pointer to converted matrix
 * 
 * The function creates a copy of a dense row-ordered matrix in CSC format.
 * Elements with absolute value smaller than <tol> are considered as zeros. 
 */
static CSC_matrix* denseToCSC(double* A,int cols,int rows,double tol){
	int i,j;
	double a;
	
	CSC_matrix* res = CSCinit(cols,rows*cols);
	int ptr = 0;
	for (j=0;j<cols;j++){
		res->col_start[j] = ptr;
		for (i=0;i<rows;i++){
			a = A[i*cols+j];
			if (fabs(a)>tol){
				res->indices[ptr] = i;
				res->elements[ptr] = a;
				ptr++;
			}
		}
	}
	res->col_start[cols] = ptr;
	CSCresizeBuffer(res,ptr);
	return res;	
}


/**
 * \brief Converts a matrix in dense column-ordered format to CSC format
 * \param double* A: pointer to the column-ordered dense matrix
 * \param int col: number of columns of matrix <A>
 * \param int rows: number of rows of matrix <A>
 * \param double tol: tolerance for considering matrix elements as zeros
 * \return CSC_matrix* : pointer to converted matrix
 * 
 * The function creates a copy of a dense column-ordered matrix in CSC format.
 * Elements with absolute value smaller than <tol> are considered as zeros. 
 */
static CSR_matrix* denseColToCSC(double* A_col,int cols,int rows,double tol){
	int i,j;
	double a;
	
	CSR_matrix* res = CSRinit(rows,rows*cols);
	int ptr = 0;
	for (j=0;j<cols;j++){
		res->row_start[j] = ptr;
		for (i=0;i<rows;i++){
			a = A_col[i*cols+j];
			if (fabs(a)>tol){
				res->indices[ptr] = i;
				res->elements[ptr] = a;
				ptr++;
			}
		}
	}
	res->row_start[rows] = ptr;
	CSRresizeBuffer(res,ptr);
	return res;	
}


/**
 * \brief Extracts a column vector from a CSR matrix
 * \param CSR_matrix* A: pointer to CSR matrix
 * \param int colIndex: index of the column to be extracted
 * \return double* : pointer to the created column vector
 * 
 * The function returns the column vector of index <colIndex> 
 * from the CSR matrix <A>.
 */
static double* CSRextractCol(CSR_matrix* A,int colIndex){
	int i,j,k;
	
	double* res = (double*)calloc(A->rows,sizeof(double));
	for (i=0;i<A->rows;i++){
		for (j=A->row_start[i];j<A->row_start[i+1];j++){
			k = A->indices[j];
			if (k==colIndex){
				res[i] = A->elements[j];
				break;
			}
		}
	}
	return res;	
}


/**
 * \brief Gets a sub block matrix from a CSR matrix
 * \param CSR_matrix* S: pointer to CSR matrix the block is extracted from
 * \param int start: smallest index of the block
 * \param int end: largest index of the block
 * \return CSR_matrix* : pointer to the created sub block
 * 
 * The function extracts a block with row indices i=start,start+1,...,end-1 and 
 * column index i=start,start+1,...,end-1 from the input and creates a new CSR matrix 
 * from the block.
 */
static CSR_matrix* CSRgetSubBlock(CSR_matrix* S,int start,int end){
	int i,j,k;
	
	int n = S->rows;
	int N = S->row_start[n];
	CSR_matrix* A = CSRinit(end-start,N);
	int ptrA = 0;
	for (i=start;i<end;i++){
		A->row_start[i-start] = ptrA;
		for (j=S->row_start[i];j<S->row_start[i+1];j++){
			k = S->indices[j];			
			if (k>=start && k<end){
				A->indices[ptrA] = k-start;
				A->elements[ptrA] = S->elements[j];
				ptrA++;
			}
		}
	}
	A->row_start[end-start] = ptrA;
	CSRresizeBuffer(A,ptrA);
	return A;
}


/**
 * \brief Inserts a square block into a square CSR matrix
 * \param CSR_matrix** S: pointer to address of the CSR matrix, in which the matrix <A> is inserted
 * \param CSR_matrix* A: pointer to the CSR matrix that is inserted into <S>
 * \param int start: insertion index
 * 
 * The function expands the (nxn)-matrix <S> by inserting the (mxm)-matrix A. The result is 
 * a (n+m x n+m)-matrix. For the new matrix it is 
 * 		Snew_{i,j}=Sold_{i,j} for 0<=i<start,0<=j<start
 * 		Snew_{start+i,start+j} = A_{i,j} for 0<=i<m,0<=j<m.
 * 		Snew_{m+i,m+j}=Sold_{i,j} for start<=i<n,start<=j<n
 * 		Snew_{i,j}=0 else. 
 */
void CSRinsertSubBlock(CSR_matrix** S,CSR_matrix* A,int start){
	int i,j;
	
	int n = (*S)->rows;
	int N = (*S)->row_start[n];
	int m = A->rows;
	int M = A->row_start[m];
	
	if (start+m>n){
		printf("error in %s: could not insert sub block -> abort\n",__func__);
		exit(EXIT_FAILURE);
	}
	
	int ptrS = (*S)->row_start[start];
	int bsize = N-ptrS;
	int* iBuffer = (int*)malloc(bsize*sizeof(int));
	double* dBuffer = (double*)malloc(bsize*sizeof(double));
	memcpy(iBuffer,&((*S)->indices[ptrS]),bsize*sizeof(int));
	memcpy(dBuffer,&((*S)->elements[ptrS]),bsize*sizeof(double));
	
	CSRresizeBuffer(*S,N+M);
	for (i=0;i<m;i++){
		(*S)->row_start[start+i] = ptrS;
		for (j=A->row_start[i];j<A->row_start[i+1];j++){
			(*S)->indices[ptrS] = A->indices[j]+start;
			(*S)->elements[ptrS] = A->elements[j];
			ptrS++;
		}
	}
	
	memcpy(&((*S)->indices[ptrS]),iBuffer,bsize*sizeof(int));
	memcpy(&((*S)->elements[ptrS]),dBuffer,bsize*sizeof(double));
	double diff = ptrS-(*S)->row_start[start+m];
	for (i=start+m;i<=n;i++) (*S)->row_start[i] += diff;	
	CSRresizeBuffer(*S,ptrS+bsize);
	if ((*S)->row_start[n]!=ptrS+bsize){
		printf("error in %s ->abort\n",__func__);
		exit(EXIT_FAILURE);
	}
	
	free(iBuffer);
	free(dBuffer);
}


/**
 * \brief Dense matrix product
 * \param double* A: first factor in dense row-ordered format
 * \param int rowsA: row number of matrix <A>
 * \param int colsA: column number of matrix <A>
 * \param double* B: second factor in dense row-ordered format
 * \param int rowsB: row number of matrix <B>
 * \param int colsB: column number of matrix <B>
 * \return double* : pointer to product <A>.<B> that is dense and row-ordered
 * 
 * The function computes the dense matrix product C=A.B.
 */
static double* denseMatrixMult(double* A,int rowsA,int colsA,double* B,int rowsB,int colsB){
	int i,j,k;
	double sum;
	
	if (colsA!=rowsB) return NULL;
	
	double* C = (double*)malloc(rowsA*colsB*sizeof(double));
	for (i=0;i<colsB;i++){
		for (j=0;j<rowsA;j++){
			sum = 0;			
			for (k=0;k<rowsB;k++) sum += A[j*colsA+k]*B[k*colsB+i];			
			C[j*colsB+i] = sum;
		}		
	}
	return C;
}


/**
 * \brief Adds a CSR matrix to another CSR matrix
 * \param CSR_matrix* A: pointer to the CSR matrix that is modified
 * \param CSR_matrix* B: pointer to the CSR matrix that is added
 * \param int cols: the column number of both matrices (that must be equal)
 * \param double factor: scaling factor
 * 
 * Thos function computes the sum A->A+factor*B of CSR matrices and stores 
 * it back in A.
 */
static void CSRmatrixAdd(CSR_matrix* A,CSR_matrix* B,int cols,double factor){
	int i,j,k,knext,lA,lB,size;
	
	int rows = A->rows;
	int* start = (int*)malloc((rows+1)*sizeof(int));
	int* list = (int*)malloc(rows*cols*sizeof(int));
	double* eles = (double*)malloc(rows*cols*sizeof(double));
	
	size = 0;
	for (i=0;i<rows;i++){		
		start[i] = size;
		k = A->row_start[i];
		knext = A->row_start[i+1];
		for (j=B->row_start[i];j<B->row_start[i+1];j++){
			lB = B->indices[j];									
			while(k<knext && (lA=A->indices[k])<lB){
				list[size] = lA;
				eles[size] = A->elements[k++];
				size++;
			}
			list[size] = lB;
			eles[size] = factor*B->elements[j];
			if (k<knext && (lA=A->indices[k])==lB) eles[size] += A->elements[k++];
			size++;
		}
		while(k<knext){
			list[size] = A->indices[k];
			eles[size] = A->elements[k];
			size++;
			k++;
		}		
	}
	start[rows] = size;
	
	list = (int*)realloc(list,size*sizeof(int));
	eles = (double*)realloc(eles,size*sizeof(double));
	free(A->row_start);
	free(A->indices);
	free(A->elements);
	A->row_start = start;
	A->indices = list;
	A->elements = eles;
}


/**
 * \brief Gets a list of indices that are decoupled
 * \param CSR_matrix* A: pointer to matrix in CSR format
 * \param int** list: pointer ot the address where the ersulting index list is stored
 * \param int* len: address where the length of the created list is stored
 * \param double tol: tolerance for considering elements as zeros
 * 
 * The function generates a list of indices that are decoupled from the rest of the matrix <A>. 
 * An index i is considered decoupled if \sum_{j}|A_{ij}|<tol and \sum_{j}|A_{ji}|<tol.
 */
static void CSRgetZeroCouplings(CSR_matrix* A,int** list,int* len,double tol){
	int i,j,k;
	double sum;
	
	int n = A->rows;
	*len = 0;
	*list = (int*)malloc(n*sizeof(int));
	for (i=0;i<n-1;i++){
		sum = 0;
		for (j=A->row_start[i];j<A->row_start[i+1];j++) if ((k=A->indices[j])!=i) sum += fabs(A->elements[j]);
		if (sum<tol){
			(*list)[*len] = i;
			(*len)++;
		}
	}
	*list = (int*)realloc(*list,(*len)*sizeof(int));
}


/**
 * \brief Helper function for "matrix_product_to_CSR"
 * \param CSR_matrix* A: first factor of matrix product
 * \param CSC_matrix* B: second factor of matrix product
 * \param int i: row index of matrix A
 * \param int j: column index of matrix B
 * 
 * This helper function computes the row x column product, in order to 
 * compute a matrix product. 
 */
static double CSRxCSC_row_times_col(CSR_matrix* A,CSC_matrix* B,int i,int j){
	
	int start_i = A->row_start[i];
	int start_j = B->col_start[j];
	int end_i = A->row_start[i+1];
	int end_j = B->col_start[j+1];
	if (start_i==end_i || start_j==end_j) return 0;
	
	int* indPtrA = &(A->indices[start_i]);
	int* indPtrB = &(B->indices[start_j]);
	int* endPtrA = &(A->indices[end_i-1]);
	int* endPtrB = &(B->indices[end_j-1]);	
	
	if ((*indPtrA)>(*endPtrB) || (*endPtrA)<(*indPtrB)) return 0;
	
	int k;
	
	int l = 0;
	double sum = 0;
	double* ptrA = &(A->elements[start_i]);
	double* ptrB = &(B->elements[start_j]);
	
	while (indPtrA<=endPtrA){
		k = *indPtrA;
		while(indPtrB<=endPtrB && (l=(*indPtrB))<k){			
			indPtrB++;
			ptrB++;
		}
		if (k==l){
			sum += (*ptrA)*(*ptrB);		
			indPtrB++;
			ptrB++;
		}
		indPtrA++;		
		ptrA++;
	}
	
	return sum;
}


/**
 * \brief Computes the matrix product
 * \param CSR_matrix* A: first factor of matrix product
 * \param CSC_matrix* B: second factor of matrix product
 * \return CSR_matrix*: pointer to the newly allocated matrix product P=A.B
 * 
 * This method computes the matrix product P of the (n x k)-matrix A and the (k x m)-matrix B. The 
 * result is the (n x m)-matrix P=A.B. For faster computation the matrix A is given in row format and 
 * B in column format. 
 */
static CSR_matrix* parallel_sqr_matrix_product_to_CSR(CSR_matrix* A,CSC_matrix* B){	
	const double tol = 1e-8;
				
	int i,j,start;
	double sum;
	int* ptrI,*ptrS;
	double* ptr;
	
	int n = A->rows;
	int m = B->cols;	
	
	
	double* buffer = (double*)malloc(n*m*sizeof(double));
	int* ind = (int*)malloc(n*m*sizeof(int));
	int* sizes = (int*)calloc(n,sizeof(int));
	
	#pragma omp parallel for schedule(dynamic,16) private(i,j,sum,ptr,ptrI,ptrS)
	for (i=0;i<n;i++){		
		ptr = &(buffer[i*m]);
		ptrI = &(ind[i*m]);
		ptrS = &(sizes[i]);
		for (j=0;j<m;j++){		
			sum = CSRxCSC_row_times_col(A,B,i,j);
			if (fabs(sum)>tol){
				*(ptr++) = sum;
				*(ptrI++) = j;
				(*ptrS)++;
			}		
		}
	}	
	
	int N = 0;
	for (i=0;i<n;i++) N += sizes[i];
	CSR_matrix* P = CSRinit(n,N);
	for (i=0;i<n;i++){
		start = P->row_start[i];
		P->row_start[i+1] = start+sizes[i];
		memcpy(&(P->indices[start]),&(ind[i*m]),sizes[i]*sizeof(int));
		memcpy(&(P->elements[start]),&(buffer[i*m]),sizes[i]*sizeof(double));
		
	}
	free(buffer);
	free(ind);
	free(sizes);
	
	return P;
}

/**
 * \brief Converts a CSR matrix to a dense matrix that is memory aligned
 * \param CSR_matrix* A: pointer to matrix in CSR format
 * \param int cols: number of columns of the matrix
 * \param double** denseA: pointer to address where the result of the conversion is stored
 * \param int* pitch: number of elements of each row in memory
 * 
 * This function converts a CSR matrix to dense row-ordered format. The allocated memory 
 * is chache-optimized for 64bit-line read and memory aligned. Therefore, the actual 
 * number of doubles in each stored row, called pitch, may be larger than the number of columns. 
 * The result is computed and stored under adress pitch.
 */
int CSRtoDensePadded(CSR_matrix* A,int cols,double** denseA,int* pitch){
	const size_t CACHE_LINE_SIZE = 64;
	const int chunk = CACHE_LINE_SIZE/sizeof(double);
	
	int i,j,err;
	
	int n = A->rows;
	int div = cols/chunk;
	int rem = cols-div*chunk;
	if (rem==0) *pitch = cols; else *pitch = (div+1)*chunk;
	size_t mem_size = (*pitch)*n*sizeof(double);
	err = posix_memalign((void**)denseA,CACHE_LINE_SIZE,mem_size);
	if (err){
		*denseA = NULL;
		return err;
	}
	memset(*denseA,0,mem_size);
	
	for (i=0;i<n;i++){
		double* ptr = &((*denseA)[i*(*pitch)]);
		for (j=A->row_start[i];j<A->row_start[i+1];j++) ptr[A->indices[j]] = A->elements[j];
	}
	return 0;
}


/**
 * \brief Converts a CSC matrix to a dense matrix that is memory aligned
 * \param CSC_matrix* A: pointer to matrix in CSR format
 * \param int rows: number of columns of the matrix
 * \param double** denseA: pointer to address where the result of the conversion is stored
 * \param int* pitch: number of elements of each column in memory
 * \return int: 0 if allocation succeeded, else: one of the error codes EINVAL and ENOMEM
 * 
 * This function converts a CSC matrix to dense column-ordered format. The allocated memory 
 * is chache-optimized for 64bit-line read and memory aligned. Therefore, the actual 
 * number of doubles in each stored column, called pitch, may be larger than the number of rows. 
 * The result is computed and stored under adress pitch.
 */
int CSCtoDenseColPadded(CSC_matrix* A,int rows,double** denseA,int* pitch){
	const size_t CACHE_LINE_SIZE = 64;
	const int chunk = CACHE_LINE_SIZE/sizeof(double);
	
	int i,j,err;
	
	int n = A->cols;
	int div = rows/chunk;
	int rem = rows-div*chunk;
	if (rem==0) *pitch = rows; else *pitch = (div+1)*chunk;
	size_t mem_size = (*pitch)*n*sizeof(double);
	err = posix_memalign((void**)denseA,CACHE_LINE_SIZE,mem_size);
	if (err){
		*denseA = NULL;
		return err;
	}
	memset(*denseA,0,mem_size);
	
	for (i=0;i<n;i++){
		double* ptr = &((*denseA)[i*(*pitch)]);
		for (j=A->col_start[i];j<A->col_start[i+1];j++) ptr[A->indices[j]] = A->elements[j];
	}
	return 0;
}


/** 
 * \brief Computes the dense matrix product
 * \param CSR_matrix* A: first factor in CSR format
 * \param int colsA: number of columns of <A>
 * \param CSC_matrix* B: second factor in CSC format
 * \param int rowsB: number of rows of <B>
 * \param double tol: tolerance to consider matrix elements of the product to be zero
 * \return CSR_matrix: product of the matrices in CSR format
 * 
 * The function computes the matrix product and is suited for small and highly filled matrices. 
 * The input matrices are converted to dense matrices to compute the matrx-matrix product 
 * for faster computation. The result is converted to CSR format afterwards. 
 */
static CSR_matrix* dense_matrix_product_to_CSR(CSR_matrix* A,int colsA,CSC_matrix* B,int rowsB,double tol){	
	int i,j,k,pitchA,pitchB;
	double *denseA,*denseB,*denseC;
	
	int rowsA = A->rows;
	int colsB = B->cols;
	CSRtoDensePadded(A,colsA,&denseA,&pitchA);
	CSCtoDenseColPadded(B,rowsB,&denseB,&pitchB);
	if (colsA!=rowsB || pitchA!=pitchB){
		printf("error in function %s: column number (%d) of first factor does not match row number of second factor (%d) -> abort\n",__func__,colsA,rowsB);
		exit(EXIT_FAILURE);
	}
	
	denseC = (double*)malloc(rowsA*colsB*sizeof(double));

	//#pragma omp parallel for schedule(static) private(i,j,k) firstprivate(rowsA,colsB)
	for (i=0;i<rowsA;i++){
		for (j=0;j<colsB;j++){
			double sum = 0;
			double* itA = &(denseA[pitchA*i]);
			double* itB = &(denseB[pitchB*j]);
			for (k=0;k<colsA;k++) sum += (*(itA++))*(*(itB++));
			denseC[i*rowsA+j] = sum;
		}
	}
	
	free(denseA);
	free(denseB);
	
	CSR_matrix* C = denseToCSR(denseC,colsB,rowsA,tol);
	free(denseC);
	return C;
}


/**
 * \brief Computes the product of two square matrices
 * \param CSR_matrix* A: first factor in CSR format
 * \param CSC_matrix* B: second factor in CSC format
 * \param double tol: tolerance to consider matrix elements of the product to be zero
 * \return CSR_matrix: product of the matrices in CSR format
 * 
 * This function choses a suitable method for a square matrix product computation depending on the 
 * size of the matrix. 
 */
static CSR_matrix* CSRxCSC_SqrMatrixProduct(CSR_matrix* A,CSC_matrix* B,double tol){
	const double dim_thres = 20;
	
	int n = A->rows;
	if (n>dim_thres){
		CSR_matrix* C = parallel_sqr_matrix_product_to_CSR(A,B);
		CSRchop(C,tol);
		return C;
	}
	else return dense_matrix_product_to_CSR(A,n,B,n,tol);
}


/**
 * \brief Right-multiplies a CSR matrix with a CSC matrix
 * \param CSR_matrix** A: pointer to address of the matrix product in CSR format
 * \param CSC_matrix* B: pointer to a CSC matrix
 * \param double tol: tolerance to consider matrix elements of the product to be zero
 * 
 * This function replaces the original matrix <A> by the product <A>.<B> where <A> 
 * is given in CSR and <B> in CSC format.
 */
static void CSRsqrMatrixRightMult(CSR_matrix** A,CSC_matrix* B,double tol){	
	CSR_matrix* P = CSRxCSC_SqrMatrixProduct(*A,B,tol);
	free_CSR_matrix(A);
	*A = P;
}


/**
 * \brief Computes the transpose of matrix
 * \param CSR_matrix* A: pointer to input matrix in CSR format
 * \param int cols: number of columns of the matrix
 * \return CSR_matrix* : transpoed of the input in CSR format
 * 
 * This function returns the transpose of a real square matrix in CSR format
 */
static CSR_matrix* CSRtranspose(CSR_matrix* A,int cols){
	int i,j;
	matrix_element* ele;
	
	RBTree_set_compare(&cmp_matr_ele);
	RBTree_set_free(&free_matr_ele);
	RBTree_set_data_size(sizeof(matrix_element));
	
	int rows = A->rows;
	int size = A->row_start[rows];
	CSR_matrix* TA = CSRinit(cols,size);
	
	rbNode* root = NULL;
	for (i=0;i<rows;i++){
		for (j=A->row_start[i];j<A->row_start[i+1];j++){
			ele = (matrix_element*)malloc(sizeof(matrix_element));
			ele->col_index = i;
			ele->row_index = A->indices[j];
			ele->value = A->elements[j];
			if (RBTinsertElement(&root,ele)!=NULL) free(ele);
		}
	}
	
	int counter = 0;
	int prev_start = -1;		
	TA->row_start[0] = 0;
	rbNode* node = RBTminNode(root);
	while(node!=NULL){
		ele = (matrix_element*)node->data;
		TA->indices[counter] = ele->col_index;
		TA->elements[counter] = ele->value;
		if (ele->row_index!=prev_start){						
			prev_start = ele->row_index;
			TA->row_start[prev_start] = counter;
		}					
		node = RBTsuccessor(node);
		counter++;					
	}
	TA->row_start[cols] = counter;	
	if (counter!=size){
		printf("error in %s!\n",__func__);
		exit(EXIT_FAILURE);
	}
	
	RBTfree(root);
	return TA;
}



/**
 * \brief Converts a CSR matrix to a CSC matrix
 * \param CSR_matrix* A: pointer to a matrix in CSR format
 * \param int cols: number of columns of the matrix
 * \return CSC_matrix* : the converted in CSC format
 * 
 * The function creates a copy of the input CSR matrix in CSC format.
 * 
 */
static CSC_matrix* CSRtoCSC(CSR_matrix* A,int cols){
	CSR_matrix* TA = CSRtranspose(A,cols);
	CSC_matrix* res = (CSC_matrix*)malloc(sizeof(CSC_matrix));
	res->cols = cols;
	res->col_start = TA->row_start;
	res->indices = TA->indices;
	res->elements = TA->elements;
	free(TA);
	return res;
}



/**
 * \brief Creates a permutation map that moves decoupled indices to the end
 * \param CSR_matrix* H: pointer to input matrix in CSR format
 * \param double tol: tolerance to consider matrix elaments a zeros
 * \param int** map: pointer to address where the index permutation map is stored
 * \param int** inv: poitner to address where the inverse map is stored
 * \return int: index where the decoupled block of the permuted matrix begins
 * 
 * This function looks for decoupled indices and creates a permutation map that 
 * moves the corresponding elemnts to the last part of the matrix. The inverse map 
 * is created, too. The index where the decoupled block starts is returned
 */
static int getPermMap(CSR_matrix* H,double tol,int** map,int** inv){
	int len = 0;
	int* list = NULL;
	CSRgetZeroCouplings(H,&list,&len,tol);
	if (len>0){
		int i,j,k,found;
		
		int n = H->rows;
		*map = (int*)malloc(n*sizeof(int));		
		*inv = (int*)malloc(n*sizeof(int));		
				
		k = 0;
		for (i=0;i<n;i++){
			found = 0;
			for (j=0;j<len;j++) if (list[j]==i) found = 1;
			if (!found){
				(*inv)[k] = i;
				(*map)[i] = k;
				k++;
			}
		}
		if (k+len!=n){
			printf("error in %s -> abort\n",__func__);
			exit(EXIT_FAILURE);
		}
		for (j=0;j<len;j++){
			(*inv)[k+j] = list[j];
			(*map)[list[j]] = k+j;
		}
		return n-len;		
	}
	else{
		*map = NULL;
		*inv = NULL;
		return -1;
	}
}



/**
 * \brief Applies a permutation map of a square matrix
 * \param CSR_matrix** A: pointer to address of the matrix to be permuted
 * \param int* P: poitner to index-permutation map as list of length <A->rows>
 * 
 * The function applies an index permutation to the input square CSR matrix.
 */
static void CSRindexPermutation(CSR_matrix** A,int* P){
	int i,j;
	matrix_element* ele;
	
	RBTree_set_compare(&cmp_matr_ele);
	RBTree_set_free(&free_matr_ele);
	RBTree_set_data_size(sizeof(matrix_element));
	
	int n = (*A)->rows;
	int N = (*A)->row_start[n];
	
	rbNode* root = NULL;
	for (i=0;i<n;i++){		
		for (j=(*A)->row_start[i];j<(*A)->row_start[i+1];j++){
			ele = (matrix_element*)malloc(sizeof(matrix_element));
			ele->row_index = P[i];
			ele->col_index = P[(*A)->indices[j]];
			ele->value = (*A)->elements[j];
			if (RBTinsertElement(&root,ele)!=NULL) free(ele);
		}
	}
	
	free_CSR_matrix(A);
	*A = CSRinit(n,N);
	
	int ptrA = 0;	
	rbNode* iterator = RBTminNode(root);
	i = -1 ;
	while(iterator!=NULL){
		ele = (matrix_element*)iterator->data;
		if (ele->row_index>i){
			i = ele->row_index;
			(*A)->row_start[i] = ptrA;
		}
		(*A)->indices[ptrA] = ele->col_index;
		(*A)->elements[ptrA] = ele->value;
		iterator = RBTsuccessor(iterator);
		ptrA++;
	}
	(*A)->row_start[n] = ptrA;
	
	RBTfree(root);
	
	if (ptrA!=N){
		printf("error in %s -> abort\n",__func__);
		exit(EXIT_FAILURE);
	}	
}

 
 /**
 * \brief Right-multiplies a dense square matrix to another
 * \param double** A: pointer to address of the dense matrix product in row-ordered format
 * \param double* B: pointer to dense matrix in column-ordered format
 * \param int dim: dimension of the vector spaces the matrices act on
 * 
 * This function replaces the original dense square matrix <A> by the product <A>.<B> where <A> 
 * is given in row-order and <B> in columns order. The product is row-ordered again.
 */
static void denseQuadMatrixMult(double* A,double* B,int dim){
	double* P = denseMatrixMult(A,dim,dim,B,dim,dim);
	memcpy(A,P,dim*dim*sizeof(double));	
	free(P);
}


/**
 * \brief Computes the Wilkinson shift
 * \param CSR_matrix* H: pointer to input matrix in CSR format
 * \param int index: diagonal index of the 2x2 block for which the Wilkinson shift is computed
 * \return double: Wilkinson shift
 * 
 * The function computes the Wilkinson shift of a 2x2 block at diagonal index <index>.
 */
static double getWilkinsonShift(CSR_matrix* H,int index){
	double res;
	
	double im = 0;
	double re1 = 0;
	double re2 = 0;
	
	double* Hsub = getSubMatrix2D(H,index-1,index);
	getEigen2D(Hsub,&re1,&re2,&im);
	
	if (fabs(re1-Hsub[3])<fabs(re2-Hsub[3])) res = re1; else res = re2;
	free(Hsub);
	return res;
}


/**
 * \brief Computes a vector to generate a Householder transformation
 * \param double* M_row: square matrix to apply the transformation on in row-ordered dense format
 * \param int dim: dimension of the vector space the matrix acts on.
 * \param int col: column index for which the Housholder reflection is computed
 * \return double* : pointer where the result of length <dim> is stored 
 * 
 * In order to create a Householder transformation T = 1-2*v.v^T for input matrix <M_row> with corresponding column 
 * index <col>, the vector v is computed.
 */
static double* getHouseholderVector(double* M_row,int dim,int col){
	int j;
	double a,b;
	
	double sum = 0;
	for (j=col+1;j<dim;j++){
		a = M_row[j*dim+col];
		sum += a*a;
	}
	
	double* v = (double*)calloc(dim,sizeof(double));
	if (sum==0) return v;
	else{
		sum = sqrt(sum);
		b = M_row[(col+1)*dim+col];
		for (j=col+1;j<dim;j++) v[j] = M_row[j*dim+col];
		v[col+1] += (b>=0) ? sum: -sum;
		normalize(v,dim);
		return v;
	}
}


/** 
 * \brief Right-multiplies with a Householder matrix
 * \param double* M_row: row-ordered dense square matrix that is right-multiplied
 * \param double* v: vector from which the Householder transformation T(v)=1-2*v.v^T is build.
 * \param int dim: dimension of the vector space <M_row> acts on
 * 
 * The function right-multiplies the dense square matrix <M_row> by a Householder 
 * transformation T(v)=1-2*v.v^T created from vector <v>. The result is stored under the same address.  
 */
static void rightMultHouseholder(double* M_row,double* v,int dim){
	
	int i;
	double s;
	double* CPM = (double*)malloc(dim*dim*sizeof(double));
	memcpy(CPM,M_row,dim*dim*sizeof(double));
	
	for (i=0;i<dim;i++){		
		s = scalar(&(CPM[i*dim]),v,dim);		
		addMultVec(&(M_row[i*dim]),v,-2.*s,dim);
	}	
	free(CPM);
}


/** 
 * \brief Left-multiplies with a Householder matrix
 * \param double* M_row: row-ordered dense square matrix that is left-multiplied
 * \param double* v: vector from which the Householder transformation T(v)=1-2*v.v^T is build.
 * \param int dim: dimension of the vector space <M_row> acts on
 * 
 * The function left-multiplies the dense square matrix <M_row> by a Householder 
 * transformation T(v)=1-2*v.v^T created from vector <v>. The result is stored under the same address.  
 */
static void leftMultHouseholder(double* M_row,double* v,int dim){
	int i,j;
	double s;
	double* CPM = (double*)malloc(dim*dim*sizeof(double));
	memcpy(CPM,M_row,dim*dim*sizeof(double));
	
	for (i=0;i<dim;i++){
		s = 0;
		for (j=0;j<dim;j++) s += CPM[j*dim+i]*v[j];
		s *= 2.;
		for (j=0;j<dim;j++) M_row[j*dim+i] -= v[j]*s;
	}
	free(CPM);
}


/**
 * \brief Hessenberg decomposition
 * \param double* M_row: dense square matrix in row-orderd format that shall be decomposed
 * \param int dim: dimension of the vector space <M_row> acts on
 * \param double** H: pointer to address under which the resulting upper-Hessenberg matrix is stored
 * \param double** Q: pointer to address under which the orthogonal transformation matrix is stored
 * 
 * The function generates a Hessenberg decomposition of the input matrix <M_row> that is 
 * M_row = Q.H.Q^T, where H is an upper-Hessenber matrix and Q is orthogonal (Q^T.Q=1).
 */
static void denseHessenbergDecomposition(double* M_row,int dim,double** H,double** Q){
	int i;
	
	if (*Q==NULL) *Q = denseIDmatrix(dim);
	if (*H==NULL){
		*H = (double*)malloc(dim*dim*sizeof(double));
		memcpy(*H,M_row,dim*dim*sizeof(double));	
	}
	
	for (i=0;i<dim-1;i++){
		double* v = getHouseholderVector(*H,dim,i);
		rightMultHouseholder(*H,v,dim);
		leftMultHouseholder(*H,v,dim);
		rightMultHouseholder(*Q,v,dim);				
		free(v);
	}
}


/**
 * \brief Givens rotation
 * \param double* M_row: input square matrix in dense row-ordered format
 * \param int dim: dimension of the vector space <M_row> is acting on.
 * \param int i: subspace index one
 * \param int j: subspace index two
 * \param int c: cosine of rotation angle
 * \param inc s: sine of rotation angle (s^2+c^2=1)
 * 
 * Applies a Givens rotation to the 2D-subspace spanned by indices <i> and <j> 
 * of the the row-ordered dense input matrix <M_row>.
 * -> multiplies by {{c,-s},{s,c}} in subspace span(e_i,e_j)
 */
static void denseGivensRotation(double* M_row,int dim,int i,int j,double c,double s){		
	double x,y;
	
	int k = 0;
	double* pi = &(M_row[dim*i]);
	double* pj = &(M_row[dim*j]);
	while(k<dim){
		x = *pi;
		y = *pj;		
		*(pi++) = c*x-s*y;
		*(pj++) = s*x+c*y;
		k++;
	}
}


/**
 * \brief QR-decomposition of an upper-Hessenberg matrix
 * \param CSR_matrix* H: pointer to an upper-Hessenberg matrix in CSR format
 * \param int dim: dimension of the vector space <H> acts on
 * \param CSC_matrix** Q: pointer to address where the orthogonal factor of the QR-decomposition is stored in CSC format
 * \param CSR_format** R: pointer to address where the triangular factor of the QR-decomposition is stored in CSR format
 * \param CSR_format** T: this is optional. If a NULL pointer is stored under the given address nothing is done. If there is a pointer to a CSR matrix, it is left-multiplied by <Q>
 * \param double tol: tolerance for considering matrix elements as zeros.
 * 
 * This function creates a QR-decomposition of the given upper-Hessenberg matrix <H>: 
 * H = Q.R, where Q is orthogonal and R triangular. One may optinally provide a pointer <*T>
 * to a transformation matrix that is left multiplied by Q if one wishes to accumulate several transformations 
 * in one matrix. 
 */
static void denseHessenbergQR(CSR_matrix* H,int dim,CSC_matrix** Q,CSR_matrix** R,CSR_matrix** T,double tol){					
	int i,j;																		
	double c,s,r,aij,ajj;
	
	double* denseR = CSRtoDense(H,dim);
	double* denseQ = denseIDmatrix(dim);
	double* denseT = (*T!=NULL) ? CSRtoDenseCol(*T,dim) : NULL;
	
	for (j=0;j<dim-1;j++){	
		i = j+1;
		ajj = denseR[j*(dim+1)];
		aij = denseR[i*dim+j];					
		if (aij!=0){
			r = sqrt(ajj*ajj+aij*aij);	
			c = ajj/r;
			s = aij/r;
			denseGivensRotation(denseR,dim,i,j,c,s);						// left multiplication, R,Q row-ordered
			denseGivensRotation(denseQ,dim,i,j,c,s);						
			if (denseT!=NULL) denseGivensRotation(denseT,dim,i,j,c,s);		// right multiplication, since T column-ordered
		}
	}
	
	DenseTranspose(&denseQ,dim);
	
	free_CSC_matrix(Q);
	free_CSR_matrix(R);
	*Q = denseToCSC(denseQ,dim,dim,tol);
	*R = denseToCSR(denseR,dim,dim,tol);
	
	if (denseT!=NULL){
		free_CSR_matrix(T);
		*T = denseColToCSC(denseT,dim,dim,tol);
		free(denseT);
	}	
	
	free(denseQ);
	free(denseR);
}	


/**
 * \brief Checks if an upper-Hessenberg Matrix has approximately block-diagonal structure and returns the separating column indices
 * \param CSR_matrix* H: pointer to an upper-Hessenberg matrix in CSR format
 * \param double tol: threshold beneath wich matrix elemnts are consiedered as zero
 * \param int last_sep: column index that confines the check to submatrix {0,...,<last_sep>}
 * \param int** separators: pointer to an address where the address of the newly created array of separator indices is stored.
 * \param int* num: address where the length of the created separator list is stored
 * 
 */
static int HessenbergHasBlockStructure(CSR_matrix* H,double tol,int last_sep,int** separators,int* num){
	int i,j,k,l,coupled,start,end;	
	double sub;
		
	int n = H->rows;	
	int N = H->row_start[n];
	int m = (last_sep>=0) ? last_sep : n;
	int size = 1;
	int* sep = (int*)malloc(n*sizeof(int));
	sep[0] = 0;
	for (i=1;i<m;i++) if ((j=H->row_start[i])<N && j<H->row_start[i+1]){
		sub = (H->indices[j]==i-1) ? H->elements[j] : 0;
		coupled = 1;
		if (fabs(sub)<tol){			
			coupled = 0;
			end = 0;
			for (k=sep[size-1];k<i;k++){
				start = end;
				end = H->row_start[k+1];
				if (end-start>0 && H->indices[end-1]>i){
					for (l=start;l<end;l++) if (H->indices[l]>i && fabs(H->elements[l])>100.*tol){
						coupled = 1;
						goto jmp_next;
					}
				}
				
			}			
			for (k=i;k<m;k++){				
				start = end;
				end = H->row_start[k+1];
				if (end-start>0 && H->indices[start]<i){
					for (l=start;l<end;l++) if (H->indices[l]<i && fabs(H->elements[l])>tol){
						coupled = 1;
						goto jmp_next;
					}
				}
				
			}					
		}
jmp_next:if (!coupled){
			if (separators==NULL){
				free(sep);
				return i;
			}
			sep[size++] = i;
		}
	}
	else{
		coupled = 0;
		goto jmp_next;
	}
	
	if (separators==NULL){
		free(sep);
		return -1;
	}
	else{
		if (last_sep>=0 && last_sep<n) sep[size++] = last_sep;
		sep[size++] = n;
		*separators = (int*)realloc(sep,size*sizeof(int));
		*num = size;
		return (size>2) ? (*separators)[size-1]: -1;
	}
}


/**
 * \brief Helper function for "QRiterations"
 * \param CSR_matrix** H: pointer to the address of the upper-Hessenberg matrix whose eigensystem is requested
 * \param int offset: column-index offset with respect to the larger supermatrix
 * \param int* eigen_ind: list of columns whose eigenvalues are found
 * \param int* eig_num: pointer to the length of the list <eigen_ind>
 * \param double tol: tolerance beneath wich matrix elements are considered as zeros
 * \param double chop_thres: cutoff threshold beneath which matrix element are deleted after matrix multiplication
 * \param int max_iter: maximum number of QR iterations within some subdiagonal should approach zero (< <tol>)
 * \param shiftStrategy st: shift strategy to apply
 * 
 * For detaila see function description of "QRiterations".
 */
static CSR_matrix* subQRiterations(CSR_matrix** H,int offset,int* eigen_ind,int* eig_num,double tol,double chop_thres,int max_iter,shiftStrategy st){	
	int i,sep;
	
	int dim = (*H)->rows;
	if (dim==1){
		double h = (*H)->elements!=NULL ? (*H)->elements[0]: 0;
		eigen_ind[(*eig_num)++] = offset;//h;
		if (verbose_level>1) printf("dim: %d finished! 0 iterations, %d eigenvalues\n",dim,*eig_num);
		if (verbose_level>2) printf("dim: 1 finished! 0 iterations\teigenvalue: %f\n",h);
		return CSRid(1);		
	}
	else if (dim==2){
		double re[2];
		double im = 0;
		double* denseH = CSRtoDense(*H,2);			
		getEigen2D(denseH,&(re[0]),&(re[1]),&im);
		double abs = (fabs(re[0])>fabs(re[1])) ? sqrt(re[0]*re[0]+im*im) : sqrt(re[1]*re[1]+im*im);
		if (abs>tol && im/abs>1000.*tol){
			printf("Error in %s: complex eigenvalues encoutered in submatrix: %f+%fi -> abort\n",__func__,re[0],im);
			exit(EXIT_FAILURE);
		}
		double* denseQ = (double*)malloc(4*sizeof(double));
		getRealEigenVectors2D(denseH,re,denseQ,tol);
		
		denseH[0] = re[0];
		denseH[1] = 0;
		denseH[2] = 0;
		denseH[3] = re[1];
		free_CSR_matrix(H);
		*H = denseToCSR(denseH,2,2,tol);
		free(denseH);		
		CSR_matrix* Q = denseToCSR(denseQ,2,2,tol);
		free(denseQ);
		
		eigen_ind[(*eig_num)++] = offset;
		if (eigen_ind[*eig_num]>=0) eigen_ind[(*eig_num)++] = offset+1;//re[1];
		if (verbose_level>1) printf("dim: 2 finished! 0 iterations, %d eigenvalues\n",*eig_num);
		if (verbose_level>2) printf("dim: 2 finished! 1 iteration\teigenvalues: %f, %f\n",re[0],re[1]);
		return Q;
	}
	else{		
		int len = 0;
		int* separators = NULL;
		int* pivot = NULL;
		int* inv_pivot = NULL;
		int pivot_sep = getPermMap(*H,tol,&pivot,&inv_pivot);
		if (pivot!=NULL){			
			CSRindexPermutation(H,pivot);
			free(pivot);
		}
		
		if (HessenbergHasBlockStructure(*H,100.*tol,pivot_sep,&separators,&len)>=0){		
			int tra_offset;
						
			CSR_matrix** Hsub = (CSR_matrix**)malloc((len-1)*sizeof(CSR_matrix*));
			CSR_matrix** Qsub =  (CSR_matrix**)malloc((len-1)*sizeof(CSR_matrix*));
			for (i=len-1;i>0;i--){				
				Hsub[i-1] = CSRgetSubBlock(*H,separators[i-1],separators[i]);
				tra_offset = (inv_pivot==NULL) ? offset+separators[i-1] : offset+inv_pivot[separators[i-1]];
				if (eigen_ind[*eig_num]>=0) Qsub[i-1] = subQRiterations(&(Hsub[i-1]),tra_offset,eigen_ind,eig_num,tol,chop_thres,max_iter,st);
				else Qsub[i-1] = CSRid(separators[i]-separators[i-1]);
			}
		
			free_CSR_matrix(H);
			*H = CSRinit(dim,0);
			CSR_matrix* Q = CSRinit(dim,0);
			for (i=0;i<len-1;i++){
				CSRinsertSubBlock(H,Hsub[i],separators[i]);
				free_CSR_matrix(&(Hsub[i]));
				CSRinsertSubBlock(&Q,Qsub[i],separators[i]);
				free_CSR_matrix(&(Qsub[i]));			
			}
			
			if (inv_pivot!=NULL){
				CSRindexPermutation(H,inv_pivot);
				CSRindexPermutation(&Q,inv_pivot);
				free(inv_pivot);
			}
							
			free(Hsub);
			free(Qsub);
			if (separators!=NULL) free(separators);
			return Q;
		}
		else{			
			CSC_matrix* Q = NULL;
			CSR_matrix* R = NULL;
			CSR_matrix* ID = CSRid(dim);
			CSR_matrix* Qtotal = CSRid(dim);
			
			double shift = 0;
			switch(st){
				case ZERO: shift = 0;break;
				case LASTDIAG: shift = CSRgetElement(*H,dim-1,dim-1);break;
				case WILKINSON: shift = getWilkinsonShift(*H,dim-1);break;				
			}			
			
			int diagFlag = 0;
			i = 0;			
			while((sep = HessenbergHasBlockStructure(*H,100.*tol,-1,NULL,NULL))<0 && !diagFlag && i<max_iter){
				CSRmatrixAdd(*H,ID,dim,-shift);		
				denseHessenbergQR(*H,dim,&Q,&R,&Qtotal,chop_thres);				
				CSR_matrix* QTHQ = CSRxCSC_SqrMatrixProduct(R,Q,chop_thres);
				CSRmatrixAdd(QTHQ,ID,dim,shift);				
				free_CSR_matrix(H);
				*H = QTHQ;				
			
				CSRchop(Qtotal,chop_thres);
				CSRchop(*H,chop_thres);
				diagFlag = CSRisDiagonal(*H,tol);										
				i++;				
			}
	
			if (i==max_iter){
				printf("Error in %s: no convergence of algorithm within %d iterations, dimension: %d, shift: %f\n",__func__,i,dim,shift);		
				exit(EXIT_FAILURE);
			}
			else{
				if (diagFlag){
					double* eig = CSRgetDiagonal(*H);
					i = 0;
					while(i<dim && eigen_ind[*eig_num]>=0){
						eigen_ind[(*eig_num)++] = (inv_pivot==NULL) ? offset+i : offset+inv_pivot[i];
						i++;
					}										
					if (verbose_level>2){
						printf("dim: %d finished! %d iterations\teigenvalues:",dim,i);
						print_CSR_diagnoals(*H);
					}
					free(eig);
				}
				if (verbose_level>1) printf("dim: %d finished! %d iterations, %d eigenvalues\n",dim,i,*eig_num);							
			}	
	
			free_CSR_matrix(&ID);
			free_CSR_matrix(&R);
			free_CSC_matrix(&Q);
			
			if (!diagFlag && eigen_ind[*eig_num]>=0){				
				CSR_matrix* Qrow = subQRiterations(H,offset,eigen_ind,eig_num,tol,chop_thres,max_iter,st);
				Q = CSRtoCSC(Qrow,dim);			
				CSRsqrMatrixRightMult(&Qtotal,Q,chop_thres);							
				free_CSC_matrix(&Q);
				free_CSR_matrix(&Qrow);
			}
			
			if (inv_pivot!=NULL){
				CSRindexPermutation(H,inv_pivot);
				CSRindexPermutation(&Qtotal,inv_pivot);
				free(inv_pivot);
			}		
			
			if (separators!=NULL) free(separators);
			return Qtotal;	
		}
		
	}
}


/**
 * \brief Get the directory where the current process is executed
 * \param char* buffer: string buffer where the result is stored
 * \param int size: size of the buffer
 * \return int: 1 when succeeded, else 0
 */
static int getExeDir(char* buffer,int size){	
	char proc[BUFF_SIZE];
	sprintf(proc,"/proc/%d/exe",getpid());
	memset(buffer,0,size*sizeof(char));
	int status = readlink(proc,buffer,size);	
	if (status>0){
		char* ptr = strrchr(buffer,'/');
		*ptr = '\0';
	}
	return status;
}


/**
 * \brief Creates a gnuplot script file to visualize data. 
 * \param char* cmdname: absolute filename of the output script file
 * \param const char* dataname: filename relative to the gnuplot script where the data is stored
 * \param int* colsToPlot: pointer to an array of length <colNum> where the indices of the columns to be plotted are listed
 * \param int colNum: number of columns to be plotted
 * \param char** labels: pointer to array of strings (of length <colNum>) whith plot labels
 * \return int: returns 1 when successful, else 0
 * 
 * This method creates a gnuplot script file that plots the data in the textfile <dataname> that is stored in columns. The 
 * column indices to be plotted must be given in the array <colsToPlot>. The scipt is output under filename <cmdname>.
 * Optionally one can give labels to each plotted column. If <labels>==NULL automatic labels or diplayed.
 */
static int createGnuplotFile(char* cmdname,const char* dataname,int* colsToPlot,int colNum,char** labels){
	char opts[BUFF_SIZE];
	
	FILE* file = fopen(cmdname,"w");
	if (file!=NULL){
		int i,index;
		
		fprintf(file,"plot ");
		for (i=0;i<colNum;i++){
			index = colsToPlot[i];
			if (labels!=NULL) sprintf(opts,"title %s w l",labels[i]); else sprintf(opts,"w l");
			if (i<colNum-1) fprintf(file,"'%s' u ($%d*$%d) %s,",dataname,index,index,opts); else fprintf(file,"'%s' u ($%d) %s\n",dataname,index,opts);
		}		
		fclose(file);
		return 1;		
	}
	else{
		printf("Warning: could not create file %s\n",cmdname);		
		return 0;
	}
}

/**
 * \brief Checks is the program named <progname> is in PATH via the shell command "which"
 * \param char* progname: Name of the executable to check
 * \return int: Returns 1 if the program is registered in PATH. If the executalbe does not exist or is not in PATH then 0 is returned
 */
static int programInPath(char* progname){
	const char* tempname = "inpathtemp";
	
	char cmd[BUFF_SIZE];
	char nam[BUFF_SIZE];
	
	sprintf(cmd,"which -a %s > %s 2> /dev/null",progname,tempname);
	system(cmd);
	getExeDir(nam,BUFF_SIZE);
	strcat(nam,"/");
	strcat(nam,tempname);
	FILE* file = fopen(nam,"r");
	if (file!=NULL){
		fseek(file,0,SEEK_SET);		
		long start = ftell(file);
		fseek(file,0,SEEK_END);
		long end = ftell(file);
		fclose(file);
		sprintf(cmd,"rm %s",nam);
		system(cmd);
		if (end==start) return 0; else return 1;
    }
    else return 0;	
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////* public functions */////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \brief Sets the verbose level of the computation
 * \param int level: detail level for terminam output
 * 
 * Depending on the verbose level more and more detailed output 
 * is written to the terminal. If <level>=0 the QR algorithm 
 * runs silently if no errors occur. For <level>=1,2 more details are written.
 */
void QRsetVerboseLevel(int level){
	verbose_level = level;
}


/**
 * \brief Computes the eigenvectors of a symmetric 2x2 matrix
 * \param double* denseH: The input matrix in dense-row format 
 * \param double* eig: if eigenvalues are known, input a pointer to an array of them. Else set to NULL.
 * \param double* denseQ: The resulting eigenvectors are stored in an allocated matrix in dense-column format
 * \param double tol: tolerance for off-diagonal elements to be treated as zeros
 * 
 * This method computes the (real) eigenvectors of a symmetric 2x2 matrix. The result is stored in <denseQ> 
 * that mus be a pointer to an allocated array of size=4*sizeof(double). When v1 and v2 are the eigenvectors, then 
 * v1=(Q[0],Q[1])^T and v2 = (Q[2],Q[3])^T.
 */
void getRealEigenVectors2D(double* denseH,double* eig,double* denseQ,double tol){
	double x,r,la1,la2;
	
	double h = denseMatrixNorm(denseH,2,2);
	double a12 = fabs(denseH[1]/h);
	double a21 = fabs(denseH[2]/h);
		
	if (a12>tol && a21>tol){
		
		if (eig==NULL){
			la1 = 0;
			la2 = 0;			
			getEigen2D(denseH,&la1,&la2,NULL);			
		}
		else{
			la1 = eig[0];
			la2 = eig[1];
		}
		
		x = (denseH[0]-la1)/denseH[1];
		r = sqrt(1.+x*x);
		denseQ[0] = 1./r;
		denseQ[1] = x/r;	
			
		x = (denseH[3]-la2)/denseH[2];
		r = sqrt(1.+x*x);
		denseQ[2] = x/r;
		denseQ[3] = 1./r;	
	}
	else if (a21>tol){
		x = (denseH[0]-denseH[2])/denseH[0];
		r = sqrt(1.+x*x);
		denseQ[0] = 1./r;
		denseQ[1] = x/r;	
			
		denseQ[2] = 0;
		denseQ[3] = 1.;	
	}
	else if (a12>tol){		
		x = (denseH[3]-denseH[1])/denseH[3];
		r = sqrt(1.+x*x);
		denseQ[2] = x/r;
		denseQ[3] = 1./r;	
			
		denseQ[0] = 1.;
		denseQ[1] = 0;	
	}
	else{
		denseQ[0] = 1.;
		denseQ[1] = 0;	
		denseQ[2] = 0;
		denseQ[3] = 1.;	
	}
}


/**
 * \brief Computes the eigenvalues of a real 2x2 matrix
 * \param double* A: matrix in dense-row format 
 * \param double* eig1_real: pointer to a double, where the real part of the smaller eigenvalues is stored
 * \param double* eig2_real: pointer to a double, where the real part of the larger eigenvalues is stored
 * \param double* eig_im: if not NULL the absolute imaginary part of both eigenvalues is stored at this address
 * 
 * This method computes the eigenvalue of a real (not necessary symmetric) matrix. Its eigenvalues 
 * are stored in the following manner: 
 * eigevalue#1 = <eig1_real>-i*<eig_im> and eigevalue#2 = <eig2_real>+i*<eig_im>
 */
void getEigen2D(double* A,double* eig1_real,double* eig2_real,double* eig_im){
	
	double tr = A[0]+A[3];
	double det = A[0]*A[3]-A[1]*A[2];
	double rad = tr*tr-4.*det;
	if (rad>=0){
		double s = sqrt(rad);
		double min = (tr-s)/2.;
		double max = (tr+s)/2.;		
		if (eig_im!=NULL) *eig_im = 0;
		if (fabs(A[0]-min)<fabs(A[0]-max)){
			*eig1_real = min;
			*eig2_real = max;
		}
		else{
			*eig1_real = max;
			*eig2_real = min;
		}		
	}
	else{
		if (eig_im!=NULL) *eig_im = sqrt(-rad);
		*eig1_real = tr/2.;
		*eig2_real = tr/2.;
	}
}


/**
 * \brief Computes the Rayleigh quotient
 * \param double* A: pointer to <dim>x<dim> square matrix in dense-row format.
 * \param double* x: pointer to a vector of length <dim>
 * \param int dim: dimension of space
 * \return double: the Rayleigh quotient
 * 
 * This method computes the Rayleigh quotient R of a vector <x> for a given matrix <A>
 * according to R = <x,A.x>/<x,x>, where <*,*> is a scalar product.
 */
double denseRayleighQuotient(double* A,double* x,int dim){
	int i,j;
	double row;
	
	double sum = 0;
	double div = 0;
	for (i=0;i<dim;i++){
		row = 0;
		for (j=0;j<dim;j++) row += A[i*dim+j]*x[j];
		sum += x[i]*row;
		div += x[i]*x[i];
	}
	if (div>0) return sum/div; else return 0;
}


/**
 * \brief QR algorithm to compute all eigenvales and vectors of a matrix
 * \param void* matrix: The matrix whose eigenvalues we are interested in. It must be stored in CSR_matrix format and its pointer cast to void* 
 * \param double* eigenValues: A pointer to an allocated array (of length <eig_num>) where the computed eigenvalues are stored
 * \param double** eigenVectors: A pointer to an allocated array of pointers (of length <eig_num>) where the computed eigenvectors are stored
 * \param int eig_num: The number of eigenvalues/vectors to be computed
 * \param double tol: floating point number tolarance. Beneath this threshold, matrix elements are approximated as zero entries
 * \param int max_iter: maximum number of QR iterations
 * \param shiftStrategy st: shift strategy for faster convergence
 * 
 * This method computes <eig_num> eigenvalues of a symmetric indefinite matrix <matrix> (where <eig_num> <= dimension(<matrix>)).
 * For this, the QR algorithm with shift and deflation is used. First the matrix is transformed to an upper-Hessenberg matrix with  
 * Householder transformations. Subsequently a number of QR iterations is applied until some subdiagonal elements approach zero. 
 * The matrix is split if it approaches a block-diagonal form. For these blocks the QR algorithm is applied recursively (deflation). 
 * To improve convergence the spectrum can be improved by using a shift the the matrix. Three options are available for <st>:
 *  (1) ZERO: zero shift (is useful if one is only interested in the smallest absolute eigenvalues)
 * 	(2) LASTDIAG: the shift is set to the last diagonal entry of the current matrix
 *	(3) WILKINSON: the Wilkinson shift is obtained by the eigenvalues of the last 2x2 block of the current Hessenber matrix. 
 * 		The eigenvalue that is closest to the last diagonal entry is chosen to be the shift. 
 */ 
void QRiterations(void* matrix,double* eigenValues,double** eigenVectors,int eig_num,double tol,int max_iter,shiftStrategy st){
	int i,k;
	
	CSR_matrix* A = (CSR_matrix*)matrix;
	int dim = A->rows;
	double* denseH = NULL;
	double* denseA = CSRtoDense(A,dim);
	double* denseQ = denseIDmatrix(dim);	
	denseHessenbergDecomposition(denseA,dim,&denseH,&denseQ);	
	CSR_matrix* Qtotal = denseToCSR(denseQ,dim,dim,tol);
	CSR_matrix* H = denseToCSR(denseH,dim,dim,tol);
	free(denseH);
	free(denseQ);
	free(denseA);
	
	int eig_ptr = 0;
	double chop_thres = 1e-3*tol;
	int* eig_ind = (int*)calloc(eig_num+1,sizeof(int));
	eig_ind[eig_num] = -1;	
	CSR_matrix* Q = subQRiterations(&H,0,eig_ind,&eig_ptr,tol,chop_thres,max_iter,st);
	CSC_matrix* Qcol = CSRtoCSC(Q,dim);
	CSRsqrMatrixRightMult(&Qtotal,Qcol,tol);	
	free_CSR_matrix(&Q);
	
	for (i=0;i<eig_ptr;i++){
		k = eig_ind[i];
		eigenValues[i] = CSRgetElement(H,k,k);
		eigenVectors[i] = CSRextractCol(Qtotal,k);
	}
	
	if (verbose_level>0){		
		denseA = CSRtoDense(A,dim);
		denseH = CSRtoDense(H,dim);
		denseQ = CSRtoDense(Qtotal,dim);
		double* denseQT = denseTranspose(denseQ,dim,dim);
				
		double* P = denseMatrixMult(denseQT,dim,dim,denseQ,dim,dim);
		double* I = denseIDmatrix(dim);
		denseMatrixAdd(P,I,dim,dim,-1.);
		
		denseQuadMatrixMult(denseQ,denseH,dim);
		denseQuadMatrixMult(denseQ,denseQT,dim);
		denseMatrixAdd(denseQ,denseA,dim,dim,-1.);
		printf("finished !\n\nrelative decomposition error: %e, orthognonality: %e\n",denseMatrixNorm(denseQ,dim,dim)/denseMatrixNorm(denseA,dim,dim),
		 denseMatrixNorm(P,dim,dim));			
		
		free(P);
		free(I);
		free(denseA);
		free(denseH);
		free(denseQ);
		free(denseQT);
	}
	
	free(eig_ind);
	free_CSC_matrix(&Qcol);
	free_CSR_matrix(&Qtotal);
	free_CSR_matrix(&H);		
}


/**
 * \brief Exemplary potential of the Schroedinger equation
 * \param double x: position in space
 * \return double: potential value
 */
static double pot(double x){
	if (fabs(x)<0.1) return -2500.; 
	else if (fabs(x-0.15)<0.05) return 2500.;
	else return 0;
}


/** 
 * \brief Test routine for the QR algorithm
 * \param int n: number of discretization nodes of the system
 * \param threads: number of threads to use
 * 
 * This method is a test routine for the QR algorithm. As an example it 
 * computes the energy levels of the one-dimensional stationary Schroedinger 
 * equation: (-d^2/dx^2+Vpot(x))psi = E*psi on interval [0,1] with a potential 
 * well in the center. A finite difference scheme with n discretization nodes is used.
 * When finished it outputs the all eigenvalues and eigenvectors 
 * to the terminal. If the system size is n>20 then the eigenvectors corresponding
 * to the five lowest energy levels are written to the file "eigendata.txt" in the
 * directory where the program is executed. If "gnuplot" is installed and in the PATH 
 * variable it is started and the eigenfunctions are shown (gnu script stored under "cmg.gnu"). 
 */
void testQR(int n,int threads){
	const double tol = 1e-10;
	const int max_iter = 200;
	const int plot_thres = 20;
	
	int i,j;
	double x;	
	
#if defined(_OPENMP)
	omp_set_dynamic(0);    
	omp_set_num_threads(threads);
#else
	threads = 1;
#endif

	QRsetVerboseLevel(0);
	
	// create matrix of 1D Schroedinger equation on domain [0,1] with n nodes
	double d = (double)1/(n-1);
	double* denseA = denseIDmatrix(n);
	scale(denseA,2.,n*n);
	for (i=0;i<n;i++){
		if (i>0) denseA[i*n+i-1] = -1.; else denseA[i*n+i] = 1.;
		if (i<n-1) denseA[i*n+i+1] = -1.; else denseA[i*n+i] = 1.;
	}
	scale(denseA,1./(d*d),n*n);
	for (i=0;i<n;i++){
		x = (double)i/(n-1);
		denseA[i*n+i] += pot(x-0.5);
	}
	
	// transform to sparse matrix format
	CSR_matrix* A = denseToCSR(denseA,n,n,tol);
	
	// QR algorithm
	double* eigenvalues = zero_vector(n);
	double** eigenvectors = (double**)malloc(n*sizeof(double*));
	printf("start QR algorithm ...\n");
#if defined(_OPENMP)	
	double tstart = omp_get_wtime();
#else 
	clock_t tstart = clock();
#endif	
	QRiterations(A,eigenvalues,eigenvectors,n,tol,max_iter,WILKINSON);
#if defined(_OPENMP)		
	double tend = omp_get_wtime();
	printf("\nexecution time with %d threads: %f s\n",threads,tend-tstart);	
#else 
	clock_t tend = clock();
	printf("\nexecution time with %d threads: %f s\n",threads,(float)(tend-tstart)/CLOCKS_PER_SEC);	
#endif	
	
	// print eigenvalues and eigenvectors
	eigen* lambda;
	printf("\neigenvalues (ascending order):\n");
	RBTree_set_compare(&cmp_eigen);
	RBTree_set_free(&free_eigen);
	RBTree_set_data_size(sizeof(eigen));
	rbNode* root = NULL;
	for(i=0;i<n;i++){		
		lambda = create_eigen(eigenvalues[i],eigenvectors[i],n);
		RBTinsertElement(&root,lambda);
	}
	
	// orders eigenvalues (and vectors) by its absolute values in descending order
	i = 0;
	rbNode* it = RBTminNode(root);
	while(it!=NULL && i<plot_thres){
		lambda = (eigen*)it->data;
		printf("%.5f,",lambda->value);
		it = RBTsuccessor(it);
		i++;
	}
	printf("...\n");
	printf("\neigenvectors:\n");
	if (n<plot_thres){		
		// print result to terminal 	
		it = RBTmaxNode(root);
		i = 1;
		while(it!=NULL){
			lambda = (eigen*)it->data;
			printf("#%d:\t(",i);
			for (j=0;j<n;j++) printf("%.4f,",lambda->vector[j]); 
			printf("\b)\n");
			it = RBTpredecessor(it);
			i++;
		}		
	}
	else{
		const char* dataFileName = "eigendata.txt";							// data filename 
		const char* cmdFileName = "cmd.gnu";								// gnuplot script filename
		const int m = 5;													// number of eigenvectors to be plotted
		
		char cmd[BUFF_SIZE];
		char datafile[BUFF_SIZE];
		char cmdfile[BUFF_SIZE];
		
		// print eigenvectors as columns, starting with closest to zero eigenvalue 
		printf("large output: print results to file %s\n\n",dataFileName);
		FILE* file = NULL;		
		if (getExeDir(datafile,BUFF_SIZE)>0){
			strcpy(cmdfile,datafile);
			strcat(datafile,"/");
			strcat(datafile,dataFileName);			
			
			file = fopen(datafile,"w");	
			if (file!=NULL){
				eigen* orderedEigs = NULL;
				RBTtoIncArray(root,(void**)&orderedEigs,&j);		
				for (i=0;i<n;i++){
					x = (double)i/(n-1);					
					for (j=0;j<m;j++){
						fprintf(file,"%e\t",orderedEigs[j].vector[i]);						
					}
					fprintf(file,"%e",1e-5*pot(x-0.5));			
					fprintf(file,"\n");
				}				
				fclose(file);
				
				// checks if gnuplot is available and if so, plots the first five eigenvectors
				if (!programInPath("gnuplot")) printf("no gnuplot installed or given in PATH -> skip diplay\n");
				else{
					int* cols = zero_int_list(m+1);
					for (i=0;i<m+1;i++) cols[i] = i+1;				
					strcat(cmdfile,"/");
					strcat(cmdfile,cmdFileName);	
					char** labels = (char**)malloc((m+1)*sizeof(char*));
					for (i=0;i<m+1;i++){
						labels[i] = (char*)malloc(BUFF_SIZE*sizeof(char));
						if (i<m) sprintf(labels[i],"'E%d=%.2f'",i,orderedEigs[i].value); else sprintf(labels[i],"'V_pot(x)' lc 9");				
					}								
					if (createGnuplotFile(cmdfile,dataFileName,cols,m+1,labels)){					
						sprintf(cmd,"gnuplot --persist %s",cmdFileName);
						system(cmd);
					}
					for (i=0;i<m+1;i++) free(labels[i]);
					free(labels);
					free(cols);				
				}
				free(orderedEigs);
			}
			else printf("Warning: could not create file %s\n",datafile);			
		}						
	
	}
	
	free(denseA);
	for (i=0;i<n;i++) free(eigenvectors[i]);		
	free(eigenvectors);
	RBTfree(root);	
	free(eigenvalues);
	free_CSR_matrix(&A);
}




