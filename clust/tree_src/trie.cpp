#include "trie.h"
#include <pybind11/pybind11.h>
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

#define PAD 5              // Position of padding nodes.
#define EOS -1             // End Of String, for 'dash()'.

// Translation table to insert nodes in the trie.
//          ' ': PAD (5)
//     'a', 'A': 1
//     'c', 'C': 2
//     'g', 'G': 3
//     't', 'T': 4
static const int translate[256] = { 
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 PAD,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
   0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
   0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
   0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
// Translation table to query a sequence in the trie.
// In the table below, non DNA letters are set to a numerical
// value of 6, which will always cause of mismatch with
// sequences translated from the table above. 'PAD' and '-' are
// the only non DNA symbols that match themselves.
//          ' ': PAD (5)
//          '-': 0
//     'a', 'A': 1
//     'c', 'C': 2
//     'g', 'G': 3
//     't', 'T': 4
static const int altranslate[256] = { 
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
 PAD,6,6,6,6,6,6,6,6,6,6,6,6,0,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,1,6,2,6,6,6,3,6,6,6,6,6,6,6,6,
   6,6,6,6,4,6,6,6,6,6,6,6,6,6,6,6,
   6,1,6,2,6,6,6,3,6,6,6,6,6,6,6,6,
   6,6,6,6,4,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
   6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
};

struct arg_t {
   char        tau;
   char        maxtau;
   int       * query;
   int         height;
   int         err;
};

void dash(node_t *, const int *, struct arg_t);
void destroy_from(node_t *, void (*)(void *), int, int, int);
int get_height(trie_t *);
void init_pebbles(node_t *);
node_t *insert(node_t *, int);
node_t *insert_wo_malloc(node_t *, int, node_t *);
node_t *new_trienode(void);
result_t *poucet(node_t *, int, struct arg_t);
int recursive_count_nodes(node_t *node, int, int);

// Globals.
int ERROR = 0;

int get_height(trie_t *trie) { return trie->info->height; }

// ------  SEARCH FUNCTIONS ------ //

result_t *search(
    trie_t *trie,
    const char *query,
    const int tau)
// SYNOPSIS:                                                              
//   Front end query of a trie with the "poucet search" algorithm. Search 
//   does not start from root. Instead, it starts from a given depth      
//   corresponding to the length of the prefix from the previous search.  
//   Since initial computations are identical for queries starting with   
//   the same prefix, the search can restart from there. Each call to     
//   the function starts ahead and seeds pebbles for the next query.      
//                                                                        
// PARAMETERS:                                                            
//   trie: the trie to query                                              
//   query: the query as an ascii string                                  
//   tau: the maximum edit distance                                       
//   hits: a hit stack to push the hits                                   
//   start_depth: the depth to start the search                           
//   seed_depth: how deep to seed pebbles                                 
//                                                                        
// RETURN:                                                                
//   A pointer to 'hits' node array.                                      
//                                                                        
// SIDE EFFECTS:                                                          
//   The non 'const' parameters of the function are modified. The node    
//   array can be resized (which is why the address of the pointers is    
//   passed as a parameter) and the trie is modified by the trailinng of  
//   effect of the search.                                                
{
   ERROR = 0;

   int height = get_height(trie);
   int length = strlen(query);


   // // Make sure the cache is allocated.
   // info_t *info = trie->info;

   // // Reset the pebbles that will be overwritten.
   // start_depth = max(start_depth, 0);
   // for (int i = start_depth+1 ; i <= min(seed_depth, height) ; i++) {
   //    info->pebbles[i]->nitems = 0;
   // }

   // Translate the query string. The first 'char' is kept to store
   // the length of the query, which shifts the array by 1 position.
   int translated[M];
   translated[0] = length;
   translated[length+1] = EOS;
   for (int i = 0 ; i < length ; i++) {
      translated[i+1] = altranslate[(int) query[i]];
   }

   // Set the search options.
   struct arg_t arg = {
        .tau     = tau,
        .query   = translated,
        .height  = height,
   };
    
   node_t *start_node = trie->root;
   return poucet(start_node, 1, arg);;
}

result_t * poucet(
    node_t *node,
    const int depth,
    struct arg_t arg)
// SYNOPSIS:                                                              
//   Back end recursive "poucet search" algorithm. Most of the time is    
//   spent in this function. The focus node sets the values of  an L-     
//   shaped section (but with the angle on the right side) of the dynamic 
//   programming table for its children. One of the arms of the L is      
//   identical for all the children and is calculated separately. The     
//   rest is classical dynamic programming computed in the 'cache' struct 
//   member of the children. The path of the last 8 nodes leading to the  
//   focus node is encoded by a 32 bit integer, which allows to perform   
//   dynamic programming without parent pointer.                          
//
//   All the leaves of the trie are at a the same depth called the   
//   "height", which is the depth at which the recursion is stopped to    
//   check for hits. If the maximum edit distance 'tau' is exceeded the   
//   search is interrupted. On the other hand, if the search has passed   
//   trailing depth and 'tau' is exactly reached, the search finishes by  
//   a 'dash()' which checks whether an exact suffix can be found.        
//   While trailing, the nodes are pushed in the g_stack 'pebbles' so     
//   they can serve as starting points for future searches.               
//
//   Since not all the sequences have the same length, they are prefixed
//   with the 'PAD' character (value 5, printed as white space) so that
//   the total length is equal to the height of the trie. This imposes
//   an important modification to the recursion, indicated by the label
//   "PAD exception" on two different lines of the code below. Without
//   this modification, "AAAAATA" and "AAAAA" would be aligned this way.
//
//                               AAAAATA
//                               x||||x|
//                               _AAAAAA
//
//   However, the best alignment is the following.
//
//                               AAAAATA
//                               |||||x|
//                               AAAAA_A
//
//   The solution is to initialize the alignment score to 0 whenever the
//   padding character is met, which in effect is equivalent to ignoring
//   starting the alignment after the PADs.
//                                                                        
// PARAMETERS:                                                            
//   node: the focus node in the trie                                     
//   depth: the depth of the children in the trie.                        
//                                                                        
// RETURN:                                                                
//   'void'.                                                              
//                                                                        
// SIDE EFFECTS:                                                          
//   Same as 'search()', it modifies nodes of the trie and the node       
//   array 'arg.hits' where hits are pushed.                              
{
   // If search faild, return -1, otherwise return the label on the leaf
   result_t * result = new result_t;
   // Reached height of the trie: it's a hit!
   if (node->isEnd) {
      result->label = node ->label;
      char *pcache = node->cache + TAU;
      result->distance = pcache[0];
      return result;
   }
   // This makes it easier to distinguish the part that goes upward,
   // with positive index and requiring the path, from the part that
   // goes horizontally, with negative index and requiring previous
   // characters of the query.
   char *pcache = node->cache + TAU;
   // Risk of overflow at depth lower than 'tau'.
   int maxa = min((depth-1), arg.tau);

   // Penalty for match/mismatch and insertion/deletion resepectively.
   unsigned char mmatch;
   unsigned char shift;

   // Part of the cache that is shared between all the children.
   char common[9] = {1,2,3,4,5,6,7,8,9};

   // The branch of the L that is identical among all children
   // is computed separately. It will be copied later.
   int32_t path = node->path;
   // Upper arm of the L (need the path).
   if (maxa > 0) {
      // Special initialization for first character. If the previous
      // character was a PAD, there is no cost to start the alignment.
      // This is the "PAD exeption" mentioned in the SYNOPSIS.
      mmatch = (arg.query[depth-1] == PAD ? 0 : pcache[maxa]) +
                  ((path >> 4*(maxa-1) & 15) != arg.query[depth]);
      shift = min(pcache[maxa-1], common[maxa]) + 1;
      common[maxa-1] = min(mmatch, shift);
      for (int a = maxa-1 ; a > 0 ; a--) {
         mmatch = pcache[a] + ((path >> 4*(a-1) & 15) != arg.query[depth]);
         shift = min(pcache[a-1], common[a]) + 1;
         common[a-1] = min(mmatch, shift);
      }
   }

   node_t *child;
   for (int i = 0 ; i < 6 ; i++) {
      // Skip if current node has no child at this position.
      if ((child = node->child[i]) == NULL) continue;

      // Same remark as for parent cache.
      char local_cache[] = {9,8,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,8,9};
      
      char *ccache = depth == arg.height ? local_cache + 9 : child->cache + TAU;
      memcpy(ccache+1, common, TAU * sizeof(char));


      // Horizontal arm of the L (need previous characters).
      if (maxa > 0) {
         // See comment above for initialization.
         // This is the "PAD exeption" mentioned in the SYNOPSIS.
         mmatch = ((path & 15) == PAD ? 0 : pcache[-maxa]) +
                     (i != arg.query[depth-maxa]);
         shift = min(pcache[1-maxa], maxa+1) + 1;
         ccache[-maxa] = min(mmatch, shift);
         for (int a = maxa-1 ; a > 0 ; a--) {
            mmatch = pcache[-a] + (i != arg.query[depth-a]);
            shift = min(pcache[1-a], ccache[-a-1]) + 1;
            ccache[-a] = min(mmatch, shift);
         }
      }
      // Center cell (need both arms to be computed).
      mmatch = pcache[0] + (i != arg.query[depth]);
      shift = min(ccache[-1], ccache[1]) + 1;
      ccache[0] = min(mmatch, shift);

      // Stop searching if 'tau' is exceeded.
      if (ccache[0] > arg.tau) continue;
      result_t * tmp_result = poucet(child, depth+1, arg);
      if(tmp_result->label != -1)return tmp_result;      
   }
   return result;
}

// ------  TRIE CONSTRUCTION AND DESTRUCTION  ------ //

trie_t *
new_trie(
    unsigned int height)
// SYNOPSIS:                                                              
//   Front end trie constructor.                                          
//                                                                        
// PARAMETERS:                                                            
//   height: the fixed depth of the leaves                                
//                                                                        
// RETURN:                                                                
//   A pointer to trie root with meta information and no children.        
{

   if (height < 1) {
      fprintf(stderr, "error: the minimum trie height is 1\n");
      ERROR = __LINE__;
      return NULL;
   }

   trie_t *trie = new trie_t;
   if (trie == NULL) {
      fprintf(stderr, "error: could not create trie\n");
      ERROR = __LINE__;
      return NULL;
   }

   node_t *root = new_trienode();
   if (root == NULL) {
      fprintf(stderr, "error: could not create root\n");
      ERROR = __LINE__;
      free(trie);
      return NULL;
   }

   info_t *info = new info_t;
   if (info == NULL) {
      fprintf(stderr, "error: could not create trie\n");
      ERROR = __LINE__;
      free(root);
      free(trie);
      return NULL;
   }

   // Set the values of the meta information.
   info->height = height;

   trie->root = root;
   trie->info = info;

   return trie;

}

node_t *
new_trienode(void)
// SYNOPSIS:                                                              
//   Back end constructor for a trie node. All values are initialized to  
//   null, except the cache for dynamic programming which is initialized  
//   as a root node.                                                      
//                                                                        
// RETURN:                                                                
//   A pointer to trie node with no data and no children.                 
{

   node_t *node = new node_t;
   if (node == NULL) {
      fprintf(stderr, "error: could not create trie node\n");
      ERROR = __LINE__;
      return NULL;
   }

   // Initialize the cache. This is important for the
   // dynamic programming algorithm.
   const char init[] = {8,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,8};
   memcpy(node->cache, init, 2*TAU+1);
   return node;

}

void
insert_string(
    trie_t *trie,
    const char *string,
    int label)
// SYNOPSIS:                                                              
//   Front end function to fill in a trie. Insert a string from root, or  
//   simply return the node at the end of the string path if it already   
//   exists.                                                              
//                                                                        
// RETURN:                                                                
//   The leaf node in case of succes, 'NULL' otherwise.                   
//                                                                        
// NB: This function is not used by 'starcode()'.
{

   int i;

   int nchar = strlen(string);
   if (nchar != get_height(trie)) {
      fprintf(stderr, "error: can only insert string of length %d\n",
            get_height(trie));
      ERROR = __LINE__;
      return;
   }
   
   // Find existing path.
   node_t *node = trie->root;
   for (i = 0; i < nchar - 1; i++)
   {
      node_t *child;
      int c = translate[(int)string[i]];
      if ((child = (node_t *)node->child[c]) == NULL)
      {
         break;
      }
      node = child;
   }

   // Append more nodes.
   for (; i < nchar - 1; i++)
   {
      int c = translate[(int)string[i]];
      node = insert(node, c);
      if (node == NULL)
      {
         fprintf(stderr, "error: could not insert string\n");
         ERROR = __LINE__;
         return;
      }
   }
   node->isEnd = true;
   node->label = label;
   return;
}

node_t *
insert(
    node_t *parent,
    int position)
// SYNOPSIS:                                                              
//   Back end function to construct tries. Append a child to an existing  
//   node at specifieid position.                                         
//   NO CHECKING IS PERFORMED to make sure that this does not overwrite   
//   an existings node child (causing a memory leak) or that 'c' is an    
//   integer less than 5. Since 'insert' is called exclusiverly by        
//   'insert_string' after a call to 'find_path', this checking is not    
//   required. If 'insert' is called in another context, this check has   
//   to be performed.                                                     
//                                                                        
// PARAMETERS:                                                            
//   parent: the parent to append the node to                             
//   position: the position of the child                                  
//                                                                        
// RETURN:                                                                
//   The appended child node in case of success, 'NULL' otherwise.        
//                                                                        
// NB: This function is not used by 'starcode()'.
{
   // Initilalize child node.
   node_t *child = new_trienode();
   if (child == NULL) {
      fprintf(stderr, "error: could not insert node\n");
      ERROR = __LINE__;
      return NULL;
   }
   // Update child path and parent pointer.
   child->path = (parent->path << 4) + position;
   // Update parent node.
   parent->child[position] = child;

   return child;

}

PYBIND11_MODULE(tree, m)
{
   m.doc() = "Python binding for tree search algorithms";
   namespace py = pybind11;
   py::class_<trie_t>(m, "Tree");
   py::class_<result_t>(m, "Result")
      .def(py::init<>())
      .def_readwrite("label", &result_t::label)
      .def_readwrite("distance", &result_t::distance);

   m.def("new_tree", &new_trie, "Construct tree");
   m.def("search", &search, "Search in the tree");
   m.def("insert", &insert_string, "Insert a string to the tree");
}
