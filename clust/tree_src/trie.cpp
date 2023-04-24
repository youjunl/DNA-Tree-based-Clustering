#include "trie.h"
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

#define PAD 5              // Position of padding nodes.
#define EOS -1             // End Of String, for 'dash()'.
using namespace std;
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
   int         tau;
   int       * query;
   int         height;
};


struct item_t
{
   int distance;
   int * query;
};

struct context_t
{
   int depth;
   long label;
   bool ins[4] = {false, false, false, false};
   bool del[4] = {false, false, false, false};
   bool sub[4] = {false, false, false, false};
   bool error = false;
};

int get_height(trie_t *);
node_t *insert(node_t *, int);
node_t *new_trienode(void);
result_t *poucet(node_t *, int, int*, char*, struct arg_t);
char * compute_ccache(char *, int *, int*, const int, const int);
context_t * context(trie_t *, int*, const int, const int);

// Globals.
int ERROR = 0;

int get_height(trie_t *trie) { return trie->height; }

// ------  SEARCH FUNCTIONS ------ //
result_t *quick_search(trie_t *trie, const char *query, const int tau, const int max_depth)
{
   int height = get_height(trie);
   int length = strlen(query);
   // Translate symbols to integers
   int translated[M];
   translated[length] = EOS;
   for (int i = 0; i < length; i++)
   {
      translated[i] = altranslate[(int)query[i]];
   }
   node_t *start_node = trie->root;
   int out[length];

   queue<item_t> q;
   item_t input = {0, translated};

   result_t *result = new result_t;

   q.push(input);
   while (true)
   {
      // If found perfect matched
      if(result->distance==0)break;
      
      // If no candidate in queue
      if(q.empty())break;

      // Get a candidate input
      item_t cur_input = q.front();
      q.pop();

      if(cur_input.distance>tau)continue;

      context_t *context_out = context(trie, cur_input.query, height, max_depth);
      if (context_out->error)
      {
         // Found unmatching during traversal
         int pos = context_out->depth;
         // Insertion fix
         for(int i=0; i<4; i++)
         {
            if(!context_out->ins[i])continue;
            int new_query[8];
            item_t new_input = {cur_input.distance+1, new_query};
            for(int j=0;j<pos-1;j++)new_query[j] = cur_input.query[j];
            for(int j=pos;j<height;j++)new_query[j-1] = cur_input.query[j];
            new_query[height] = 0; // Add a PAD to the end
            q.push(new_input);          
         }
         // Deletion fix
         for(int i=0; i<4; i++)
         {
            if(!context_out->del[i])continue;
            int new_query[8];
            item_t new_input = {cur_input.distance+1, new_query};
            for(int j=0;j<pos;j++)new_query[j] = cur_input.query[j];
            new_query[pos] = i+1;
            for(int j=pos+1;j<height;j++)new_query[j] = cur_input.query[j-1];
            q.push(new_input);
         }
         // Substitution fix
         for(int i=0; i<4; i++)
         {
            if(!context_out->sub[i])continue;
            int new_query[8];
            item_t new_input = {cur_input.distance+1, new_query};
            for(int j=0;j<pos;j++)new_query[j] = cur_input.query[j];
            new_query[pos] = i+1;
            for(int j=pos+1;j<height;j++)new_query[j] = cur_input.query[j];
            q.push(new_input);
         }
      }
      else
      {
         // Found an output
         int cur_distance = cur_input.distance;
         long cur_label = context_out->label;
         if(result->distance>cur_distance)
         {
            result->label=cur_label;
            result->distance=cur_distance;
         }
      }
   }
   return result;
}

context_t * context(trie_t *trie, int * query, const int height, const int max_depth)
{
   node_t *node = trie->root;
   context_t *out = new context_t;
   for(int i=0; i < height; i++)
   {
      int ch_num = query[i];
      out->depth = i;
      if(node->child[ch_num] == NULL)
      {
         out->error = true;
         // Pretend overflow
         int depth;
         int search_depth = min(height-i-1, max_depth);
         for(int node_num = 1; node_num<5; node_num++)
         {
            node_t *tmp_node;
            if(node->child[node_num] != NULL)
            {
               // Deletion
               depth = 0;
               tmp_node = node->child[node_num];
               while (depth<search_depth)
               {
                  if((tmp_node=tmp_node->child[query[i+depth]]) == NULL)break;
                  depth++;
               }
               if (depth == search_depth)out->del[node_num - 1] = true;

               // Insertion
               depth = 0;
               tmp_node = node;
               while (depth<search_depth-1)
               {
                  if((tmp_node=tmp_node->child[query[i+depth+1]]) == NULL)break;
                  depth++;
               }
               if (depth == search_depth-1)out->ins[node_num - 1] = true;

               // Substitution
               depth = 0;
               tmp_node = node->child[node_num];
               while (depth<search_depth)
               {
                  if((tmp_node=tmp_node->child[query[i+depth+1]]) == NULL)break;
                  depth++;
               }
               if (depth == search_depth)out->sub[node_num - 1] = true;               
            }
         }
         return out;
      }
      // If node exist, continue
      else node = node->child[ch_num];
   }
   if(node->isEnd)out->label=node->label;
   else if(!out->error)
   {
      // Error not found, but terminated early because of insertions
      // Add a deletion flag to the end
      out->error = true;
      for(int node_num = 1; node_num<5; node_num++)
      {
         if(node->child[node_num] != NULL)out->del[node_num-1]=true;
      } 
   }
   
   return out;
}

result_t *search(trie_t *trie, const char *query, const int tau)
{
   int height = get_height(trie);
   int length = strlen(query);

   // Translate symbols to integers
   int translated[M];
   translated[length] = EOS;
   for (int i = 0; i < length; i++)
   {
      translated[i] = altranslate[(int)query[i]];
   }

   // Set the search options.
   struct arg_t arg;
   arg.tau = tau;
   arg.query = translated;
   arg.height = height;

   node_t *start_node = trie->root;
   int out[length];
   
   // Initialize cache
   char cache[2 * tau + 1];
   for (int i = 0; i < 2 * tau + 1; i++)
      cache[i] = MAX_TAU;
   cache[tau] = 0;
   return poucet(start_node, 0, out, cache, arg);
}

char *compute_ccache(char *pcache, int *inp, int *out, const int depth, const int tau)
{
   const int n = 2 * tau + 1;
   char *ccache = (char *)malloc(n);
   int delta;
   for (int i = 0; i < 2 * tau + 1; i++)
      ccache[i] = pcache[i];
   if (depth + 1 > tau)
   {
      ccache[0] = min(pcache[0] + (out[depth] != inp[0]), pcache[1] + 1);
      ccache[n - 1] = min(pcache[n - 1] + (out[depth - tau] != inp[tau - 1]), pcache[n - 2] + 1);
      for (int i = 1; i < tau; i++)
      {
         // Horizental arm
         int hshift = min(ccache[i - 1] + 1, pcache[i + 1] + 1);
         delta = out[depth] != inp[i - 1];
         ccache[i] = min(pcache[i] + delta, hshift);
         // Vertical arm
         int j = n - i - 1;
         int vshift = min(ccache[j + 1] + 1, pcache[j - 1] + 1);
         delta = out[depth - tau + i] != inp[tau - 1];
         ccache[j] = min(pcache[j] + delta, vshift);
      }
   }
   else
   {
      int offset = tau - depth - 1;
      for (int i = offset + 1; i < tau; i++)
      {
         // Horizental arm
         int hshift = min(ccache[i - 1] + 1, pcache[i + 1] + 1);
         delta = out[depth] != inp[i - offset - 1];
         ccache[i] = min(pcache[i] + delta, hshift);
         // Vertical arm
         int j = n - i - 1;
         int vshift = min(ccache[j + 1] + 1, pcache[j - 1] + 1);
         delta = out[i - offset - 1] != inp[depth];
         ccache[j] = min(pcache[j] + delta, vshift);
      }
   }
   // Center cell
   int shift = min(ccache[tau - 1] + 1, ccache[tau + 1] + 1);
   if (depth + 1 > tau)
      ccache[tau] = min(pcache[tau] + (out[depth] != inp[tau - 1]), shift);
   else
      ccache[tau] = min(pcache[tau] + (out[depth] != inp[depth]), shift);
   return ccache;
}

result_t *poucet(node_t *node, const int depth, int *out, char *pcache, struct arg_t arg)
{
   // This makes it easier to distinguish the part that goes upward,
   // with positive index and requiring the path, from the part that
   // goes horizontally, with negative index and requiring previous
   // characters of the query.
   result_t *result = new result_t;
   const int tau = arg.tau;

   // Exceed threshold
   if (pcache[tau] > tau || depth > arg.height)
   {
      result->distance = pcache[tau];
      return result;
   }
 
   if (node->isEnd)
   {
      result->label = node->label;
      result->distance = pcache[tau];
      return result;
   }

   // Update cache edge
   if (depth < tau)
   {
      pcache[tau - depth - 1] = depth + 1;
      pcache[tau + depth + 1] = depth + 1;
   }

   node_t *child;
   int *inp = arg.query + max(0, depth + 1 - tau);
   int next_ind = arg.query[depth];

   // If node exist for the query[depth]
   if ((child = node->child[next_ind]) != NULL)
   {
      out[depth] = next_ind;
      char * ccache = compute_ccache(pcache, inp, out, depth, tau);
      result_t *tmp_result = poucet(child, depth + 1, out, ccache, arg);
      if (tmp_result->label != -1)
         return tmp_result;
   }

   // If traversal fail or node not exists, try other nodes
   for (int i = 1; i < 5; i++)
   {
      // Skip if current node has no child at this position.
      if (i == next_ind || (child = node->child[i]) == NULL)
         continue;
      out[depth] = i;
      char * ccache = compute_ccache(pcache, inp, out, depth, tau);
      result_t *tmp_result = poucet(child, depth + 1, out, ccache, arg);
      if (tmp_result->label != -1)
         return tmp_result;
   }
   // There is no more node for traversal
   return result;
}

// ------  TRIE CONSTRUCTION AND DESTRUCTION  ------ //

trie_t * new_trie(unsigned int height)
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

   // Set the values of the meta information.
   trie->root = root;
   trie->height = height;

   return trie;

}

node_t * new_trienode(void)
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

void insert_string(trie_t *trie, const char *string, long label)
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
   for (i = 0; i < nchar; i++)
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
   for (; i < nchar; i++)
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

node_t * insert(node_t *parent, int position)
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
   if (child == NULL)
   {
      fprintf(stderr, "error: could not insert node\n");
      ERROR = __LINE__;
      return NULL;
   }
   // Update parent node.
   child->ch_num = position;
   parent->child[position] = child;
   return child;
}

void delete_string(trie_t *trie, const char*string)                                                            
{
   int i;
   int nchar = strlen(string);   
   // Find existing path.
   node_t *node = trie->root;
   for (i = 0; i < nchar; i++)
   {
      node_t *child;
      int c = translate[(int)string[i]];
      if ((child = (node_t *)node->child[c]) == NULL)
      {
         break;
      }
      node = child;
   }
   if(node->isEnd)
   {
      node->isEnd = false;
   }
   return;
}