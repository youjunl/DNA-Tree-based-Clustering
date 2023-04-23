#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#ifndef _TRIE_HEADER
#define _TRIE_HEADER

#define DESTROY_NODES_YES 1
#define DESTROY_NODES_NO 0
#define MAX_TAU 255
static const char BASES[8] = "ACGTN";

struct info_t;
struct node_t;
struct trie_t;
struct result_t;

typedef struct info_t info_t;
typedef struct node_t node_t;
typedef struct trie_t trie_t;
typedef struct result_t result_t;
// Global constants.
#define TAU 8           // Max Levenshtein distance.
#define M 1024          // MAXBRCDLEN + 1, for short.
#define MAXBRCDLEN 1023 // Maximum barcode length.

void insert_string(trie_t *, const char *, unsigned int);
trie_t *new_trie(unsigned int);
result_t *search(trie_t *, const char *, int);
struct trie_t
{
   node_t *root;
   unsigned int height;
};

struct node_t
{
   node_t *child[6] = {nullptr}; // Array of 6 children pointers.
   char cache[2 * TAU + 1];      // Dynamic programming space.
   bool isEnd = false;           // Indentifier of leaf
   unsigned int label = 0;       // Lable on the leaf
   char ch_num;
};

struct result_t
{
   int label = -1;      // return label
   int distance = MAX_TAU; // edit distance
};

#endif
