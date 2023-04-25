#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <queue>
#include <vector>

#ifndef _TRIE_HEADER
#define _TRIE_HEADER

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
#define MAX_GREEDY_DEPTH 2
void insert_string(trie_t *, const char *, long);
void delete_string(trie_t *, const char *);
trie_t *new_trie(unsigned int);
result_t *search(trie_t *, const char *, const int);
result_t *quick_search(trie_t *, const char *, const int, const int);
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
   long label;       // Lable on the leaf
   char ch_num;
};

struct result_t
{
   long label = -1;      // return label
   int distance = MAX_TAU; // edit distance
};

#endif
