#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <queue>
#include <vector>
#include <map>
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
#define M 300           // Max length of input size
#define TAU 8           // Max Levenshtein distance.
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
   node_t *child[4] = {nullptr}; // Array of 6 children pointers.
   bool isEnd = false;           // Indentifier of leaf
   unsigned long label;       // Lable on the leaf
   char ch_num;
};

struct result_t
{
   unsigned long label = 0;      // return label
   int distance = MAX_TAU; // edit distance
};

#endif
