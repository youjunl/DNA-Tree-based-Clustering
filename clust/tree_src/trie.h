#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <queue>
#include <limits>
#include <vector>
#include "darts.h"

#ifndef _TRIE_HEADER
#define _TRIE_HEADER

#define MAX_TAU 127

struct info_t;
struct node_t;
struct trie_t;
struct result_t;

typedef struct info_t info_t;
typedef struct node_t node_t;
typedef struct trie_t trie_t;
typedef struct result_t result_t;
// Global constants.
#define M 200           // Max length of input size
#define TAU 5           // Max Levenshtein distance.
void insert_string(trie_t *, const char *, const long);
void delete_string(trie_t *, const char *);
trie_t *new_trie(const unsigned int);
result_t *poucet_search(trie_t *, const char *, const int);
result_t *quick_search(trie_t *, const char *, const int, const int);

using namespace std;

struct node_t
{
   // node_t *child[4] = {nullptr}; // Array of 6 children pointers.
   vector<node_t*> child = vector<node_t*>(4);
   bool isEnd = false;           // Indentifier of leaf
   unsigned long label;       // Lable on the leaf
};

struct trie_t
{
   node_t *root;
   unsigned int height;
};

struct result_t
{
   unsigned long label = 0;      // return label
   int distance = MAX_TAU; // edit distance
};

#endif
