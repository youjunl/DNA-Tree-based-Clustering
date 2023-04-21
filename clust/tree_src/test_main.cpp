#include<iostream>
#include<stdlib.h>
#include<vector>
#include<string>
#include "trie.h"
using namespace std;

int main()
{
    cout << "Demo" << endl;
    vector<string> inp;
    inp.push_back("ATTGCCCA");
    inp.push_back("ATTGCCCA");
    inp.push_back("CAGCAAAA");
    inp.push_back("CAGCAAAA");
    inp.push_back("ATTGACCA");
    inp.push_back("ATCGCCCA");
    
    // Create tree
    trie_t *tree = new_trie(8);
    const int tau = 2;
    int clust_ind = 1;
    for(size_t i=0; i<inp.size(); i++)
    {
        const char *seq = inp[i].c_str();
        result_t *tmp = search(tree, seq, tau);
        if (tmp->label < 0)
        {
            insert_string(tree, seq, clust_ind);
            clust_ind++;
        }
        cout << seq << " " << to_string(tmp->label) << " " << to_string(tmp->distance) << endl;
    }
    return 0;
}