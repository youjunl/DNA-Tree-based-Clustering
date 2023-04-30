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
    inp.push_back("ATTGCATAAT");
    inp.push_back("ATTGCATAAT");
    inp.push_back("ATTCGATTTC"); // New sequence
    inp.push_back("ATTGGATTTC"); // New sequence with error

    // Create tree
    trie_t *tree = new_trie(10);
    const int tau = 3;
    int clust_ind = 1;
    for(size_t i=0; i<inp.size(); i++)
    {
        const char *seq = inp[i].c_str();
        result_t *tmp = quick_search(tree, seq, tau, 2);
        result_t *tmp_slow = poucet_search(tree, seq, tau);
        if (tmp->label == 0)
        {
            cout << "Poucet" << to_string(tmp_slow->label) << " " << to_string(tmp_slow->distance) << endl;
            insert_string(tree, seq, clust_ind);
            clust_ind++;
        }
        cout << seq << " " << to_string(tmp->label) << " " << to_string(tmp->distance) << endl;
    }
    return 0;
}