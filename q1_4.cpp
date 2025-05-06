#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <numeric>
#include <mpi.h>
#include "functions.h"
using namespace std;
const int n = 11;
const int BATCH_SIZE = 10000;
class permNode{
    public:
        vector<unsigned char> perm;
        vector<vector<unsigned char>> parents;
        vector<vector<unsigned char>> parentindex;
    permNode (){
        parentindex.resize(n-1);
        parents.resize(n-1);
    }
};
// Swap the symbol x with the symbol at its adjacent position
vector<unsigned char> Swap(const vector<unsigned char>& v, int x, const vector<int>& inv) {
    int n = v.size();
    int i = inv[x];            // correct 0-based lookup
    if (i >= 0 && i+1 < n) {
        vector<unsigned char> p = v;
        std::swap(p[i], p[i+1]);
        return p;
    }
    return v;
}

vector<int> inverse(const vector<unsigned char>& v){
    int n = v.size();
    vector<int> inv(n);
    for (int i = 0; i < n; ++i) {
        inv[v[i] - 1] = i;
    }
    return inv;
}
int pos(const vector<unsigned char>& v){
    int r = 0;
    int n=v.size();
    for (int i = n - 1; i >= 0; --i) {
        if (v[i] != i + 1) {
            r = i + 1; // 1-based rightmost position where v_i != i
            return r;
        }
    }
    return r;
}

// FindPosition function as per the paper's algorithm
vector<unsigned char> FindPosition(const vector<unsigned char>& v, int t, 
                                   const vector<unsigned char>& identity, int r, const vector<int>& inv) {
    int n = v.size();
    if (t == 2) {
        vector<unsigned char> swapped = Swap(v, t-1,inv);
        if (swapped == identity) {
            return Swap(v, t - 2,inv);
        }
    }
    if (v[n - 2] == t || v[n - 2] == n - 1) { // v_{n-1} in 1-based indexing
        if (r > 0) {
            int k = r - 1; // Convert to 0-based
            return Swap(v, k,inv);
        }
        return v;
    }
    return Swap(v, t-1,inv);
}

// Parent1 function: computes parent of v in tree t
vector<unsigned char> parent1(const vector<unsigned char>& v, int t, 
                              const vector<unsigned char>& identity) {
    int n = v.size();
    vector<int> inv = inverse(v);
    int r= pos(v);
    
    // CASE A: v_n == n
    if (v[n - 1] == n) {
        if (t != n - 1) {
            return FindPosition(v, t, identity, r, inv);
        } else {
            return Swap(v, v[n - 2]-1,inv); // Swap v_{n-1}
        }
    }
    // CASE B: v_n == n-1 and v_{n-1} == n and Swap(v,n) != identity
    else if (v[n - 1] == n - 1 && v[n - 2] == n && Swap(v, n-1,inv) != identity) {
        if (t == 1) {
            return Swap(v, n-1,inv);
        } else {
            return Swap(v, t - 2,inv);
        }
    }
    // Otherwise
    else {
        if (v[n - 1] == t) {
            return Swap(v, n-1,inv);
        } else {
            return Swap(v, t-1,inv);
        }
    }
}
// Convert vector<unsigned char> to string
string perm_to_string(const vector<unsigned char>& p) {
    string s;
    for (unsigned char val : p)
        s += to_string((int)val);
    return s;
}

// Compute the k-th permutation (0-based index) for n elements
vector<unsigned char> get_kth_permutation(long long k, int n,long long fact) {
    vector<unsigned char> numbers(n);
    iota(numbers.begin(), numbers.end(), 1);
    vector<unsigned char> perm;
    // long long fact = 1;
    // for (int i = 1; i < n; ++i) fact *= i; // (n-1)!
    
    k = k % (fact * n);
    for (int i = n - 1; i >= 0; --i) {
        if (i == 0) {
            perm.push_back(numbers[0]);
            break;
        }
        int index = k / fact;
        k = k % fact;
        fact /= i;
        perm.push_back(numbers[index]);
        numbers.erase(numbers.begin() + index);
    }
    return perm;
}

// Compute factorial
long long factorial(int n) {
    long long res = 1;
    for (int i = 1; i <= n; ++i) res *= i;
    return res;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    MPI_Get_processor_name(hostname, &namelen);
    
    const long long total_perms = factorial(n); // 24
    const int perms_per_process = total_perms / size;
    long long k_start = rank * perms_per_process;
    long long k_end = (rank == size - 1) ? total_perms : k_start + perms_per_process;
    const long long fact = factorial(n-1);
    // Ensure 12 perms each for 2 processes
    
     
    k_start = (total_perms/float(6))*rank;
    k_end = (total_perms/float(6))*(rank+1);
    cout << "rank: " << rank << " kStart: " << k_start<< " kEnd: " << k_end << endl;
    
    permNode node;
    double start_time = MPI_Wtime();


    // Open output file
    ofstream out("../ist_parents_n" + to_string(n) + "_rank" + to_string(rank) + ".txt");
    if (!out) {
        cerr << "Error opening output file for rank " << rank << ".\n";
        MPI_Finalize();
        return 1;
    }

    ostringstream buffer;
    
    long long count = 0;
    vector<unsigned char> identity(n);
    for (int i = 0; i < n; ++i) identity[i] = i + 1;
    // Process assigned permutations
    for (long long k = k_start; k < k_end; ++k) {
        vector<unsigned char> perm = get_kth_permutation(k, n, fact);
        //auto index = permutation_to_index_naive(perm);
        // for (unsigned char p:perm)
        //     cout << int(p);
        // cout << endl;
        // Write permutation
        for (unsigned char val : perm)
            buffer << (int)val;
        buffer << "  ";

        node.perm = perm;
        // Compute parents for t = 1 to n-1
        for (int t = 1; t <= n - 1; ++t) {
            node.parents[t-1] = parent1(perm, t,identity);
            node.parentindex[t-1] = ll_to_char_digits(permutation_to_index_naive(node.parents[t-1]));
        }
        for (vector<unsigned char> n : node.parentindex){
            for (unsigned char val : n)
                buffer << val;
            buffer << "  ";
        }
        buffer << "\n";
        ++count;

        if (count % BATCH_SIZE == 0) {
            out << buffer.str();
            buffer.str("");
            buffer.clear();
        }
    }

    // Write remaining buffer
    out << buffer.str();
    out.close();
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    cout << "⏱️ Rank " << rank <<", hostname: " << hostname << " finished in " << elapsed_time << " seconds." << endl;

    if (rank == 0) {
        cout << "✅ Done. Total permutations: " << total_perms << "\n";
        cout << "📁 Results saved to: ist_parents_n" << n << "_rankX.txt\n";
    }

    MPI_Finalize();
    return 0;
}