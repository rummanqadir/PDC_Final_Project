#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <mpi.h>

using namespace std;

static constexpr int n = 11;
struct Record {
    array<unsigned char, n>             node;
    array<array<unsigned char,n>, n-1> parents;
};

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
    const vector<unsigned char>& identity , vector<int> inv, int r) {
int n = v.size();
// vector<int> inv = inverse(v);
// int r= pos(v);

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
vector<unsigned char> get_kth_permutation(long long k, int n) {
    vector<unsigned char> numbers(n);
    iota(numbers.begin(), numbers.end(), 1);
    vector<unsigned char> perm;
    long long fact = 1;
    for (int i = 1; i < n; ++i) fact *= i; // (n-1)!
    //cout << fact << endl;
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

int main(int argc,char**argv){
    MPI_Init(&argc,&argv);
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    MPI_Get_processor_name(hostname, &namelen);
    // total = n!
    long long total = factorial(n);
    
    long long chunk = total/size;
    long long lo = rank*chunk;
    long long hi = (rank+1==size ? total : lo+chunk);
    cout << "rank: " << rank <<", hostname: " << hostname << ", kStart: " << lo<< " kEnd: " << hi << endl;
    // vector<vector<unsigned char>> perms(hi-lo);
    vector<Record> table(hi-lo);
    // identity
    vector<unsigned char> id(n);
    iota(id.begin(),id.end(),1);

    //cout << "ID : " <<  perm_to_string(id)<<endl;
    double t0 = MPI_Wtime();
    // generate and store parents
    // for(long long k=lo;k<hi;++k){
    //     auto v = get_kth_permutation(k,n);
    //     // perms[k-lo] = v;
    //     auto &R = table[k-lo];
    //     R.node = v;
    //     //cout << perm_to_string(R.node)<< " ";
    //     for(int t=1;t<n;++t){
    //         R.parents[t-1] = parent1(v,t,id);
    //         //cout<< perm_to_string(R.parents[t-1]) << " ";
    //     }
    //     //cout << endl;
    // }
    for (long long k = lo; k < hi; ++k) {
        auto perm = get_kth_permutation(k, n);      // returns vector<uchar> of length n
        Record &R = table[k - lo];
        // copy into fixed array:
        vector<int> inv = inverse(perm);
        int r= pos(perm);
        for (int i = 0; i < n; ++i) R.node[i] = perm[i];
        for (int t = 1; t < n; ++t) {
          auto p = parent1(perm, t, id,inv,r);
          for (int i = 0; i < n; ++i) R.parents[t-1][i] = p[i];
        }
      }

    double t1 = MPI_Wtime();
    cout << "Rank " << rank << " did " << (hi - lo) << " records in " << (t1 - t0) << " s\n";
    // if(rank==0){
    //     cout<<"Rank 0 built ["<<lo<<","<<hi<<") in "<<(t1-t0)<<"s\n";
    //     // interactive lookup
    //     long long idx;
    //     cout<<"Enter lex index (0.."<<total-1<<"): ";
    //     cin>>idx;
    //     if(idx<lo||idx>=hi){
    //         cout<<"Index on another rank.\n";
    //     } else {
    //         auto &perm = perms[idx-lo];
    //         cout<<"Perm: ";
    //         for(auto c:perm) cout<<int(c);
    //         cout<<"\nParents:\n";
    //         auto &R = table[idx-lo];
    //         for(int t=1;t<n;++t){
    //             cout<<" t="<<t<<": ";
    //             for(auto c:R.parents[t-1]) cout<<int(c);
    //             cout<<"\n";
    //         }
    //     }
    // }
    MPI_Finalize();
    return 0;
}
