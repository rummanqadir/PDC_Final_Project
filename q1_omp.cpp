#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cstdint>
#include <numeric>
#include <mpi.h>
#include <map>
#include <string>
#include <sstream>
#include <unordered_map>
#include <omp.h>


using namespace std;

const int n = 11;
const int BATCH_SIZE = 10000;
struct Record {
    uint32_t idx;
    // vector<unsigned char> perm;
    array<uint32_t, n - 1> parents;
};




//— a simple hash for fixed-length vector<unsigned char> permutations —
struct PermHash {
  size_t operator()(vector<unsigned char> const& v) const noexcept {
    // e.g. a rolling‐XOR/FNV‐like mix
    size_t h = 146527;
    for ( auto c : v ) h = (h ^ c) * 0x9e3779b97f4a7c15ULL;
    return h;
  }
};

unordered_map<vector<unsigned char>,Record,PermHash> permToRecord;



vector<unsigned char> Swap(const vector<unsigned char>& v, int x, const vector<int>& inv) {
    int n = v.size();
    int i = inv[x];
    if (i >= 0 && i + 1 < n) {
        auto p = v;
        std::swap(p[i], p[i + 1]);
        return p;
    }
    return v;
}

vector<int> inverse(const vector<unsigned char>& v) {
    int n = v.size();
    vector<int> inv(n);
    for (int i = 0; i < n; ++i)
        inv[v[i] - 1] = i;
    return inv;
}

int pos(const vector<unsigned char>& v) {
    int n = v.size();
    for (int i = n - 1; i >= 0; --i) {
        if (v[i] != i + 1) return i + 1;
    }
    return 0;
}

vector<unsigned char> FindPosition(const vector<unsigned char>& v, int t,
                                  const vector<unsigned char>& identity, int r, const vector<int>& inv) {
    int n = v.size();
    if (t == 2) {
        auto swapped = Swap(v, t - 1, inv);
        if (swapped == identity)
            return Swap(v, t - 2, inv);
    }
    if (v[n - 2] == t || v[n - 2] == n - 1) {
        if (r > 0) return Swap(v, r - 1, inv);
        return v;
    }
    return Swap(v, t - 1, inv);
}

vector<unsigned char> parent1(const vector<unsigned char>& v, int t,
                              const vector<unsigned char>& identity) {
    int n = v.size();
    auto inv = inverse(v);
    int r = pos(v);

    if (v[n - 1] == n) {
        if (t != n - 1)
            return FindPosition(v, t, identity, r, inv);
        else
            return Swap(v, v[n - 2] - 1, inv);
    }
    else if (v[n - 1] == n - 1 && v[n - 2] == n && Swap(v, n - 1, inv) != identity) {
        if (t == 1)
            return Swap(v, n - 1, inv);
        else
            return Swap(v, t - 2, inv);
    }
    else {
        if (v[n - 1] == t)
            return Swap(v, n - 1, inv);
        else
            return Swap(v, t - 1, inv);
    }
}

long long factorial(int n) {
    long long f = 1;
    for (int i = 1; i <= n; ++i) f *= i;
    return f;
}

vector<unsigned char> get_kth_permutation(long long k, int n) {
    vector<unsigned char> nums(n);
    iota(nums.begin(), nums.end(), 1);
    vector<unsigned char> perm;
    perm.reserve(n);

    long long fact = 1;
    for (int i = 1; i < n; ++i) fact *= i;

    k = k % (fact * n);
    for (int i = n - 1; i >= 0; --i) {
        if (i == 0) {
            perm.push_back(nums[0]);
            break;
        }
        int idx = k / fact;
        k %= fact;
        fact /= i;
        perm.push_back(nums[idx]);
        nums.erase(nums.begin() + idx);
    }
    return perm;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    MPI_Get_processor_name(hostname, &namelen);
    uint32_t total = factorial(n);
    uint32_t perRank = total / size;
    uint32_t start = perRank * rank;
    uint32_t end = (rank == size - 1 ? total : start + perRank);

    //for only my pc run only on rank0 and 1
    // uint32_t total = factorial(n);
    // uint32_t perRank = total / 6;
    // uint32_t start = perRank * rank;
    // uint32_t end = (rank == 6 - 1 ? total : start + perRank);
    cout << "rank: " << rank  << " HOSTNAME: " << hostname << " kStart: " << start<< " kEnd: " << end << endl;
    vector<unsigned char> identity(n);
    iota(identity.begin(), identity.end(), 1);
    vector<vector<unsigned char>> permutations(end-start);
    double t0 = MPI_Wtime();
    
    // Open output file
    // ofstream out("../ist_parents_n" + to_string(n) + "_rank" + to_string(rank) + ".txt");
    // if (!out) {
    //     cerr << "Error opening output file for rank " << rank << ".\n";
    //     MPI_Finalize();
    //     return 1;
    // }

    // ostringstream buffer;
    
    // long long count = 0;
    // before the parallel region
    omp_lock_t map_lock;
    omp_init_lock(&map_lock);

    int nprocs = omp_get_num_procs();
    // if(hostname == "client1")
    //     nprocs/=2;   
    double kLoopTimeStart = MPI_Wtime();
    #pragma omp parallel for num_threads(nprocs)
    for (uint32_t k = start; k < end; ++k) {
        // 1) build the permutation and record
        auto perm = get_kth_permutation(k, n);
        permutations[k - start] = perm;
        Record rec; rec.idx = k;

        // 2) protect the map insertion
        omp_set_lock(&map_lock);
        permToRecord.emplace(std::move(perm), rec);
        omp_unset_lock(&map_lock);
    }
    double kLoopTimeEnd = MPI_Wtime();

// after you’re done
    omp_destroy_lock(&map_lock);



        // … after you’ve built permToRecord …

        double permtoRecordTime = MPI_Wtime();

        // 1) Build a flat vector of (key*, record*) pairs
        vector<pair<const vector<unsigned char>*, Record*>> entries;
        entries.reserve(permToRecord.size());
        for (auto & kv : permToRecord) {
            entries.emplace_back(&kv.first, &kv.second);
        }
    
        // 2) Parallel parent‐finding
        #pragma omp parallel for num_threads(nprocs)
        for (size_t i = 0; i < entries.size(); ++i) {
            auto key   = *entries[i].first;    // make a local copy of the permutation
            auto rec   = entries[i].second;    // pointer to the Record to fill
            for (int t = 1; t < n; ++t) {
                auto p        = parent1(key, t, identity);
                auto itParent = permToRecord.find(p);
                rec->parents[t-1] = (itParent != permToRecord.end()
                                      ? itParent->second.idx
                                      : 0u);
            }
        }
    
        double permToRecordEndTime = MPI_Wtime();
        cout << "Parent loop time: "
             << (permToRecordEndTime - permtoRecordTime)
             << " s\n";
    

    // out << buffer.str();
    // out.close();
    double t1 = MPI_Wtime();
    cout << "Rank " << rank << " did " << (end - start) << " records in " << (t1 - t0) << " s with : "<< nprocs<< " processors k loop time : " << kLoopTimeEnd - kLoopTimeStart << " s " << " perm time : " << permToRecordEndTime - permtoRecordTime  << " s\n";
//     if(rank==0){
//     for (uint32_t i = start; i<end; i++){
//         auto perm = permutations[i-start];
//         auto it = permToRecord.find(perm);
//         for (auto p : perm)
//             cout << int(p);
//         cout << " ";
//         for(int t=0; t<n-1;t++){
//             uint32_t parent_idx = it->second.parents[t];
//             perm = permutations[parent_idx];
//             for (auto p : perm)
//                 cout << int(p);
//             cout << " ";
//         }
//         cout << endl;
//     }
// }
// int a;
// if(rank == 0){
//     cin >> a;
// }
    MPI_Finalize();
    return 0;
}
