// #include <iostream>
// #include <fstream>
// #include <array>
// #include <vector>
// #include <cstdint>
// #include <numeric>
// #include <mpi.h>

// using namespace std;

// const int n = 10;

// struct Record {
//     uint32_t idx;
//     array<vector<unsigned char>, n - 1> parents;  // Store permutations directly
// };

// vector<unsigned char> Swap(const vector<unsigned char>& v, int x, const vector<int>& inv) {
//     int n = v.size();
//     int i = inv[x];
//     if (i >= 0 && i + 1 < n) {
//         auto p = v;
//         std::swap(p[i], p[i + 1]);
//         return p;
//     }
//     return v;
// }

// vector<int> inverse(const vector<unsigned char>& v) {
//     int n = v.size();
//     vector<int> inv(n);
//     for (int i = 0; i < n; ++i)
//         inv[v[i] - 1] = i;
//     return inv;
// }

// int pos(const vector<unsigned char>& v) {
//     int n = v.size();
//     for (int i = n - 1; i >= 0; --i) {
//         if (v[i] != i + 1) return i + 1;
//     }
//     return 0;
// }

// vector<unsigned char> FindPosition(const vector<unsigned char>& v, int t,
//                                   const vector<unsigned char>& identity, int r, const vector<int>& inv) {
//     int n = v.size();
//     if (t == 2) {
//         auto swapped = Swap(v, t - 1, inv);
//         if (swapped == identity)
//             return Swap(v, t - 2, inv);
//     }
//     if (v[n - 2] == t || v[n - 2] == n - 1) {
//         if (r > 0) return Swap(v, r - 1, inv);
//         return v;
//     }
//     return Swap(v, t - 1, inv);
// }

// vector<unsigned char> parent1(const vector<unsigned char>& v, int t,
//                               const vector<unsigned char>& identity) {
//     int n = v.size();
//     auto inv = inverse(v);
//     int r = pos(v);

//     if (v[n - 1] == n) {
//         if (t != n - 1)
//             return FindPosition(v, t, identity, r, inv);
//         else
//             return Swap(v, v[n - 2] - 1, inv);
//     }
//     else if (v[n - 1] == n - 1 && v[n - 2] == n && Swap(v, n - 1, inv) != identity) {
//         if (t == 1)
//             return Swap(v, n - 1, inv);
//         else
//             return Swap(v, t - 2, inv);
//     }
//     else {
//         if (v[n - 1] == t)
//             return Swap(v, n - 1, inv);
//         else
//             return Swap(v, t - 1, inv);
//     }
// }

// long long factorial(int n) {
//     long long f = 1;
//     for (int i = 1; i <= n; ++i) f *= i;
//     return f;
// }

// vector<unsigned char> get_kth_permutation(long long k, int n) {
//     vector<unsigned char> nums(n);
//     iota(nums.begin(), nums.end(), 1);
//     vector<unsigned char> perm;
//     perm.reserve(n);

//     long long fact = 1;
//     for (int i = 1; i < n; ++i) fact *= i;

//     k = k % (fact * n);
//     for (int i = n - 1; i >= 0; --i) {
//         if (i == 0) {
//             perm.push_back(nums[0]);
//             break;
//         }
//         int idx = k / fact;
//         k %= fact;
//         fact /= i;
//         perm.push_back(nums[idx]);
//         nums.erase(nums.begin() + idx);
//     }
//     return perm;
// }

// int main(int argc, char** argv) {
//     MPI_Init(&argc, &argv);
//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     // uint32_t total = factorial(n);
//     // uint32_t perRank = total / size;
//     // uint32_t start = perRank * rank;
//     // uint32_t end = (rank == size - 1 ? total : start + perRank);

//     //for only my pc run only on rank0 and 1
//     uint32_t total = factorial(n);
//     uint32_t perRank = total / 6;
//     uint32_t start = perRank * rank;
//     uint32_t end = (rank == 6 - 1 ? total : start + perRank);
//     cout << "rank: " << rank << " kStart: " << start<< " kEnd: " << end << endl;
//     vector<Record> table(end - start);

//     vector<unsigned char> identity(n);
//     iota(identity.begin(), identity.end(), 1);

//     double t0 = MPI_Wtime();
//     for (uint32_t k = start; k < end; ++k) {
//         auto perm = get_kth_permutation(k, n);
//         auto &rec = table[k - start];
//         rec.idx = k;
//         for (int t = 1; t < n; ++t) {
//             rec.parents[t - 1] = parent1(perm, t, identity);  // Store permutation directly
//         }
//     }
//     double t1 = MPI_Wtime();

//     cout << "Rank " << rank << " did " << (end - start) << " records in " << (t1 - t0) << " s\n";

//     MPI_Finalize();
//     return 0;
// }


#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <array>
#include <numeric>
#include <mpi.h>

using namespace std;

const int n = 10;

struct Record {
    uint32_t idx;
    array<uint32_t, n - 1> parents;
};

vector<unsigned char> Swap(const vector<unsigned char>& v, int x, const vector<int>& inv, vector<unsigned char>& result) {
    int n = v.size();
    int i = inv[x];
    if (i >= 0 && i + 1 < n) {
        result = v;
        std::swap(result[i], result[i + 1]);
        return result;
    }
    result = v;
    return result;
}

vector<int> inverse(const vector<unsigned char>& v, vector<int>& inv) {
    int n = v.size();
    inv.clear();
    inv.resize(n);
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
                                  const vector<unsigned char>& identity, int r, const vector<int>& inv,
                                  vector<unsigned char>& result) {
    int n = v.size();
    if (t == 2) {
        Swap(v, t - 1, inv, result);
        if (result == identity)
            return Swap(v, t - 2, inv, result);
    }
    if (v[n - 2] == t || v[n - 2] == n - 1) {
        if (r > 0) return Swap(v, r - 1, inv, result);
        result = v;
        return result;
    }
    return Swap(v, t - 1, inv, result);
}

vector<unsigned char> parent1(const vector<unsigned char>& v, int t,
                              const vector<unsigned char>& identity, vector<unsigned char>& result,
                              vector<int>& inv) {
    int n = v.size();
    inverse(v, inv);
    int r = pos(v);

    if (v[n - 1] == n) {
        if (t != n - 1)
            return FindPosition(v, t, identity, r, inv, result);
        else
            return Swap(v, v[n - 2] - 1, inv, result);
    }
    else if (v[n - 1] == n - 1 && v[n - 2] == n && Swap(v, n - 1, inv, result) != identity) {
        if (t == 1)
            return Swap(v, n - 1, inv, result);
        else
            return Swap(v, t - 2, inv, result);
    }
    else {
        if (v[n - 1] == t)
            return Swap(v, n - 1, inv, result);
        else
            return Swap(v, t - 1, inv, result);
    }
}

long long factorial(int n) {
    long long f = 1;
    for (int i = 1; i <= n; ++i) f *= i;
    return f;
}

vector<unsigned char> get_kth_permutation(long long k, int n, vector<unsigned char>& perm,
                                         vector<unsigned char>& nums) {
    nums.clear();
    nums.resize(n);
    iota(nums.begin(), nums.end(), 1);
    perm.clear();
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

uint32_t perm_to_index(const vector<unsigned char>& perm) {
    int n = perm.size();
    vector<int> fenwick(n + 1, 0);
    auto update = [&](int i) { for (; i <= n; i += i & -i) fenwick[i]++; };
    auto query = [&](int i) { int sum = 0; for (; i > 0; i -= i & -i) sum += fenwick[i]; return sum; };
    uint32_t idx = 0;
    vector<uint32_t> fact(n, 1);
    for (int i = 1; i < n; ++i) fact[i] = fact[i - 1] * i;
    for (int i = 0; i < n; ++i) {
        int v = perm[i];
        idx += (v - 1 - query(v)) * fact[n - 1 - i];
        update(v);
    }
    return idx;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint32_t total = factorial(n);
    uint32_t perRank = total / size;
    uint32_t start = perRank * rank;
    uint32_t end = (rank == size - 1 ? total : start + perRank);

    vector<Record> table;
    table.reserve(perRank);

    vector<unsigned char> identity(n);
    iota(identity.begin(), identity.end(), 1);
    vector<unsigned char> perm, result, nums;
    vector<int> inv;
    perm.reserve(n);
    result.reserve(n);
    nums.reserve(n);
    inv.reserve(n);

    double t0 = MPI_Wtime();
    for (uint32_t k = start; k < end; ++k) {
        get_kth_permutation(k, n, perm, nums);
        Record rec;
        rec.idx = k;
        for (int t = 1; t < n; ++t) {
            parent1(perm, t, identity, result, inv);
            rec.parents[t - 1] = perm_to_index(result);
        }
        table.push_back(rec);
        // if (table.size() >= 1000000) { // Write to disk every 1M records
        //     ofstream bin("../parents_n" + to_string(n) + "_rank" + to_string(rank) + ".bin", ios::binary | ios::app);
        //     bin.write((char*)table.data(), table.size() * sizeof(Record));
        //     bin.close();
        //     table.clear();
        //     table.reserve(perRank);
        // }
    }
    // if (!table.empty()) {
    //     ofstream bin("parents_n" + to_string(n) + "_rank" + to_string(rank) + ".bin", ios::binary | ios::app);
    //     bin.write((char*)table.data(), table.size() * sizeof(Record));
    //     bin.close();
    // }

    double t1 = MPI_Wtime();

    cout << "Rank " << rank << " did " << (end - start) << " records in " << (t1 - t0) << " s\n";

    MPI_Finalize();
    return 0;
}