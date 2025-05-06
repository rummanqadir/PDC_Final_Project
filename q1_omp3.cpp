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

const int n = 4;
const int BATCH_SIZE = 10000;

struct Record {
    uint32_t idx;
    array<uint32_t, n - 1> parents;
};

unordered_map<string, Record> permToRecord;

vector<unsigned char> Swap(const vector<unsigned char>& v, int x, const vector<int>& inv) {
    int n = v.size();
    int i = inv[x];
    if (i >= 0 && i+1 < n) {
        vector<unsigned char> p = v;
        swap(p[i], p[i+1]);
        return p;
    }
    return v;
}

vector<int> inverse(const vector<unsigned char>& v) {
    int n = v.size();
    vector<int> inv(n);
    for (int i = 0; i < n; ++i) {
        inv[v[i] - 1] = i;
    }
    return inv;
}

int pos(const vector<unsigned char>& v) {
    int r = 0, n = v.size();
    for (int i = n - 1; i >= 0; --i) {
        if (v[i] != i + 1) {
            return i + 1;
        }
    }
    return r;
}

vector<unsigned char> FindPosition(const vector<unsigned char>& v, int t, 
                                   const vector<unsigned char>& identity, int r, const vector<int>& inv) {
    int n = v.size();
    if (t == 2) {
        auto swapped = Swap(v, t-1, inv);
        if (swapped == identity) {
            return Swap(v, t-2, inv);
        }
    }
    if (v[n-2] == t || v[n-2] == n-1) {
        if (r > 0) return Swap(v, r-1, inv);
        return v;
    }
    return Swap(v, t-1, inv);
}

vector<unsigned char> parent1(const vector<unsigned char>& v, int t, 
                              const vector<unsigned char>& identity) {
    int n = v.size();
    auto inv = inverse(v);
    int r = pos(v);

    if (v[n - 1] == n) {
        return (t != n - 1) ? FindPosition(v, t, identity, r, inv) 
                            : Swap(v, v[n - 2] - 1, inv);
    } else if (v[n - 1] == n - 1 && v[n - 2] == n && Swap(v, n-1, inv) != identity) {
        return (t == 1) ? Swap(v, n - 1, inv)
                        : Swap(v, t - 2, inv);
    } else {
        return (v[n - 1] == t) ? Swap(v, n - 1, inv)
                               : Swap(v, t - 1, inv);
    }
}

string perm_to_string(const vector<unsigned char>& p) {
    string s;
    for (unsigned char val : p) s += to_string((int)val);
    return s;
}

vector<unsigned char> get_kth_permutation(long long k, int n) {
    vector<unsigned char> numbers(n);
    iota(numbers.begin(), numbers.end(), 1);
    vector<unsigned char> perm;
    long long fact = 1;
    for (int i = 1; i < n; ++i) fact *= i;
    k = k % (fact * n);
    for (int i = n - 1; i >= 0; --i) {
        if (i == 0) {
            perm.push_back(numbers[0]);
            break;
        }
        int index = k / fact;
        k %= fact;
        fact /= i;
        perm.push_back(numbers[index]);
        numbers.erase(numbers.begin() + index);
    }
    return perm;
}

long long factorial(int n) {
    long long res = 1;
    for (int i = 1; i <= n; ++i) res *= i;
    return res;
}

void print_perm(const vector<unsigned char>& p) {
    for (auto c : p) cout << int(c);
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

    cout << "Rank " << rank << " on " << hostname << " Start: " << start << " End: " << end << endl;

    vector<unsigned char> identity(n);
    iota(identity.begin(), identity.end(), 1);
    vector<vector<unsigned char>> permutations(end - start);

    double t0 = MPI_Wtime();
    vector<vector<unsigned char>> perms(end - start);
    for (long long k = start; k < end; ++k) {
        auto perm = get_kth_permutation(k, n);
        perms[k - start] = perm;
        Record r; r.idx = k;
        permToRecord[perm_to_string(perm)] = r;
    }


    for (uint32_t k = start; k < end; ++k) {
        auto perm = get_kth_permutation(k, n);
        permutations[k - start] = perm;
        Record rec;
        rec.idx = k;
        permToRecord[perm_to_string(perm)] = rec;
    }

    double buildTime = MPI_Wtime();

    vector<pair<string, Record*>> entries;
    entries.reserve(permToRecord.size());
    for (auto &kv : permToRecord) {
        entries.emplace_back(kv.first, &kv.second);
    }

    for (size_t i = 0; i < entries.size(); ++i) {
        auto keyStr = entries[i].first;
        auto keyVec = vector<unsigned char>();
        for (char c : keyStr) keyVec.push_back(c - '0');
        auto rec = entries[i].second;
        cout << keyStr << " " ;
        for (int t = 1; t < n; ++t) {
            auto p = parent1(keyVec, t, identity);
            string pStr = perm_to_string(p);
            auto it = permToRecord.find(pStr);
            cout << pStr << " ";
            rec->parents[t - 1] = (it != permToRecord.end() ? it->second.idx : 0u);
        }
        cout << endl;
    }

    double endTime = MPI_Wtime();

    cout << "Rank " << rank << " built map in " << (buildTime - t0) 
         << "s, computed parents in " << (endTime - buildTime) << "s.\n";

    // ✅ Print child and its parents
    for (long long k = start; k < end; ++k) {
        auto& child = perms[k - start];
        cout << "child k=" << k << ": ";
        print_perm(child);
        for (int t = 1; t < n; ++t) {
            auto pk = permToRecord[perm_to_string(child)].parents[t - 1];
            auto parent = get_kth_permutation(pk, n);
            cout << " | parent T_" << t << ": ";
            print_perm(parent);
        }
        cout << "\n";
    }

    MPI_Finalize();
    return 0;
}

