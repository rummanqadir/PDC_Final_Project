#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <numeric>
#include <mpi.h>
#include <unordered_map>
#include <omp.h>
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