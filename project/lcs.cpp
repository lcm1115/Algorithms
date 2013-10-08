#include "lcs.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>

using std::cerr;
using std::cin;
using std::clock;
using std::cout;
using std::ifstream;
using std::endl;
using std::max;
using std::string;

const char* kNaive = "-naive";
const char* kMemo = "-memo";
const char* kDP = "-dp";
const char* kUsage =
    "Usage: lcs [-naive] [-memo] [-dp] num_tests sequence1 sequence2";
const int kUp = 1;
const int kUpLeft = 2;
const int kLeft = 3;

int lcs_length_dp(const string seq1, const string seq2) {
    int m = seq1.length();
    int n = seq2.length();
    int result = 0;
    int* b[m + 1];
    int* c[m + 1];

    for (int i = 0; i < m + 1; ++i) {
        c[i] = new int[n + 1];
        b[i] = new int[n + 1];
    }

    for (int i = 1; i <= m; ++i) {
        c[i][0] = 0;
        b[i][0] = 0;
    }

    for (int i = 0; i <= n; ++i) {
        c[0][i] = 0;
        b[0][i] = 0;
    }

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (seq1.at(i - 1) == seq2.at(j - 1)) {
                c[i][j] = c[i - 1][j - 1] + 1;
                b[i][j] = kUpLeft;
            } else if (c[i - 1][j] >= c[i][j - 1]) {
                c[i][j] = c[i - 1][j];
                b[i][j] = kUp;
            } else {
                c[i][j] = c[i][j - 1];
                b[i][j] = kLeft;
            }
        }
    }

    return c[m][n];
}

void print_lcs(const string seq1, int** b, int i, int j) {
    if (i == 0 || j == 0) {
        return;
    } else if (b[i][j] == kUpLeft) {
        print_lcs(seq1, b, i - 1, j - 1);
        cout << seq1.at(i - 1);
    } else if (b[i][j] == kUp) {
        print_lcs(seq1, b, i - 1, j);
    } else {
        print_lcs(seq1, b, i, j - 1);
    }
}

int lcs_length_memo(string seq1, string seq2, int i, int j) {
    int m = seq1.length();
    int n = seq2.length();
    int result = 0;
    int* lengths[m];

    // Initialize memoization table.
    for (int k = 0; k < m; ++k) {
        lengths[k] = new int[n];
        for (int l = 0; l < n; ++l) {
            lengths[k][l] = -1;
        }
    }

    result = lcs_length_memo_sub(seq1, seq2, i, j, lengths);

    for (int k = 0; k < m; ++k) {
        delete[] lengths[k];
    }

    return result;
}

int lcs_length_memo_sub(
    string seq1, string seq2, int i, int j, int** lengths) {
    if (i == seq1.length() || j == seq2.length()) {
        return 0;
    } else if (lengths[i][j] != -1) {
        return lengths[i][j];
    } else if (seq1.at(i) == seq2.at(j)) {
        lengths[i][j] =
            1 + lcs_length_memo_sub(seq1, seq2, i + 1, j + 1, lengths);
        return lengths[i][j];
    } else {
        lengths[i][j] = 
            max(lcs_length_memo_sub(seq1, seq2, i, j + 1, lengths),
                lcs_length_memo_sub(seq1, seq2, i + 1, j, lengths));
        return lengths[i][j];
    }
}

int lcs_length_naive(string seq1, string seq2, int i, int j) {
    if (i == seq1.length() || j == seq2.length()) {
        return 0;
    } else if (seq1.at(i) == seq2.at(j)) {
        return 1 + lcs_length_naive(seq1, seq2, i + 1, j + 1);
    } else {
        return max(lcs_length_naive(seq1, seq2, i, j + 1),
                   lcs_length_naive(seq1, seq2, i + 1, j));
    }
}

int main(int argc, char** argv) {
    bool use_naive = false;
    bool use_memo = false;
    bool use_dp = false;

    int num_args = 0;

    for (int i = 1; i < argc; ++i) {
        if (i < (argc - 2) && !strchr(argv[i], '-')) {
            cerr << kUsage << endl;
            return 1;
        }
        if (strcmp(argv[i], kNaive) == 0) {
            use_naive = true;
            ++num_args;
        } else if (strcmp(argv[i], kMemo) == 0) {
            use_memo = true;
            ++num_args;
        } else if (strcmp(argv[i], kDP) == 0) {
            use_dp = true;
            ++num_args;
        }
    }

    if (argc - num_args < 3) {
        cerr << kUsage << endl;
        return 1;
    }

    int num_tests = atoi(argv[argc - 2]);

    ifstream input_file(argv[argc - 1]);
    string seq1, seq2;
    getline(input_file, seq1);
    getline(input_file, seq2);
    int solution;
    clock_t start, end;

    // Run naive solution.
    if (use_naive) {
        start = clock();
        solution = lcs_length_naive(seq1, seq2, 0, 0);
        end = clock();
        cout << "Naive solution: " << solution << endl
             << "Running time: " << (end - start) << endl;
    }

    // Run memoization solution.
    if (use_memo) {
        start = clock();
        for (int i = 0; i < 10000; ++i)
        solution = lcs_length_memo(seq1, seq2, 0, 0);
        end = clock();
        cout << "Memoization solution: " << solution << endl
             << "Running time: " << (end - start) << endl;
    }

    // Run dynamic programming solution.
    if (use_dp) {
        start = clock();
        for (int i = 0; i < num_tests; ++i) {
            solution = lcs_length_dp(seq1, seq2);
        }
        end = clock();
        cout << "DP solution: " << solution << endl
             << "Running time: " << (end - start) << endl;
    }

    return 0;
}
