// File:        lcs.cpp
// Author:      Liam Morris
// Description: Implements the functions contained in lcs.h and contains the
//              main function for executing the program.

#include "lcs.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sys/resource.h>
#include <unistd.h>

using std::cerr;
using std::cin;
using std::clock;
using std::cout;
using std::ifstream;
using std::endl;
using std::max;
using std::string;

const char* kDP = "-dp";
const char* kHirsch = "-hirsch";
const char* kNaive = "-naive";
const char* kMemo = "-memo";
const char* kUsage =
    "Usage: lcs [-naive] [-memo] [-dp] num_tests input_file";
const int kUp = 1;
const int kUpLeft = 2;
const int kLeft = 3;

string reverse(const string& s) {
    string reversed;

    for (int i = 0; i < s.length(); ++i) {
        reversed.push_back(s.at(s.length() - i - 1));
    }

    return reversed;
}

int* alg_b(int m, int n, const string& seq1, const string& seq2, int** k) {
    int* k[2];
    k[0] = new int[n + 1];
    k[1] = new int[n + 1];

    for (int i = 0; i < n + 1; ++i) {
        k[1][i] = 0;
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            k[0][j] = k[1][j];
        }

        for (int j = 0; j < n; ++j) {
            if (seq1.at(i) == seq2.at(j)) {
                k[1][j + 1] = k[0][j] + 1;
            } else {
                k[1][j + 1] = max(k[1][j], k[0][j + 1]);
            }
        }
    }

    delete[] k[0];

    return k[1];
}

string alg_c(int m, int n, const string& seq1, const string& seq2) {
    string c;

    if (n == 0) {
        c == "";
    } else if (m == 1) {
        if (seq2.find(seq1.at(0)) != string::npos) {
            c = seq1.at(0);
        } else {
            c = "";
        }
    } else {
        int i = m / 2;

        int* L1 = alg_b(i, n, seq1.substr(0, i), seq2);
        int* L2 = alg_b(m - i, n, reverse(seq1.substr(i)), reverse(seq2));

        int k = 0;
        int max = 0;
        for (int j = 0; j < n + 1; ++j) {
            if (max < L1[j] + L2[n - j]) {
                max = L1[j] + L2[n - j];
                k = j;
            }
        }

        delete[] L1;
        delete[] L2;

        c =  alg_c(i, k, seq1.substr(0, i), seq2.substr(0, k)) +
             alg_c(m - i, n - k, seq1.substr(i), seq2.substr(k));
    }

    return c;
}

string hirschberg(const string& seq1, const string& seq2) {
    return alg_c(seq1.length(), seq2.length(), seq1, seq2);
}

string lcs_dp(const string& seq1, const string& seq2) {
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


    string lcs = build_lcs(seq1, b, seq1.length(), seq2.length());
    
    for (int i = 0; i < m + 1; ++i) {
        delete b[i];
        delete c[i];
    }

    return lcs;
}

string build_lcs(const string& seq1, int** b, int i, int j) {
    if (i == 0 || j == 0) {
        return "";
    } else if (b[i][j] == kUpLeft) {
        return build_lcs(seq1, b, i - 1, j - 1) + seq1.at(i - 1);
    } else if (b[i][j] == kUp) {
        return build_lcs(seq1, b, i - 1, j);
    } else {
        return build_lcs(seq1, b, i, j - 1);
    }
}

int lcs_length_memo(const string& seq1, const string& seq2, int i, int j) {
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
    const string& seq1, const string& seq2, int i, int j, int** lengths) {
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

int lcs_length_naive(const string& seq1, const string& seq2, int i, int j) {
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
    bool use_hirschberg = false;

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
        } else if (strcmp(argv[i], kHirsch) == 0) {
            use_hirschberg = true;
            ++num_args;
        }
    }

    if (argc - num_args < 3) {
        cerr << kUsage << endl;
        return 1;
    }

    int num_tests = atoi(argv[argc - 2]);

    for (int i = 0; i < num_tests; ++i) {
        ifstream input_file(argv[argc - 1]);
        string seq1, seq2;
        getline(input_file, seq1);
        getline(input_file, seq2);
        input_file.close();
        clock_t start, end;

        // Run naive solution.
        if (use_naive) {
            start = clock();
            int solution = lcs_length_naive(seq1, seq2, 0, 0);
            end = clock();
            cout << "Naive solution: " << solution << endl
                 << "Running time: " << (end - start) << endl;
        }

        // Run memoization solution.
        if (use_memo) {
            start = clock();
            int solution = lcs_length_memo(seq1, seq2, 0, 0);
            end = clock();
            cout << "Memoization solution: " << solution << endl
                 << "Running time: " << (end - start) << endl;
        }

        // Run dynamic programming solution.
        if (use_dp) {
            start = clock();
            string solution = lcs_dp(seq1, seq2);
            end = clock();
            cout << "DP solution: " << solution << endl
                 << "Length: " << solution.length() << endl
                 << "Running time: " << (end - start) << endl;
        }

        // Run Hirschberg algorithm
        if (use_hirschberg) {
            start = clock();
            string solution = hirschberg(seq1, seq2);
            end = clock();
            cout << "Hirschberg algorithm solution: " << solution << endl
                 << "Length: " << solution.length() << endl
                 << "Running time: " << (end - start) << endl;
        }
    }

    return 0;
}
