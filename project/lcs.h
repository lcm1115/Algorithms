// File:        lcs.h
// Author:      Liam Morris
// Description: Contains algorithms for solving the length of and (for some
//              algorithms) construction of the longest common subsequence of
//              two strings.
//
//              The algorithms currently supported are:
//                  Naive recursive approach
//                  Naive recursive approach with memoization
//                  Dynamic programming approach
//
//              The algorithms that will be supported in the future are:
//                  Hirschberg's quadratic time linear space algorithm

#include <string>

// Returns the reverse of a given string.
// 's' is the string to be reversed.
// Returns the reverse of 's'.
std::string reverse(const std::string& s);

// Computes one row of the LCS matrix for the Hirschberg algorithm.
// 'm' is the length of seq1
// 'n' is the length of seq2
// 'seq1' is the first string for which LCS is being found
// 'seq2' is the second string for which LCS is being found
// Returns the row of the LCS matrix that was computed.
int* alg_b(int m, int n, const std::string& seq1, const std::string& seq2);

// Finds the LCS string of two given strings.
// 'm' is the length of seq1
// 'n' is the length of seq2
// 'seq1' is the first string for which LCS is being found
// 'seq2' is the second string for which LCS is being found
// Returns the LCS between 'seq1' and 'seq2'.
std::string alg_c(
        int m, int n, const std::string& seq1, const std::string& seq2);

// Finds the LCS of two given strings using the Hirschberg algorithm.
// 'seq1' is the first string for which LCS is being found
// 'seq2' is the second string for which LCS is being found
// Returns the LCS between 'seq1' and 'seq2'
std::string hirschberg(const std::string& seq1, const std::string& seq2);

// Finds the length of the longest common sequence of two strings by dynamic
// programming method.
// 'seq1' is the first string for which LCS is being found
// 'seq2' is the second string for which LCS is being found
// Returns the length of the LCS between seq1 and seq2.
std::string lcs_dp(const std::string& seq1, const std::string& seq2);

// Given a b table generated by lcs_dp, print out the LCS that was found.
// 'seq1' is the first string for which LCS was being found
// 'b' must be non-NULL and is the table generated by lcs_dp
// 'i' is the row at which to start searching
// 'j' is the column at which to start searching
std::string build_lcs(const std::string& seq1, int** b, int i, int j);

// Finds the length of the longest common sequence of two strings recursively
// using memoization.
// 'seq1' is the first sequence
// 'seq2' is the second sequence
// Returns the length of the longest common subsequence between the two
// sequences.
int lcs_length_memo(const std::string& seq1, const std::string& seq2);

// Helper function for lcs_length_memo.
// 'seq1' is the first sequence
// 'seq2' is the second sequence
// 'i' is the current index within the first sequence
// 'j' is the current index within the second sequence
// 'lengths' must be non-NULL and is the table in which values are saved for
// memoization.
int lcs_length_memo_sub(const std::string& seq1,
                        const std::string& seq2,
                        int i,
                        int j,
                        int** lengths);

// Finds the length of the longest common sequence of two strings recursively
// using the naive solution.
// 'seq1' is the first sequence
// 'seq2' is the second sequence
// 'i' is the current index within the first sequence
// 'j' is the current index within the second sequence
int lcs_length_naive(const std::string& seq1, const std::string& seq2, int i, int j);
