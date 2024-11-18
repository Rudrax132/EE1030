#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to multiply two matrices A and B, result stored in C
void multiply_matrices(double* A, double* B, double* C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = 0.0;
            for (int k = 0; k < n; k++) {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

// Function to compute the QR decomposition of a matrix A (result in Q and R)
void qr_decomposition(double* A, double* Q, double* R, int n) {
    for (int i = 0; i < n; i++) {
        // Copy the i-th column of A into the i-th column of Q
        for (int j = 0; j < n; j++) {
            Q[j * n + i] = A[j * n + i];
        }

        // Compute the norm of the i-th column of Q
        double norm = 0.0;
        for (int j = 0; j < n; j++) {
            norm += Q[j * n + i] * Q[j * n + i];
        }
        norm = sqrt(norm);

        // Normalize the i-th column of Q
        for (int j = 0; j < n; j++) {
            Q[j * n + i] /= norm;
        }

        // Set R[i][i] to the norm of the i-th column of A
        R[i * n + i] = norm;

        // Orthogonalize the remaining columns of Q
        for (int j = i + 1; j < n; j++) {
            double dot_product = 0.0;
            for (int k = 0; k < n; k++) {
                dot_product += Q[k * n + i] * A[k * n + j];
            }
            R[i * n + j] = dot_product;

            // Subtract the projection of the j-th column of A onto the i-th column of Q
            for (int k = 0; k < n; k++) {
                A[k * n + j] -= dot_product * Q[k * n + i];
            }
        }
    }
}

// Function to perform the QR algorithm to find eigenvalues
void qr_algorithm(double* A, int n, int max_iterations) {
    double* Q = (double*)malloc(n * n * sizeof(double));
    double* R = (double*)malloc(n * n * sizeof(double));
    double* A_next = (double*)malloc(n * n * sizeof(double));

    for (int iter = 0; iter < max_iterations; iter++) {
        // Perform QR decomposition of A
        qr_decomposition(A, Q, R, n);

        // Compute A_next = R * Q
        multiply_matrices(R, Q, A_next, n);

        // Copy A_next to A
        for (int i = 0; i < n * n; i++) {
            A[i] = A_next[i];
        }
    }

    // Print the eigenvalues (diagonal elements of A)
    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%.2f\n", A[i * n + i]);
    }

    free(Q);
    free(R);
    free(A_next);
}

int main() {
    int n;

    // Input matrix size
    printf("Enter the size of the matrix (n): ");
    scanf("%d", &n);

    // Allocate memory for the matrix (single dimensional array)
    double* A = (double*)malloc(n * n * sizeof(double));

    // Input matrix elements
    printf("Enter the matrix elements (row by row):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i * n + j]);
        }
    }
    //Input the number of iterations for QR algorithm
    int itr;
    printf("Enter the number of iterations for QR algorithm: ");
    scanf("%d",&itr);
    // Perform QR algorithm with a maximum of given number of iterations
    qr_algorithm(A, n, itr);

    // Free the allocated memory
    free(A);

    return 0;
}


