public class Jacobi {
    static final int N = 4;
    static double det = 1;
    static double[] X = new double[N];
    static double[] prevX = new double[N];
    static double epsilon = 0.0001;
    static double[][] A = {{4, 0, 1, 1}, {0, 3, 0, 1}, {1, 0, 2, 0}, {1, 1, 0, 5}};
    static double[] B = {11, 10, 7, 23};
    static double[][] AB = new double[N][N + 1];
    static int[][] P = new int[N][N];
    static double[][] M = new double[N][N];
    public static void main(String[] args) {
        System.out.println("\n\nМетод Якобі:");
        calculateJacobi();
    }

    private static void calculateJacobi() {
        X = new double[N];
        prevX = new double[N];

        double norm = 0;
        int iterations = 0;

        do {
            iterations++;
            for (int i = 0; i < N; i++) {
                double sum = 0;
                for (int j = 0; j < N; j++) {
                    if (j != i) {
                        sum += A[i][j] * prevX[j];
                    }
                }
                X[i] = (B[i] - sum) / A[i][i];
            }
            norm = calculateNorm(X, prevX);
            if (iterations == 1) System.out.println("Норма на першій ітерації: " + norm);
            System.arraycopy(X, 0, prevX, 0, X.length);

        } while (norm > epsilon);

        System.out.println("Розв'язки системи:");
        for (int i = 0; i < X.length; i++) {
            System.out.println("x" + (i + 1) + " = " + X[i]);
        }

        System.out.println("Норма на останній ітерації: " + norm);
    }


    public static double calculateNorm(double[] X, double[] prevX) {
        double norm = 0;
        for (int i = 0; i < X.length; i++) {
            norm += Math.pow(X[i] - prevX[i], 2);
        }
        return Math.pow(norm,  (double) 1 / X.length);
    }


    private static void calculateNumberOf() {
        System.out.println("\n\n");
        double conditionNumber = conditionNumber(A);
        System.out.println("Число обумовленості: " + conditionNumber);
    }

    public static double conditionNumber(double[][] A) {
        double normA = matrixNorm(A, "A");
        double normAInv = matrixNorm(inverse(A), "invA");
        return normA * normAInv;
    }


    public static double matrixNorm(double[][] A, String name) {
        double maxEigenvalue = maxEigenvalue(A);
        double norm = Math.sqrt(maxEigenvalue);
        System.out.println("Норма матриці " + name + ": " + norm);
        return norm;
    }


    public static double[][] inverse(double[][] A) {
        int n = A.length;
        double[][] AInv = new double[n][n];
        double det = determinant(A);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                AInv[i][j] = Math.pow(-1, i + j) * determinant(minor(A, i, j));
            }
        }
        AInv = transpose(AInv);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                AInv[i][j] /= det;
            }
        }
        System.out.println("Обернена до А матриця: ");
        printMatrix(AInv);
        return AInv;
    }


    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double element : row) {
                System.out.printf("%10.5f", element);
            }
            System.out.println();
        }
    }


    public static double determinant(double[][] A) {
        int n = A.length;
        if (n == 1) {
            return A[0][0];
        } else if (n == 2) {
            return A[0][0] * A[1][1] - A[0][1] * A[1][0];
        } else {
            double det = 0.0;
            for (int j = 0; j < n; j++) {
                det += Math.pow(-1, j) * A[0][j] * determinant(minor(A, 0, j));
            }
            return det;
        }
    }

    public static double[][] minor(double[][] A, int i, int j) {
        int n = A.length;
        double[][] M = new double[n - 1][n - 1];
        for (int k = 0, x = 0; k < n; k++) {
            if (k == i) {
                continue;
            }
            for (int l = 0, y = 0; l < n; l++) {
                if (l == j) {
                    continue;
                }
                M[x][y] = A[k][l];
                y++;
            }
            x++;
        }
        return M;
    }


    public static double maxEigenvalue(double[][] A) {
        int n = A.length;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = Math.random();
        }
        double lambda = 0.0;
        double eps = 0.0001;
        while (true) {
            for (int i = 0; i < n; i++) {
                y[i] = 0.0;
                for (int j = 0; j < n; j++) {
                    y[i] += A[i][j] * x[j];
                }
            }
            double maxAbsValue = Math.abs(y[0]);
            for (int i = 1; i < n; i++) {
                if (Math.abs(y[i]) > maxAbsValue) {
                    maxAbsValue = Math.abs(y[i]);
                }
            }
            for (int i = 0; i < n; i++) {
                x[i] = y[i] / maxAbsValue;
            }
            double newLambda = 0.0;
            for (int i = 0; i < n; i++) {
                newLambda += y[i] * x[i];
            }
            if (Math.abs(newLambda - lambda) < eps) {
                break;
            }
            lambda = newLambda;
        }
        return lambda;
    }

    public static double[][] transpose(double[][] A) {
        int numRows = A.length;
        int numCols = A[0].length;
        double[][] result = new double[numCols][numRows];

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                result[j][i] = A[i][j];
            }
        }

        return result;
    }

}
