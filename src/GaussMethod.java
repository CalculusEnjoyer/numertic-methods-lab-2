import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

class GaussMethod {

    static final int N = 4;
    static double det = 1;
    static double[] X = new double[N];
    static List<Double> X_correct = new ArrayList<>();
    static double[][] A = {{5, 2, 1, 0}, {1, 3, 2, 8}, {4, -6, 1, 0}, {5, 0, 3, 2}};
    static double[] B = {14, 65, -3, 32};
    static double[][] AB = new double[N][N + 1];
    static int[][] P = new int[N][N];
    static double[][] M = new double[N][N];

    public static void main(String[] args) {
        calculateGaus();
    }

    public static void calculateGaus() {
        generateAbFromAandB(A, B);
        printMatrix(AB);
        initializeX();
        for (int k = 0; k < N; k++) {
            initializeMmatrix();
            initializePmatrix();
            System.out.println("=================== " + (k + 1) + " ітерація ===================");
            int maxRow = k;
            int maxCol = k;
            double maxElement = AB[k][k];
            for (int i = k; i < N; i++) {
                for (int j = k; j < N; j++) {
                    if (Math.abs(AB[i][j]) > Math.abs(AB[maxRow][maxCol])) {
                        maxRow = i;
                        maxCol = j;
                        maxElement = AB[i][j];
                    }
                }
            }

            // Перестановка рядків
            if (maxRow != k) {
                swapRows(AB, maxRow, k);
                swapRows(P, maxRow, k);

                System.out.println("Матриця перестановки рядків: ");
                printMatrix(P);

                det *= -1;
            }

            // Перестановка стовпців
            if (maxCol != k) {
                swapCols(AB, maxCol, k);
                swapCols(P, maxCol, k);

                swapXToCorrectPosition(k, maxCol);

                System.out.println("Матриця перестановки стовпців: ");
                printMatrix(P);

                det *= -1;
            }

            //Обраховування М матриці
            M[k][k] = 1 / AB[k][k];
            for (int i = k + 1; i < N; i++) {
                M[i][k] = -AB[i][k] * M[k][k];
            }

            AB = matrixMultiplication(M, AB);

            showMatrixM();
            System.out.println("Матриця AB");
            printMatrix(AB);

            det *= maxElement;
        }

        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += AB[i][j] * X[j];
            }
            X[i] = (AB[i][N] - sum) / AB[i][i];
        }

        System.out.println("Розв'язки системи: ");
        for (int i = 0; i < N; i++) {
            System.out.println("x" + (i + 1) + " = " + X[X_correct.indexOf((double) i + 1)]);
        }

        System.out.println("Детермінант матриці коефіцієнтів обрахований прямим шляхом: " + det);
        System.out.println("Обернена матриця: ");
        inverse(A);
    }


    private static void swapXToCorrectPosition(int k, int maxCol) {
        double temp = X_correct.get(maxCol);
        X_correct.set(maxCol, X_correct.get(k));
        X_correct.set(k, temp);
    }

    private static void initializeX() {
        for (int i = 0; i < N; i++) {
            X_correct.add(i, (double) (i + 1));
        }
    }

    private static void swapRows(double[][] A, int row1, int row2) {
        double[] temp = A[row1];
        A[row1] = A[row2];
        A[row2] = temp;
    }

    private static void swapRows(int[][] A, int row1, int row2) {
        int[] temp = A[row1];
        A[row1] = A[row2];
        A[row2] = temp;
    }


    public static void swapCols(double[][] A, int col1, int col2) {
        for (int i = 0; i < A.length; i++) {
            double temp = A[i][col1];
            A[i][col1] = A[i][col2];
            A[i][col2] = temp;
        }
    }

    public static void swapCols(int[][] A, int col1, int col2) {
        for (int i = 0; i < A.length; i++) {
            int temp = A[i][col1];
            A[i][col1] = A[i][col2];
            A[i][col2] = temp;
        }
    }


    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double element : row) {
                System.out.printf("%10.5f", element);
            }
            System.out.println();
        }
    }

    public static void printMatrix(int[][] matrix) {
        for (int[] row : matrix) {
            for (int element : row) {
                System.out.print(element + " ");
            }
            System.out.println();
        }
    }

    private static void showMatrixB() {
        System.out.println("Матриця B");
        System.out.println(Arrays.toString(B));
    }

    private static void showMartixA() {
        System.out.println("Матриця A");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.print(String.format("%7.2f ", A[i][j]));
            }
            System.out.println();
        }
        System.out.println();
    }

    private static void generateAbFromAandB(double[][] A, double[] B) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, AB[i], 0, n);
            AB[i][n] = B[i];
        }
    }

    private static void initializePmatrix() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                P[i][j] = i == j ? 1 : 0;
            }
        }
    }

    private static void initializeMmatrix() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                M[i][j] = i == j ? 1 : 0;
            }
        }
    }


    public static double[][] matrixMultiplication(double[][] a, double[][] b) {
        int n = a.length;
        int m = b[0].length;
        int p = b.length;

        double[][] c = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                for (int k = 0; k < p; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return c;
    }


    private static void showMatrixM() {
        System.out.println("Матриця обнулення M");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.print(String.format("%7.2f ", M[i][j]));
            }
            System.out.println();
        }
    }


    public static double[][] inverse(double[][] A) {
        int n = A.length;
        double[][] AInv = new double[n][n];
        double det = Jacobi.determinant(A);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                AInv[i][j] = Math.pow(-1, i + j) * Jacobi.determinant(Jacobi.minor(A, i, j));
            }
        }
        AInv = Jacobi.transpose(AInv);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                AInv[i][j] /= det;
            }
        }
        System.out.println("Обернена до А матриця: ");
        printMatrix(AInv);
        return AInv;
    }
}