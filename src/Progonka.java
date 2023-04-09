public class Progonka {

    public static void main(String[] args) {
        double[][] matrix = {{2, 4, 0}, {4, 1, 5}, {0, 5, 2}};
        double[] constants = {20, 37, 30};
        double[] solutions = solveTridiagonalMatrix(matrix, constants);

        System.out.println("Solutions:");
        for (int i = 0; i < solutions.length; i++) {
            System.out.println("x" + (i+1) + " = " + solutions[i]);
        }
    }

    public static double[] solveTridiagonalMatrix(double[][] matrix, double[] constants) {
        int n = constants.length;
        double[] x = new double[n];
        double[] a = new double[n-1];
        double[] b = new double[n];
        double[] c = new double[n-1];
        double[] d = new double[n];

        // Decompose the matrix into tridiagonal form
        for (int i = 0; i < n-1; i++) {
            a[i] = matrix[i+1][i];
            b[i] = matrix[i][i];
            c[i] = matrix[i][i+1];
        }
        b[n-1] = matrix[n-1][n-1];

        for (int i = 1; i < n; i++) {
            double m = a[i-1] / b[i-1];
            b[i] -= m * c[i-1];
            constants[i] -= m * constants[i-1];
        }

        x[n-1] = constants[n-1] / b[n-1];
        for (int i = n-2; i >= 0; i--) {
            x[i] = (constants[i] - c[i] * x[i+1]) / b[i];
        }

        return x;
    }
}
