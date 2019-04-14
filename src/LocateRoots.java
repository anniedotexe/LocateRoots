/**
 * @author      Annie Wu
 * @date        April 13, 2019
 *
 * Class:       CS 3010 - Numerical Methods
 * Project:     2
 *
 * Description: This program will locate roots of an equation using Bisection,
 *              Newton-Raphson, Secant, False-Position, and Modified Secant methods.
 *
 *              Function #1 f(x) = 2x3 – 11.7x2 + 17.7x – 5 has 3 roots in range [0, 4]
 *                  - actual roots: 0.365, 1.922, 3.563
 *              Function #2 f(x) = x + 10 – xcosh(50/x) has 1 root in range [120, 130]
 *                  - actual root: 126.632
 */

import java.io.*;

public class LocateRoots {

    private static final int MAX = 100; // max iterations
    private static final double delta = 0.01; // delta for modified secant
    private static final double ERROR = 0.01; // 1% desired approx error
    private static final double DIVERGING_ERROR = 20; // 2000% error for checking divergence
    private static PrintStream output;

    /**
     * Write iteration number and the approximate error to output file
     * @param n iteration
     * @param error approximate error
     */
    private static void writeToFile(int n, double error) {
        output.println(n + "," + error);
    }

    /**
     * Write header for Iteration and Error to output file
     */
    private static void header() {
        output.println("Iteration, Error");
    }

    /**
     * Print formatted table for bracketing methods (Bisection and False-Position) to console
     * @param n iteration
     * @param a left / lower value
     * @param b right / upper value
     * @param c current root guess
     * @param fa f(a)
     * @param fb f(b)
     * @param fc f(c)
     * @param error approximate error
     */
    private static void printBracketing(int n, double a, double b, double c, double fa, double fb,
                                        double fc, double error) {
        System.out.printf("   %d \t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\n",
                        n, a, b, c, fa, fb, fc, error);
    }

    /**
     * Print formatted table for Newton-Raphson method to console
     * @param n iteration
     * @param x current root x
     * @param fx f(x) value
     * @param fPrimeX f'(x) value
     * @param error approximate error
     */
    private static void printNewton(int n, double x, double fx, double fPrimeX, double error) {
        System.out.printf("   %d \t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f \n", n, x, fx, fPrimeX, error);
    }

    /**
     * Print formatted table for Secant and Modified Secant methods to console
     * @param n iteration
     * @param previousX xn-1
     * @param x xn
     * @param fPreviousX f(xn-1)
     * @param fx f(x)
     * @param fPrimeX f'(x)
     * @param error approximate error
     */
    private static void printSecant(int n, double previousX, double x, double fPreviousX, double fx,
                                    double fPrimeX, double error) {
        System.out.printf("   %d \t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\n",
                        n, previousX, x, fPreviousX, fx, fPrimeX, error);
    }

    private static void printModSecant(int n, double previous, double x, double fPreviousX, double fx,
                                       double fxAndDelta, double fPrimeX, double currentError) {
        System.out.printf("   %d \t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\t|   %.3f\n",
                n, previous, x, fPreviousX, fx, fxAndDelta, fPrimeX, currentError);
    }

    /**
     * Check if we found a root by checking if f(x) is close to 0
     * @param fx current x value
     * @return True if root
     */
    private static boolean isRoot(double fx){
        if (fx > -.5 && fx < .5) {
            return true;
        }
        return false;
    }

     /**
     * Get f(x) at this x value
     * @param functionNumber function 1 or 2
     * @param x current value
     * @return f(x) value
     */
    // #1. f(x) = 2x3 – 11.7x2 + 17.7x – 5
    // #2. f(x) = x + 10 – xcosh(50/x)
     private static double getFX(int functionNumber, double x) {
        if (functionNumber == 1) { // function #1
            return ((2 * x * x * x) - (11.7 * x * x) + (17.7 * x) - 5);
        }
        else { // function #2
            if (x == 0) {
                throw new IllegalArgumentException("Cannot Divide By Zero!!");
            }
            return (x + 10 - (x * Math.cosh(50/x)));
        }
    }

    /**
     * Get f'(x) at this x value
     * @param functionNumber function 1 or 2
     * @param x current value
     * @return f'(x) value
     */
    // #1. f'(x) = 6x^2 - 23.4x + 17.7
    // #2. f'(x) = 1 - cosh(50/x) + (50sinh(50/x))/x
    private static double getFPrimeX(int functionNumber, double x) {
        if (functionNumber == 1) { // function #1
            return ((6 * x * x) - (23.4 * x) + 17.7);
        }
        else { // function #2
            if (x == 0) {
                throw new IllegalArgumentException("\nERROR: Cannot Divide By Zero!!");
            }
            return (1 - Math.cosh(50/x) + ((50*Math.sinh(50/x))/x));

        }
    }

    /**
     * Check if we have ran the method for the maximum iterations (100)
     * @param n iterations
     * @return True if we hit 100
     */
    private static boolean maxIterations(int n) {
        if (n == MAX) {
            System.out.println("\nRoot Has Not Been Found after 100 iterations.");
            return true;
        }
        return false;
    }

    /**
     * Get absolute value of the approximate error for current value
     * @param current current value
     * @param previous previous value
     * @return error
     */
    private static double getError(double current, double previous) {
        return Math.abs((current - previous) / current);
    }

    /**
     * Bisection Method
     * @param functionNumber function 1 or 2
     * @param a left / lower value
     * @param b right / upper value
     */
    private static void bisection(int functionNumber, double a, double b) {
        double initialA = a;
        double initialB = b;
        double c = 0;
        double fa = 0;
        double fb = 0;
        double fc = 0;
        double previous = 0;
        double currentError = 1;
        int n = 0; // iterations

        System.out.println("\n   n \t|     a  \t|     b  \t|     c  \t|    f(a)\t|    f(b)\t|    f(c)\t|   Error ");
        System.out.println("-------------------------------------------------------------------------------------------------------------------");

        while (!maxIterations(n)) {
            c = (a + b) / 2;
            fa = getFX(functionNumber, a);
            fb = getFX(functionNumber, b);
            fc = getFX(functionNumber, c);
            currentError = getError(c, previous);

            // print to console and output file
            printBracketing(n, a, b, c, fa, fb, fc, currentError);
            writeToFile(n, currentError);

            if (currentError > DIVERGING_ERROR) {
                System.out.println("ERROR: This equation is diverging.");
                output.println("ERROR: This equation is diverging.");
                return;
            }

            if (fc == 0 || currentError < ERROR) { // found root
                // check if it is an actual root
                if (isRoot(fc)) {
                    System.out.printf("BISECTION - The root %.3f has been found in between %.0f and %.0f for " +
                            "function #%d in %d iterations.\n", c, initialA, initialB, functionNumber, n);
                    output.println("root = " + c);
                }
                else {
                    System.out.printf("BISECTION - There are no roots in between %.0f and %.0f for " +
                            "function #%d.\n", initialA, initialB, functionNumber);
                    output.println("root = DNE");
                }
                return;
            }
            else { // fc != 0
                if (fa * fc < 0) {
                    b = c;
                }
                else { // fa * fc > 0
                    a = c;
                }
            }

            previous = c;
            n++;
        }
    }

    /**
     * Newton-Raphson Method
     * @param functionNumber function 1 or 2
     * @param x initial root guess
     */
    private static void newtonRaphson(int functionNumber, double x) {
        double initialX = x;
        double next = 0;
        double fx = 0;
        double fPrimeX = 0;
        double currentError = 1;
        double previous = 0;
        int n = 0; // iterations

        System.out.println("\n   n \t|     xn  \t|   f(xn) \t|   f'(xn) \t|   Error ");
        System.out.println("------------------------------------------------------------------");

        while (!maxIterations(n)) {
            fx = getFX(functionNumber, x);
            fPrimeX = getFPrimeX(functionNumber, x);
            currentError = getError(x, previous);

            // xn+1 = xn - f(xn)/f'(xn)
            next = x - (fx / fPrimeX);

            // print to console and output file
            printNewton(n, x, fx, fPrimeX, currentError);
            writeToFile(n, currentError);

            if (currentError > DIVERGING_ERROR) {
                System.out.println("ERROR: This equation is diverging.");
                output.println("ERROR: This equation is diverging.");
                return;
            }

            if (fx == 0 || currentError < ERROR) { // found root
                // check if it is an actual root
                if (isRoot(fx)) {
                    System.out.printf("NEWTON - The root %.3f has been found for function #%d starting at x = %.0f in " +
                            "%d iterations.\n", x, functionNumber, initialX, n);
                    output.println("root = " + x);
                }
                else {
                    System.out.printf("NEWTON - There is no root for function #%d starting at " +
                            "x = %.0f.\n", functionNumber, initialX);
                    output.println("root = DNE");
                }
                return;
            }

            if (fPrimeX == 0) {
                System.out.println("ERROR: f'(xn) = 0, cannot continue finding the root.");
                output.println("ERROR: f'(xn) = 0.");
                return;
            }

            previous = x;
            x = next;
            n++;
        }
    }

    /**
     * Secant Method.
     * @param functionNumber function 1 or 2
     * @param previous previous value
     * @param x current value
     */
    private static void secant(int functionNumber, double previous, double x) {
        double initialX = x;
        double next = 0;
        double fPreviousX = 0;
        double fx = 0;
        double fPrimeX = 0;
        double currentError = 1;
        int n = 0; // iterations

        System.out.println("\n   n \t|    xn-1 \t|     xn \t|   f(xn-1)\t|   f(xn)\t|   f'(xn)\t|   Error ");
        System.out.println("----------------------------------------------------------------------------------------------------");

        while (!maxIterations(n)) {
            fx = getFX(functionNumber, x);
            fPreviousX = getFX(functionNumber, previous);
            fPrimeX = getFPrimeX(functionNumber, x);
            currentError = getError(x, previous);

            // xn+1 = xn - f(xn) * (xn - xn-1) / (f(xn) - f(xn-1))
            next = x - fx * (x - previous) / (fx - fPreviousX);

            // print to console and output file
            printSecant(n, previous, x, fPreviousX, fx, fPrimeX, currentError);
            writeToFile(n, currentError);

            if (currentError > DIVERGING_ERROR) {
                System.out.println("ERROR: This equation is diverging.");
                output.println("ERROR: This equation is diverging.");
                return;
            }

            if (fx == 0 || currentError < ERROR) { // found root
                // check if it is an actual root
                if (isRoot(fx)) {
                    System.out.printf("SECANT - The root %.3f has been found for function #%d starting at x = %.0f in " +
                            "%d iterations.\n", x, functionNumber, initialX, n);
                    output.println("root = " + x);
                }
                else {
                    System.out.printf("SECANT - There is no root for function #%d starting at " +
                            "x = %.0f.\n", functionNumber, initialX);
                    output.println("root = DNE");
                }
                return;
            }

            if (fPrimeX == 0) {
                System.out.println("ERROR: f'(xn) = 0, cannot continue finding the root.");
                output.println("ERROR: f'(xn) = 0");
                return;
            }

            previous = x;
            x = next;
            n++;
        }
    }

    /**
     * False-Position Method
     * @param functionNumber function 1 or 2
     * @param a left / lower value
     * @param b right / upper value
     */
    private static void falsePosition(int functionNumber, double a, double b) {
        double initialA = a;
        double initialB = b;
        double c = 0;
        double fa = 0;
        double fb = 0;
        double fc = 0;
        double previous = 0;
        double currentError = 1;
        int n = 0; // iterations

        System.out.println("\n   n \t|     a  \t|     b  \t|     c  \t|    f(a)\t|    f(b)\t|    f(c)\t|   Error ");
        System.out.println("-------------------------------------------------------------------------------------------------------------------");

        while (!maxIterations(n)) {
            fa = getFX(functionNumber, a);
            fb = getFX(functionNumber, b);

            // c = (af(b) - bf(a)) / (f(b) - f(a))
            c = (a*fb - b*fa) / (fb - fa);
            fc = getFX(functionNumber, c);
            currentError = getError(c, previous);

            // print to console and output file
            printBracketing(n, a, b, c, fa, fb, fc, currentError);
            writeToFile(n, currentError);

            if (currentError > DIVERGING_ERROR) {
                System.out.println("ERROR: This equation is diverging.");
                output.println("ERROR: This equation is diverging.");
                return;
            }

            if (fc == 0 || currentError < ERROR) { // found root
                // check if it is an actual root
                if (isRoot(fc)) { // a root
                    System.out.printf("FALSE-POSITION - The root %.3f has been found in between %.0f and %.0f for " +
                            "function #%d in %d iterations.\n", c, initialA, initialB, functionNumber, n);
                    output.println("root = " + c);
                } else { // not a root
                    System.out.printf("FALSE-POSITION - There is no root in between %.0f and %.0f for " +
                            "function #%d.\n", initialA, initialB, functionNumber);
                    output.println("root = DNE");
                }
                return;
            }
            else { // fc != 0
                if (fa * fc < 0) {
                    b = c;
                }
                else { // fa * fc > 0
                    a = c;
                }
            }

            previous = c;
            n++;
        }
    }

    /**
     * Modified Secant Method, with delta = 0.01
     * @param functionNumber function 1 or 2
     * @param previous previous value
     * @param x current value
     */
    private static void modifiedSecant(int functionNumber, double previous, double x) {
        double initialX = x;
        double next = 0;
        double fx = 0;
        double fPreviousX = 0;
        double fPrimeX = 0;
        double fxAndDelta = 0;
        double currentError = 1;
        int n = 0; // iterations

        System.out.println("\n   n \t|    xn-1 \t|     xn \t|   f(xn-1)\t|   f(xn)\t|  f(x+delta*x)\t|   f'(xn)\t|   Error ");
        System.out.println("-------------------------------------------------------------------------------------------------------------------");

        while (!maxIterations(n)) {
            fx = getFX(functionNumber, x);
            fPreviousX = getFX(functionNumber, previous);
            fPrimeX = getFPrimeX(functionNumber, x);
            fxAndDelta = getFX(functionNumber, (x + (delta * x)));
            currentError = getError(x, previous);

            // xn+1 = xn - (f(xn) * ((delta * xn)) / (f(xn + (delta * xn)) - f(xn)))
            next = x - (fx*(delta * x))/(fxAndDelta - fx);

            // print to console and output file
            printModSecant(n, previous, x, fPreviousX, fx, fxAndDelta, fPrimeX, currentError);
            writeToFile(n, currentError);

            if (currentError > DIVERGING_ERROR) {
                System.out.println("ERROR: This equation is diverging.");
                output.println("ERROR: This equation is diverging.");
                return;
            }

            if (fx == 0 || currentError < ERROR) { // found root
                // check if it is an actual root
                if (isRoot(fx)) {
                    System.out.printf("MODIFIED SECANT - The root %.3f has been found for function #%d starting at x = %.0f in " +
                            "%d iterations.\n", x, functionNumber, initialX, n);
                    output.println("root = " + x);
                }
                else {
                    System.out.printf("MODIFIED SECANT - There is no root for function #%d starting at " +
                            "x = %.0f.\n", functionNumber, initialX);
                    output.println("root = DNE");
                }
                return;
            }

            if (fPrimeX == 0) {
                System.out.println("ERROR: f'(xn) = 0, cannot continue finding the root.");
                output.println("ERROR: f'(xn) = 0");
                return;
            }

            previous = x;
            x = next;
            n++;
        }
    }

    /**
     * Main method
     * @param args the input arguments
     */
    public static void main(String[] args) throws FileNotFoundException {
        String fileName = "output.csv";
        FileOutputStream out = new FileOutputStream(fileName);
        output = new PrintStream(out);

        System.out.println("\n Program 2: Locating Roots of a Function - Annie Wu\n");

        // Function #1
        // Bisection
        output.println("\nBisection Function #1");
        header();
        bisection(1, 0, 1);
        header();
        bisection(1, 1, 2);
        header();
        bisection(1, 2, 3);
        header();
        bisection(1, 3, 4);


        // Newton-Raphson
        output.println("\nNewton-Raphson Function #1");
        header();
        newtonRaphson(1, 1);
        header();
        newtonRaphson(1, 2);
        header();
        newtonRaphson(1, 3);
        header();
        newtonRaphson(1, 4);

        // Secant
        output.println("\nSecant Function #1");
        header();
        secant(1, 0, 1);
        header();
        secant(1, 1, 2);
        header();
        secant(1, 2, 3);
        header();
        secant(1, 3, 4);

        // False-Position
        output.println("\nFalse-Position Function #1");
        header();
        falsePosition(1, 0, 1);
        header();
        falsePosition(1, 1, 2);
        header();
        falsePosition(1, 2, 3);
        header();
        falsePosition(1, 3, 4);

        //Modified Secant
        output.println("\nModified Secant Function #1");
        header();
        modifiedSecant(1, 0, 1);
        header();
        modifiedSecant(1, 1, 2);
        header();
        modifiedSecant(1, 2, 3);
        header();
        modifiedSecant(1, 3, 4);

        // Function #2
        output.println("\nBisection Function #2");
        header();
        bisection(2, 120, 130);
        output.println("\nNewton-Raphson Function #2");
        header();
        newtonRaphson(2, 130);
        output.println("\nSecant Function #2");
        header();
        secant(2, 120, 130);
        output.println("\nFalse Position Function #2");
        header();
        falsePosition(2, 120, 130);
        output.println("\nModified Secant Function #2");
        header();
        modifiedSecant(2, 120, 130);

        output.close();

        System.out.println("\nEnd of Program\n");
    }
}
