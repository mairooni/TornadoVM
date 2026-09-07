package uk.ac.manchester.tornado.examples;

public class DebuggingClass {

    public static int[] poopyadd(int[] a, int[] b, int[] c) {
        for(int i = 0; i < c.length; i++) {
            c[i] = a[i] + b[i];
        }
        return c;
    }

    public static void main() {
        int[] a = new int[10];
        int[] b = new int[10];
        int[] c = new int[10];

        for(int i = 0; i < 10; i++) {
            a[i] = i;
            b[i] = 2;
        }

        c = poopyadd(a, b, c);

        for(int i = 0; i < 10; i++) {
            System.out.println("c[" + i + "] = " + c[i]);

        }
    }

}
