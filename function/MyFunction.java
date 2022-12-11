package function;

abstract public class MyFunction {
    abstract public double[] objectives(double[] X);
    abstract public double[] constraints(double[] X);
}
