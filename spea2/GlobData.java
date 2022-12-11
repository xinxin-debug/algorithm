package spea2;

public class GlobData {
    // Spea2 算法参数
    public static int nvar;
    public static int nfunc;
    public static int popsize;
    public static int maxgen;
    public static int di;
    public static int dim;
    public static double pcross;
    public static double pmut;
    public static double[][] lim;
    public static int inf;

    public static double[][] dis;       // 种群距离数组
    public static double[][] archdis;   // 归档集距离数组
    public static double[][] disrow;    // 截断操作的辅助数组，第一行为下标，第二行为距离
    public static double[][] discol;

    // 设置参数
    public static void setPara(int fun_input, int popSize, int Gen, int varNum, int objNum) {
        // 算法参数
        popsize = popSize;
        maxgen = Gen;
        nvar = varNum;
        nfunc = objNum;
        di = 20;
        dim = 20;
        pcross = 0.9;
        inf = Integer.MAX_VALUE;

        dis = new double[2 * popsize][2 * popsize];
        archdis = new double[2 * popsize][2 * popsize];
        disrow = new double[2][2 * popsize];
        discol = new double[2][2 * popsize];

        lim = new double[nvar][2];
        for (int i = 0; i < nvar; i++) {
            lim[i][0] = 0.0;
            lim[i][1] = 1.0;
        }
        pmut = 1.0 / nvar;
    }

    // 构造方法
    private GlobData() {
    }

}
