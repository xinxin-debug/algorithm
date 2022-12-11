package moead;

public class GlobData {
    // 算法参数
    public static int nvar;
    public static int nfunc;
    public static int popsize;
    public static int maxgen;
    public static int di;
    public static int dim;
    public static double pcross;
    public static double pmut;
    public static double[][] lim;
    public static int T;                // 邻近权向量的个数，即每个权向量领域内权向量的个数

    public static int[][] B;            // 每个权向量的相邻矩阵 对 i=1,2...N, B(i) = {i1, ... iT}
    public static double[][] lamda;     // N个用于定义切比雪夫子问题的均匀分布的权重向量，维度为目标函数的个数  即：为每个决策变量附上一个权值
    public static double[] Z;           // 当前目标函数的最优值（最大值）(参考点）

    public static int[] disindex;       // 计算两个权向量之间的距离时用到的临时数组的下标矩阵
    public static double[] dis;         // 计算两个权向量之间的距离时用到的临时数组


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
        T = 20;

        disindex = new int[popsize];
        dis = new double[popsize];
        lim = new double[nvar][2];
        for (int i = 0; i < nvar; i++) {
            lim[i][0] = 0.0;
            lim[i][1] = 1.0;
        }

        pmut = 1.0 / nvar;
        B = new int[popsize][T];
        lamda = new double[popsize][nfunc];
        Z = new double[nfunc];

    }

    // 构造方法
    private GlobData() {
    }

}
