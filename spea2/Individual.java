package spea2;

// 个体类
public class Individual {

    double[] xreal;     //  决策变量
    double[] fvalue;    // 函数值
    double fitness;     // 适应度

    public Individual() {
        xreal = new double[GlobData.nvar];     // 决策变量维度
        fvalue = new double[GlobData.nfunc];  // 目标函数个数
    }
}
