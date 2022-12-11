package tadf;

import spea2.GlobData;

import java.util.*;
import java.util.ArrayList;
import java.util.Random;


public class individual{
    int dimension;//维度
    int M;//目标数
    double fitness;//适应度
    double distance;//拥挤度
    int ri;//在UpdateCA里个体绑定的最近权重向量序号
    int rj;//在UpdateDA里个体绑定的最近权重向量序号
    int S;//存储SPEA2算法的S值。   R和D值没存在个体里，单独存储了一下
    double[] n_var;//存储每一维排序的信息
    int nq;//个体被支配的次数
    ArrayList<individual> Sp=new ArrayList<>();//个体支配个体的集合
    int rank;//个体的帕累托的等级
    double[] dec;//二维决策空间 dec[0] dec[1]
    double[] obj;//二维目标空间
    int N;

    int number;
    double high=1.0;
    double low=0.0;
    String testFunName;

    public individual(int i, String tFN, int n, int oNum, int vNum){
        this.testFunName = tFN;
        this.N = n;
        this.M = oNum;
        this.dimension = vNum;
        this.init(i);
    }


    //*******个体初始化*******
    public void init(int k){
        this.fitness=0;
        this.distance=0;
        this.ri=-1;
        this.rj=-1;
        this.S=0;
        this.n_var=new double[N];//暂定为有100个决策变量

        double[][] lim = new double[N][2];
        for (int i = 0; i < N; i++) {
            lim[i][0] = 0.0;
            lim[i][1] = 1.0;
        }
        this.dec=new double[dimension];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < dimension; j++) {
                double rd = Math.random();
                this.dec[j] = rd * (lim[j][1] - lim[j][0]) + lim[j][0];
            }
        }
        this.number=0;
        this.calculate_MMF1(testFunName);//计算该个体的目标空间值
    }


    //深拷贝构造函数
    public individual(individual that){
        this.M = that.M;
        this.dimension = that.dimension;
        this.N = that.N;
        this.testFunName = that.testFunName;
        this.fitness=that.fitness;
        this.distance=that.distance;
        this.ri=that.ri;
        this.rj=that.rj;
        this.S=that.S;
        this.n_var=new double[that.n_var.length];
        for(int i=0;i<that.n_var.length;i++) {
            this.n_var[i] = that.n_var[i];
        }
        this.nq=that.nq;
        this.Sp=that.Sp;
        this.rank=that.rank;
        this.dec=new double[that.dimension];
        this.obj=new double[that.M];
        for(int i=0;i<that.dimension;i++) {
            this.dec[i] = that.dec[i];
        }
        for(int j = 0; j < that.M; j++){
            this.obj[j] = that.obj[j];
        }
    }


    //********** MMF1 *******
    //计算函数目标值
    public void calculate_MMF1(String testFunName){ //根据决策空间计算MMF1问题的目标空间值,这里还不是适应度
        //修复Dec里无效的解(即保证决策空间的值在指定范围内)
        double sum, g;
        //double[] x = new double[dimension];
        double[] f = new double[M];

        switch (testFunName){
            case "ZDT1":
                f[0] = this.dec[0];
                sum = 0.0;

                for (int j = 1; j < dimension; j++) {
                    sum += this.dec[j];
                }
                g = 1 + 9 * (sum / (dimension - 1));

                f[1] = g * (1 - (Math.pow((double)f[0]/g, 0.5)));
                break;
            case "ZDT2":
                f[0] = this.dec[0];
                sum = 0.0;

                for (int j = 1; j < dimension; j++) {
                    sum += this.dec[j];
                }
                g = 1 + 9 * (sum / (dimension - 1));
                f[1] = g * (1 - (Math.pow(f[0] / g, 2)));
                break;
            case "ZDT3":
                f[0] = this.dec[0];
                sum = 0.0;

                for (int j = 1; j < dimension; j++) {
                    sum += this.dec[j];
                }
                g = 1 + 9 * (sum / (dimension - 1));
                f[1] = g *(1 -Math.pow((f[0]/g),0.5)-(f[0]/g)*Math.sin((10*Math.PI*f[0])));
                break;
            case "DTLZ1":
                sum = 0.0;
                for (int j = M; j < dimension; j++) {
                    sum += (Math.pow(this.dec[j] - 0.5, 2) - Math.cos(20 * Math.PI * (this.dec[j] - 0.5)));
                }

                double g1 = 100.0 * (dimension - M) + 100 * sum;
                if(M <= 2){
                    f[0] = 0.5 * (1 + g1) * this.dec[0];
                    f[1] = 0.5 * (1 + g1) * (1 - this.dec[0]);
                }else {
                    f[0] = 1;
                    for (int i1 = 0; i1 < M - 1; i1++) {
                        f[0] = f[0] * this.dec[i1];
                    }
                    f[0] = 0.5 * (1 + g1) * f[0];
                    for (int i1 = 1; i1 < M; i1++) {
                        f[i1] = 1;
                        for (int j = 0; j < M - i1 - 1; j++) {
                            f[i1] = f[i1] * this.dec[j];
                        }
                        f[i1] = f[i1] * (1 - this.dec[M - i1 - 1]) * 0.5 * (1 + g1);
                    }
                    // Y[objNum-1] = (1 + g) * (1 - X[0]);
                }
                break;
            case "DTLZ3":
                double sum2 = 0.0;
                for (int j = M; j < dimension; j++) {
                    sum2 += (Math.pow(this.dec[j] - 0.5, 2) - Math.cos(20 * Math.PI * (this.dec[j] - 0.5)));
                }

                g = 100.0 * (dimension - M) + 100 * sum2;
                if(M <= 2){
                    f[0] = (1 + g) * Math.cos( this.dec[0] * Math.PI/2);
                    f[1] = (1 + g) * Math.sin( this.dec[0] * Math.PI/2);
                }else{
                    f[0] = 1;
                    for(int i1 = 0; i1 < M - 1; i1++){
                        f[0] = f[0] * Math.cos( this.dec[i1] * Math.PI/2);
                    }
                    f[0] = (1 + g) * f[0];
                    for(int i2 = 1; i2 < M; i2++){
                        f[i2] = 1;
                        for(int j = 0; j < M - i2 - 1; j++) {
                            f[i2] = f[i2] * Math.cos( this.dec[j] * Math.PI/2);
                        }
                        f[i2] = f[i2] * (1 + g) * Math.sin( this.dec[M - i2 - 1] * Math.PI/2);
                    }
                    // Y[objNum-1] = (1 + g) * (1 - X[0]);
                }
                break;
            case "DTLZ7":
                double sum7 = 0;
                for (int j = M; j < dimension; j++) {
                    sum7 += this.dec[j];
                }
                double g7 = 1 + 9 / (dimension - M) * sum7;
                double sumH = 0;
                for(int i3 = 0; i3 < M - 1; i3++){
                    f[i3] = this.dec[i3];
                    sumH += f[i3]/(1+g7)*(1+Math.sin(3*Math.PI*f[i3]));
                }
                f[M - 1] = (1 + g7) * (M-sumH);
                break;
        }

        //计算目标空间值
        this.obj=new double[M];
        this.obj[0]=f[0];
        this.obj[1]=f[1];
    }
}
