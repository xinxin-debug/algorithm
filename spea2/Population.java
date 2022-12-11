package spea2;

import java.util.Arrays;

// 种群类
public class Population {
    int[] S;    // 种群个体的强度
    int[] R;    // 种群个体的原始适应度
    double[] D; //
    double[] fitness;

    int[] domnum;       // 下标表示个体编号，数组的值为支配该个体的个体总数
    int[] bedomnum;     // 下标表示个体编号，数组的值表示被其他个体支配的个数
    int[][] bedomarr;   // 行表示个体编号，列表示支配该个体的所有个体的编号

    Individual[] ind;   // 种群个体

    // 构造方法
    public Population() {
        S = new int[GlobData.popsize];
        R = new int[GlobData.popsize];
        D = new double[GlobData.popsize];
        fitness = new double[GlobData.popsize];

        domnum = new int[GlobData.popsize];
        bedomarr = new int[GlobData.popsize][GlobData.popsize];
        bedomnum = new int[GlobData.popsize];

        ind = new Individual[GlobData.popsize];

        for (int i = 0; i < GlobData.popsize; i++) {
            ind[i] = new Individual();
        }
    }

    // 种群的初始化
    public void initPop() {

        for (int i = 0; i < GlobData.popsize; i++) {
            for (int j = 0; j < GlobData.nvar; j++) {
                double rd = Math.random();
                ind[i].xreal[j] = rd * (GlobData.lim[j][1] - GlobData.lim[j][0]) + GlobData.lim[j][0];
            }
        }
    }

    // 计算函数值——————目标函数在这里修改
    public void func(int fun_input) {

        double sum, g;
        double[] x = new double[GlobData.nvar];
        double[] f = new double[GlobData.nfunc];

        // 测试函数 ————ZDT1
        for (int i = 0; i < GlobData.popsize; i++) {
            System.arraycopy(ind[i].xreal, 0, x, 0, GlobData.nvar);

            // 计算函数值
            if (fun_input == 1) {       // ZDT1
                f[0] = x[0];
                sum = 0.0;

                for (int j = 1; j < GlobData.nvar; j++) {
                    sum += x[j];
                }
                g = 1 + 9 * (sum / (GlobData.nvar - 1));
                f[1] = g * (1 - (Math.pow(f[0] / g, 0.5)));
            } else if(fun_input == 2){ //ZDT2
                f[0] = x[0];
                sum = 0.0;

                for (int j = 1; j < GlobData.nvar; j++) {
                    sum += x[j];
                }
                g = 1 + 9 * (sum / (GlobData.nvar - 1));
                f[1] = g * (1 - (Math.pow(f[0] / g, 2)));
            }else if (fun_input == 3){     //ZDT3
                f[0] = x[0];
                sum = 0.0;

                for (int j = 1; j < GlobData.nvar; j++) {
                    sum += x[j];
                }
                g = 1 + 9 * (sum / (GlobData.nvar - 1));
                f[1] = g *(1 -Math.pow((f[0]/g),0.5)-(f[0]/g)*Math.sin((10*Math.PI*f[0])));
            }else if (fun_input == 4) {    // DTLZ1
                sum = 0.0;
                for (int j = GlobData.nfunc; j < GlobData.nvar; j++) {
                    sum += (Math.pow(x[j] - 0.5, 2) - Math.cos(20 * Math.PI * (x[j] - 0.5)));
                }

                double g1 = 100.0 * (GlobData.nvar - GlobData.nfunc) + 100 * sum;
                if(GlobData.nfunc <= 2){
                    f[0] = 0.5 * (1 + g1) * x[0];
                    f[1] = 0.5 * (1 + g1) * (1 - x[0]);
                }else{
                    f[0] = 1;
                    for(int i1 = 0 ; i1 < GlobData.nfunc - 1; i1++){
                        f[0] = f[0] * x[i1];
                    }
                    f[0] = 0.5 * (1 + g1) * f[0];
                    for(int i1 = 1; i1 < GlobData.nfunc; i1++){
                        f[i1] = 1;
                        for(int j = 0; j < GlobData.nfunc - i1 - 1;j++) {
                            f[i1] = f[i1] * x[j];
                        }
                        f[i1] = f[i1] * (1 - x[GlobData.nfunc - i1 - 1]) * 0.5 * (1 + g1);
                    }
                    // Y[objNum-1] = (1 + g) * (1 - X[0]);
                }
            } else if (fun_input == 5) { //DTLZ3
                double sum2 = 0.0;
                for (int j = GlobData.nfunc; j < GlobData.nvar; j++) {
                    sum2 += (Math.pow(x[j] - 0.5, 2) - Math.cos(20 * Math.PI * (x[j] - 0.5)));
                }

                g = 100.0 * (GlobData.nvar - GlobData.nfunc) + 100 * sum2;
                if(GlobData.nfunc <= 2){
                    f[0] = (1 + g) * Math.cos( x[0] * Math.PI/2);
                    f[1] = (1 + g) * Math.sin( x[0] * Math.PI/2);
                }else{
                    f[0] = 1;
                    for(int i1 = 0 ; i1 < GlobData.nfunc - 1; i1++){
                        f[0] = f[0] * Math.cos( x[i1] * Math.PI/2);
                    }
                    f[0] = (1 + g) * f[0];
                    for(int i2 = 1; i2 < GlobData.nfunc; i2++){
                        f[i2] = 1;
                        for(int j = 0; j < GlobData.nfunc - i2 - 1;j++) {
                            f[i2] = f[i2] * Math.cos( x[j] * Math.PI/2);
                        }
                        f[i2] = f[i2] * (1 + g) * Math.sin( x[GlobData.nfunc - i2 - 1] * Math.PI/2);
                    }
                    // Y[objNum-1] = (1 + g) * (1 - X[0]);
                }
            }else if (fun_input == 6){   //DTLZ7
                double sum7 = 0;
                for (int j = GlobData.nfunc; j < GlobData.nvar; j++) {
                    sum7 += x[j];
                }
                double g7 = 1 + 9 / (GlobData.nvar - GlobData.nfunc) * sum7;
                double sumH = 0;
                for(int i3 = 0 ; i3 < GlobData.nfunc - 1; i3++){
                    f[i3] = x[i3];
                    sumH += f[i3]/(1+g7)*(1+Math.sin(3*Math.PI*f[i3]));
                }
                f[GlobData.nfunc - 1] = (1 + g7) * (GlobData.nfunc-sumH);
            }

            // 更新个体的函数值
            System.arraycopy(f, 0, ind[i].fvalue, 0, GlobData.nfunc);
        }
    }

    // 计算初始种群个体的适应度
    public void calFitness() {
        int ret;
        int q;
        int ind;
        int k;
        double sigma;

        for (int i = 0; i < GlobData.popsize; i++) {
            S[i] = 0;
            R[i] = 0;
            D[i] = 0.0;
        }

        // 首先需要计算个体的强度 S(i)，即受 i 支配的解的数量
        for (int i = 0; i < GlobData.popsize; i++) {
            q = 0;
            for (int j = 0; j < GlobData.popsize; j++) {
                ret = dominate(this.ind[i], this.ind[j]);
                if (ret == 1)       // i 支配 j
                    S[i]++;
                else if (ret == 2) {
                    bedomarr[i][q] = j;
                    q++;
                }
            }
            bedomnum[i] = q;
        }

        // 然后计算个体的原始适应度 R(i)
        for (int i = 0; i < GlobData.popsize; i++) {
            for (int j = 0; j < bedomnum[i]; j++) {
                ind = this.bedomarr[i][j];
                R[i] += S[ind];
            }
        }

        // 接着计算个体的密度值 D(i) 和 适应度 F(i)
        for (int i = 0; i < GlobData.popsize; i++) {
            for (int j = i; j < GlobData.popsize; j++) {
                if (i == j)
                    GlobData.dis[i][j] = GlobData.inf;
                else {
                    GlobData.dis[i][j] = 0;
                    for (int l = 0; l < GlobData.nfunc; l++) {
                        GlobData.dis[i][j] += Math.pow((this.ind[i].fvalue[l] - this.ind[j].fvalue[l]), 2.0);
                    }
                    GlobData.dis[i][j] = Math.pow(GlobData.dis[i][j], 0.5);
                    GlobData.dis[j][i] = GlobData.dis[i][j];
                }
            }
        }

        // 最后对距离数组进行排序，获取第 k 近邻个体的距离 Sigma，再计算 D(i)、F(i)
        k = (int) Math.sqrt(GlobData.popsize * 2);
        for (int i = 0; i < GlobData.popsize; i++) {
            Arrays.sort(GlobData.dis[i]);
            sigma = GlobData.dis[i][k - 1];
            D[i] = 1.0 / (sigma + 2.0);
            fitness[i] = R[i] + D[i];
        }

        for (int i = 0; i < GlobData.popsize; i++) {
            this.ind[i].fitness = fitness[i];
        }

    }

    // 个体间的支配关系
    public int dominate(Individual y1, Individual y2) {
        int greater = 0, equal = 0, less = 0;

        for (int i = 0; i < GlobData.nfunc; i++) {
            if (y1.fvalue[i] < y2.fvalue[i])
                less++;
            else if (Math.abs(y1.fvalue[i] - y2.fvalue[i]) <= 1e-6)
                equal++;
            else
                greater++;
        }

        if (greater == 0 && equal != GlobData.nfunc)
            return 1;   // y1 支配 y2
        else if (less == 0 && equal != GlobData.nfunc)
            return 2;   // y2 支配 y1
        else
            return 3;   // 没有支配关系
    }

}
