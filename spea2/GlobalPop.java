package spea2;

import java.util.Arrays;

public class GlobalPop {
    int[] S;        // 种群个体的强度
    int[] R;        // 种群个体的原始适应度
    double[] D;     // 密度
    double[][] fitness;     // 第一行存放下标，第二行存放距离值

    int[] domnum;       // 下标表示个体编号，数组的值为支配该个体的个体总数
    int[] bedomnum;     // 下标表示个体编号，数组的值表示被其他个体支配的个数
    int[][] bedomarr;   // 行表示个体编号，列表示支配该个体的所有个体的编号
    double[][] globdis; // 组合种群距离数组

    int[] flag;         // 标记哪些个体加入了新种群
    Individual[] ind;   // 种群个体

    // 构造方法
    public GlobalPop() {
        S = new int[2 * GlobData.popsize];
        R = new int[2 * GlobData.popsize];
        D = new double[2 * GlobData.popsize];
        fitness = new double[2][2 * GlobData.popsize];

        domnum = new int[2 * GlobData.popsize];
        bedomnum = new int[2 * GlobData.popsize];
        bedomarr = new int[2 * GlobData.popsize][2 * GlobData.popsize];
        globdis = new double[2 * GlobData.popsize][2 * GlobData.popsize];

        flag = new int[2 * GlobData.popsize];
        ind = new Individual[2 * GlobData.popsize];

        for (int i = 0; i < 2 * GlobData.popsize; i++) {
            ind[i] = new Individual();
        }
    }

    // 计算组合种群的适应度
    public void globFitness() {
        int ret = 0;
        int q = 0;
        int ind;
        int k;            // 第 k 近邻个体
        double sigma;
        int gpopsize = 2 * GlobData.popsize;

        for (int i = 0; i < gpopsize; i++) {
            S[i] = 0;
            R[i] = 0;
            D[i] = 0.0;
        }

        // 首先计算个体的强度 S(i)，即受 i 支配的解的数量
        for (int i = 0; i < gpopsize; i++) {
            q = 0;
            for (int j = 0; j < gpopsize; j++) {
                ret = dominate(this.ind[i], this.ind[j]);
                if (ret == 1)
                    S[i]++;
                else if (ret == 2) {
                    bedomarr[i][q] = j;
                    q++;
                }
            }
            bedomnum[i] = q;
        }

        // 然后计算个体的原始适应度 R(i)
        for (int i = 0; i < gpopsize; i++) {
            for (int j = 0; j < bedomnum[i]; j++) {
                ind = bedomarr[i][j];
                R[i] += S[ind];
            }
        }

        // 接着计算个体的密度值 D(i) 和适应度 F(i)
        for (int i = 0; i < gpopsize; i++) {
            for (int j = i; j < gpopsize; j++) {
                if (i == j)
                    globdis[i][j] = GlobData.inf;
                else {
                    globdis[i][j] = 0.0;
                    for (int l = 0; l < GlobData.nfunc; l++) {
                        globdis[i][j] += Math.pow((this.ind[i].fvalue[l] - this.ind[j].fvalue[l]), 2.0);
                    }
                    globdis[i][j] = Math.pow(globdis[i][j], 0.5);
                    globdis[j][i] = globdis[i][j];
                }
            }
        }

        // 距离数组复制给存放距离数组archdis
        for (int i = 0; i < gpopsize; i++) {
            for (int j = 0; j < gpopsize; j++) {
                GlobData.archdis[i][j] = globdis[i][j];
            }
//            System.arraycopy(globdis[i], 0, GlobData.archdis[i], 0);
        }

        // 然后对距离数组进行排序，获取第 k 近邻个体的距离Sigma，再计算 D(i)、F(i)
        k = (int) Math.sqrt(gpopsize);
        for (int i = 0; i < gpopsize; i++) {
            Arrays.sort(globdis[i]);
            sigma = globdis[i][k - 1];
            D[i] = 1.0 / (sigma + 2.0);
            fitness[0][i] = i;
            fitness[1][i] = R[i] + D[i];
        }

        for (int i = 0; i < gpopsize; i++) {
            this.ind[i].fitness = fitness[1][i];
        }
    }

    // 对组合种群的适应度连同个体编号进行升序排序
    public void globFitSort() {
        double tmp1, tmp2;
        boolean flag;

        for (int i = 0; i < 2 * GlobData.popsize - 1; i++) {
            flag = true;
            for (int j = i + 1; j < 2 * GlobData.popsize; j++) {
                if (fitness[1][j] < fitness[1][i]) {
                    tmp1 = fitness[1][j];
                    fitness[1][j] = fitness[1][i];
                    fitness[1][i] = tmp1;

                    tmp2 = fitness[0][j];
                    fitness[0][j] = fitness[0][i];
                    fitness[0][i] = tmp2;

                    flag = false;
                }
            }
            if (flag)
                break;
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
