package moead;


import function.MyFunction;

import java.io.*;
import java.util.Random;

// 算法类
public class MyAlgorithm{
    int fun_input;
    FileWriter fw = new FileWriter("fit.txt");
    BufferedWriter fp = new BufferedWriter(fw);

    // 种群声明
    Population pop;      // 种群
    Individual par1;
    Individual par2;
    Individual besty;

    // 构造方法
    public MyAlgorithm(MyFunction func, int fun_input, int popSize, int Gen, int varNum, int objNum) throws IOException {
        // 首先进行参数设置
        this.fun_input = fun_input;
        GlobData.setPara(fun_input,popSize,Gen,varNum,objNum);

        initialize();       // 初始化
        run();              // 种群迭代

        // 最佳一代种群
        for (int i = 0; i < GlobData.popsize; i++) {
            for (int j = 0; j < GlobData.nfunc; j++) {
                fp.write(String.format("%.6f", pop.ind[i].fitness[j]) + "\t");
            }
            fp.newLine();
        }

        // 打印最佳一代种群
        String dataStr = "";
        for (int i = 0; i < moead.GlobData.popsize; i++) {
            for (int j = 0; j < moead.GlobData.nfunc; j++) {
                System.out.print(String.format("%.6f", pop.ind[i].fitness[j]) + "\t");
                dataStr += pop.ind[i].fitness[j] + "\t";
            }
            dataStr += "\n";
            System.out.println();
        }
        try {
            File file = new File("E:\\algorithmTest\\res\\moeadPF_"+ fun_input +"obj"+ moead.GlobData.nfunc + "var" + moead.GlobData.nvar+".txt");
            if(!file.exists()) {
                file.createNewFile();
            }
            FileOutputStream fos = new FileOutputStream(file,true);
            OutputStreamWriter osw = new OutputStreamWriter(fos);
            BufferedWriter bw = new BufferedWriter(osw);
            bw.write(dataStr);
            bw.newLine();
            bw.flush();
            bw.close();
            osw.close();
            fos.close();
        }catch (FileNotFoundException e1) {
            e1.printStackTrace();
        } catch (IOException e2) {
            e2.printStackTrace();
        }
        fp.close();
    }

    // 初始化
    private void initialize() {
        pop = new Population();
        par1 = new Individual();
        par2 = new Individual();
        besty = new Individual();

        pop.initPop();
        pop.func(fun_input);

        initLamda();            // 产生N个均匀分布的权向量
        calculateDis();         // 计算任意两个权向量之间的距离，得到每个权向量最近的T个权向量矩阵
        pop.calBestZ();         // 计算当前种群的最优P值 Z
    }

    // 产生 N 个均匀分布的权向量
    private void initLamda() {
        for (int i = 0; i < GlobData.popsize; i++) {
            if (GlobData.nfunc == 2) {      // 目标函数个数为 2，因此权重系数维度为 2
                GlobData.lamda[i][0] = ((double) i) / GlobData.popsize;
                GlobData.lamda[i][1] = 1.0 - ((double) i) / GlobData.popsize;
            }
        }
    }

    // 计算任意两个权向量之间的距离，得到每个权向量最近的T个权向量矩阵
    private void calculateDis() {
        double d;       // 暂存距离

        // 初始化下标矩阵
        for (int i = 0; i < GlobData.popsize; i++) {
            GlobData.disindex[i] = i;
        }

        for (int i = 0; i < GlobData.popsize; i++) {
            for (int j = 0; j < GlobData.popsize; j++) {
                d = 0.0;
                for (int k = 0; k < GlobData.nfunc; k++) {      // 海明公式
                    d += Math.pow((GlobData.lamda[i][k] - GlobData.lamda[j][k]), 2);
                }
                d = Math.pow(d, 0.5);
                GlobData.dis[j] = d;
            }

            // 到这里就获得了第 i 个权向量到其余N个权向量之间的欧氏距离
            disSort();
            System.arraycopy(GlobData.disindex, 0, GlobData.B[i], 0, GlobData.T);
        }
    }

    // 距离数组升序排序
    private void disSort() {
        double tmpd;
        int tmpi;
        boolean flag;

        for (int i = 1; i < GlobData.popsize; i++) {
            flag = true;
            for (int j = 0; j < GlobData.popsize - i; j++) {
                if (GlobData.dis[j + 1] < GlobData.dis[j]) {
                    tmpd = GlobData.dis[j + 1];
                    GlobData.dis[j + 1] = GlobData.dis[j];
                    GlobData.dis[j] = tmpd;

                    tmpi = GlobData.disindex[j + 1];
                    GlobData.disindex[j + 1] = GlobData.disindex[j];
                    GlobData.disindex[j] = tmpi;

                    flag = false;
                }
            }
            if (flag)
                break;
        }
    }

    // 种群迭代
    private void run() {
        double gte_xi, gte_y;

        for (int gen = 0; gen < GlobData.maxgen; gen++) {
            //System.out.println("Generation = " + (gen + 1));

            for (int j = 0; j < GlobData.popsize; j++) {
                select();           // 选择两个个体 par1 和 par2
                binCross();         // 模拟二进制交叉 生成两个新个体 par1 par2
                mutate(par1);       // 新个体 par1 进行多项式变异
                mutate(par2);

                // 判断个体间的支配关系，获取优秀子代 besty
                besty = (dominate(par1, par2) == 1) ? par1 : par2;
//                besty.indFunc(fun_input);

                // 更新 Z
                for (int k = 0; k < GlobData.nfunc; k++) {
                    if (besty.fitness[k] < GlobData.Z[k])
                        GlobData.Z[k] = besty.fitness[k];
                }

                // 更新邻域解，遍历每个个体的所有邻近权向量，并将其切比雪夫值与孩子个体 besty 作比较，看是否需要用新个体 besty 来取代原来的个体
                for (int k = 0; k < GlobData.T; k++) {
                    gte_xi = Tchebycheff(pop.ind[GlobData.B[j][k]], GlobData.B[j][k]);
                    gte_y = Tchebycheff(besty, GlobData.B[j][k]);
                    if (gte_y <= gte_xi){   // 若满足条件，则由更好的 besty 个体取代种群中原先的个体
                        // besty 取代种群个体
                        System.arraycopy(besty.xreal, 0, pop.ind[GlobData.B[j][k]].xreal, 0, GlobData.nvar);
                        System.arraycopy(besty.fitness, 0, pop.ind[GlobData.B[j][k]].fitness, 0, GlobData.nfunc);
                    }
                }
            }
        }

    }

    // 选择两个进行交叉的个体
    private void select() {
        int rnd1, rnd2;
        Random random = new Random();

        rnd1 = random.nextInt(GlobData.popsize) + 1;        // 选择下标
        rnd2 = random.nextInt(GlobData.popsize) + 1;

        if (rnd1 == GlobData.popsize)
            rnd1 = (GlobData.popsize - 2) / 2;
        if (rnd2 == GlobData.popsize)
            rnd2 = (GlobData.popsize - 4) / 2;

        // 个体 par1 复制
        System.arraycopy(pop.ind[rnd1 - 1].xreal, 0, par1.xreal, 0, GlobData.nvar);
        System.arraycopy(pop.ind[rnd1 - 1].fitness, 0, par1.fitness, 0, GlobData.nfunc);

        // 个体 par2 复制
        System.arraycopy(pop.ind[rnd2 - 1].xreal, 0, par2.xreal, 0, GlobData.nvar);
        System.arraycopy(pop.ind[rnd2 - 1].fitness, 0, par2.fitness, 0, GlobData.nfunc);

    }

    // 模拟二进制交叉
    private void binCross() {
        double p1, p2, beta, betaq, alpha, u, expp;
        double x1, x2, child1, child2, xl, xu;

        if (Math.random() <= GlobData.pcross) {      // 若个体发生交叉
            for (int j = 0; j < GlobData.nvar; j++) {   // 两个个体的所有变量发生交叉
                p1 = par1.xreal[j];
                p2 = par2.xreal[j];

                xl = GlobData.lim[j][0];
                xu = GlobData.lim[j][1];

                if (Math.random() <= 0.5) {
                    if (Math.abs(p1 - p2) > 1e-6) {  // 若 p1 != p2
                        if (p1 < p2) {
                            x1 = p1;
                            x2 = p2;
                        } else {
                            x2 = p1;
                            x1 = p2;
                        }

                        // 计算 beta
                        if ((x1 - xl) < (xu - x2))
                            beta = 1 + 2.0 * (x1 - xl) / (x2 - x1);
                        else
                            beta = 1 + 2.0 * (xu - x2) / (x2 - x1);

                        // 计算 alpha
                        expp = GlobData.di + 1.0;
                        alpha = 2.0 - Math.pow(beta, -expp);

                        // 计算 betaq
                        u = Math.random();
                        expp = 1.0 / (GlobData.di + 1.0);
                        if (u <= (1.0 / alpha))
                            betaq = Math.pow(alpha * u, expp);
                        else
                            betaq = Math.pow(1.0 / (2.0 - alpha * u), expp);

                        // 计算child
                        child1 = 0.5 * ((x1 + x2) - betaq * (x2 - x1));
                        child2 = 0.5 * ((x1 + x2) + betaq * (x2 - x1));

                    } else {            // 若 p1 == p2
                        betaq = 1.0;
                        x1 = p1;
                        x2 = p2;
                        child1 = 0.5 * ((x1 + x2) - betaq * (x2 - x1));
                        child2 = 0.5 * ((x1 + x2) + betaq * (x2 - x1));
                    }
                    if (child1 < xl) child1 = xl;            // 变量越界 即约束条件
                    if (child1 > xu) child1 = xu;
                    if (child2 < xl) child2 = xl;
                    if (child2 > xu) child2 = xu;

                } else {    // 若变量不发生交叉
                    child1 = p1;
                    child2 = p2;
                }
                par1.xreal[j] = child1;
                par2.xreal[j] = child2;
            }
        }
    }

    // 多项式变异
    private void mutate(Individual child) {
        double x, xl, xu, delta, deltaq, u, max, expp, val1, val2;

        for (int j = 0; j < GlobData.nvar; j++) {
            if (Math.random() <= GlobData.pmut) {    // 若发生变异
                x = child.xreal[j];
                xl = GlobData.lim[j][0];
                xu = GlobData.lim[j][1];

                if (x > xl) {
                    // 计算 delta
                    if ((x - xl) < (xu - x))
                        delta = (x - xl) / (xu - xl);
                    else
                        delta = (xu - x) / (xu - xl);

                    // 计算 deltaq
                    u = Math.random();
                    expp = 1.0 / (GlobData.dim + 1.0);
                    if (u <= 0.5) {
                        val1 = Math.pow((1.0 - delta), (1.0 + GlobData.dim));
                        val2 = 2 * u + (1 - 2 * u) * val1;
                        deltaq = Math.pow(val2, expp) - 1;
                    } else {
                        val1 = 2 * (u - 0.5) * Math.pow((1.0 - delta), (1.0 + GlobData.dim));
                        val2 = 2 * (1 - u) + val1;
                        deltaq = 1.0 - Math.pow(val2, expp);
                    }

                    // 计算变异之后的子代 x
                    max = xu - xl;
                    x = x + deltaq * max;

                    if (x < xl) x = xl;
                    if (x > xu) x = xu;

                    child.xreal[j] = x;

                } else {        // 变量越界
                    val1 = Math.random();
                    child.xreal[j] = val1 * (xu - xl) + xl;
                }
            }
        }
    }

    // 个体间的支配关系
    public int dominate(Individual y1, Individual y2) {
        int greater = 0, equal = 0, less = 0;

        for (int i = 0; i < GlobData.nfunc; i++) {
            if (y1.fitness[i] < y2.fitness[i])
                less++;
            else if (Math.abs(y1.fitness[i] - y2.fitness[i]) <= 1e-6)
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

    // 切比雪夫聚合函数，返回的是表达式的最大值
    private double Tchebycheff(Individual y, int index) {
        double[] temp = new double[GlobData.nfunc];
        double maxtemp = 0;

        for (int i = 0; i < GlobData.nfunc; i++) {
            temp[i] = GlobData.lamda[index][i] * Math.abs(y.fitness[i] - GlobData.Z[i]);
            if (i == 0)
                maxtemp = temp[0];
            if (temp[i] > maxtemp)
                maxtemp = temp[i];
        }

        return maxtemp;
    }
    public double calIGD(double[][] pf){
        double sumD = 0;
        int samplePointNum = pf.length;
        for (int i = 0 ; i < samplePointNum; i ++){
            double minDis = Double.MAX_VALUE;
            for(int j = 0; j < moead.GlobData.popsize; j ++){
                double dis = 0;
                for (int k = 0; k < moead.GlobData.nfunc; k++){
                    dis += Math.pow(pop.ind[j].fitness[k] - pf[i][k],2);
                }
                dis = Math.sqrt(dis);
                if (dis < minDis){
                    minDis = dis;
                }
            }
            sumD += minDis;
        }
        return sumD/samplePointNum;
    }

}
