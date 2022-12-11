package spea2;

import java.io.*;
import java.util.Random;

public class MyAlgorithm {
    int fun_input;

    // 种群声明
    Population oldpop;      // 旧一代种群
    Population newpop;      // 新一代种群
    Population matepop;     // 交配池
    GlobalPop globpop;      // 组合种群
    GlobalPop archive;      // 归档集

    public MyAlgorithm(int fun_input,int popSize, int Gen, int varNum, int objNum) throws IOException {
        // 首先进行参数设置
        this.fun_input = fun_input;
        GlobData.setPara(fun_input,popSize,Gen,varNum,objNum);        // 1 / 2

        initialize();       // 初始化
        run();              // 种群迭代

        // 保存最佳一代种群
        String dataStr = "";
        for (int i = 0; i < GlobData.popsize; i++) {
            for (int j = 0; j < GlobData.nfunc; j++) {
                System.out.print(String.format("%.6f", matepop.ind[i].fvalue[j]) + "\t");
                dataStr += matepop.ind[i].fvalue[j] + "\t";
            }
            dataStr += "\n";
            System.out.println();
        }
        try {
            File file = new File("E:\\algorithmTest\\res\\spea2PF_"+ fun_input +"obj"+GlobData.nfunc + "var" + GlobData.nvar+".txt");
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
    }

    // 种群初始化
    private void initialize() {
        oldpop = new Population();
        newpop = new Population();
        matepop = new Population();
        globpop = new GlobalPop();
        archive = new GlobalPop();

        oldpop.initPop();           // 初始化种群
        oldpop.func(fun_input);     // 计算初始种群的函数值
        oldpop.calFitness();        // 计算初始种群的适应度

    }

    // 种群迭代
    private void run() {
        for (int gen = 0; gen < GlobData.maxgen; gen++) {
            //System.out.println("Generation = " + (gen + 1));

            tSelect(oldpop, matepop);       // 二进制锦标赛选择
            binCross(newpop, matepop);      // 模拟二进制交叉
            mutate(newpop);                 // 多项式变异

            newpop.func(fun_input);         // 计算新一代种群的函数值
            newpop.calFitness();            // 计算新一代种群的适应度
            enviSelect(oldpop, newpop, matepop);    // 环境选择

            // 新旧种群交替
            for (int i = 0; i < GlobData.popsize; i++) {
                // 变量复制
                System.arraycopy(matepop.ind[i].xreal, 0, oldpop.ind[i].xreal, 0, GlobData.nvar);
                // 函数值复制
                System.arraycopy(matepop.ind[i].fvalue, 0, oldpop.ind[i].fvalue, 0, GlobData.nfunc);
                // 适应度复制
                oldpop.ind[i].fitness = matepop.ind[i].fitness;
            }
        }

    }

    // 锦标赛选择
    private void tSelect(Population oldpop, Population matepop) {
        int rnd1, rnd2;
        Individual par1, par2;
        Random random = new Random();

        // 选择 popsize 个个体加入交配池
        for (int i = 0; i < GlobData.popsize; i++) {
            rnd1 = random.nextInt(GlobData.popsize) + 1;    // [1, popsize]
            rnd2 = random.nextInt(GlobData.popsize) + 1;

            if (rnd1 == GlobData.popsize)
                rnd1 = (GlobData.popsize - 2) / 2;
            if (rnd2 == GlobData.popsize)
                rnd2 = (GlobData.popsize - 4) / 2;

            par1 = oldpop.ind[rnd1 - 1];
            par2 = oldpop.ind[rnd2 - 1];

            if (par1.fitness < par2.fitness) {
                System.arraycopy(par1.xreal, 0, matepop.ind[i].xreal, 0, GlobData.nvar);
            } else {
                System.arraycopy(par2.xreal, 0, matepop.ind[i].xreal, 0, GlobData.nvar);
            }
        }
    }

    // 模拟二进制交叉
    private void binCross(Population newpop, Population matepop) {
        int k1 = 0, k2 = 0;     // k1 定位旧种群下标，k2 定位新种群下标
        double p1, p2, beta, betaq, alpha, u, expp;
        double x1, x2, child1, child2, xl, xu;

        for (int i = 0; i < GlobData.popsize / 2; i++) {
            if (Math.random() <= GlobData.pcross) {      // 若个体发生交叉
                for (int j = 0; j < GlobData.nvar; j++) {
                    p1 = matepop.ind[k1].xreal[j];
                    p2 = matepop.ind[k1 + 1].xreal[j];

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
                            expp = GlobData.di + 1;
                            alpha = 2 - Math.pow(beta, -expp);

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
                    newpop.ind[k2].xreal[j] = child1;
                    newpop.ind[k2 + 1].xreal[j] = child2;
                }
            } else {        // 若个体不发生交叉
                for (int j = 0; j < GlobData.nvar; j++) {
                    child1 = matepop.ind[k1].xreal[j];
                    child2 = matepop.ind[k1 + 1].xreal[j];

                    newpop.ind[k2].xreal[j] = child1;
                    newpop.ind[k2 + 1].xreal[j] = child2;
                }
            }
            k1 += 2;
            k2 += 2;
        }
    }

    // 个体多项式变异
    private void mutate(Population newpop) {
        double x, xl, xu, delta, deltaq, u, max, expp, val1, val2;

        for (int i = 0; i < GlobData.popsize; i++) {        // 种群所有个体参加变异
            for (int j = 0; j < GlobData.nvar; j++) {
                if (Math.random() <= GlobData.pmut) {    // 若发生变异
                    x = newpop.ind[i].xreal[j];
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

                        newpop.ind[i].xreal[j] = x;

                    } else {        // 变量越界
                        val1 = Math.random();
                        newpop.ind[i].xreal[j] = val1 * (xu - xl) + xl;
                    }
                }
            }
        }
    }

    // 环境选择
    private void enviSelect(Population oldpop, Population newpop, Population matepop) {
        int cnt;        // 当前加入到存档种群中的个体数
        int gpopsize;
        int arcind;     // 追踪加入到存档的个体下标
        int row = 0, col = 0, minind, ind1, ind2;
        double min;

        // 合并 oldpop 和 newpop，产生全部的交配池，初始化存档的种群标志
        gpopsize = 2 * GlobData.popsize;
        for (int i = 0; i < GlobData.popsize; i++) {
            // 旧种群变量 + 新种群变量
            for (int j = 0; j < GlobData.nvar; j++) {
                globpop.ind[i].xreal[j] = oldpop.ind[i].xreal[j];
                globpop.ind[i + GlobData.popsize].xreal[j] = newpop.ind[i].xreal[j];
            }

            // 旧种群函数值 + 新种群函数值
            for (int j = 0; j < GlobData.nfunc; j++) {
                globpop.ind[i].fvalue[j] = oldpop.ind[i].fvalue[j];
                globpop.ind[i + GlobData.popsize].fvalue[j] = newpop.ind[i].fvalue[j];
            }

            // 初始化存档的个体标志
            archive.flag[i] = 0;
            archive.flag[i + GlobData.popsize] = 0;
        }

        // 重计算组合种群的适应度
        globpop.globFitness();

        // 初始化归档的个体标志
        for (int i = 0; i < gpopsize; i++) {
            archive.flag[i] = 0;    // 标志为 0 表示还没有任何个体加入到存档中
        }

        // 第一步：先将组合种群中所有适应度都低于 1 的个体加入到下一代存档种群中
        cnt = 0;
        for (int i = 0; i < gpopsize; i++) {
            if (globpop.fitness[1][i] < 1) {
                cnt++;
                archive.flag[i] = 1;    // 标志个体 i 加入了下一代存档
            }
        }

        // 第二步：若下一代存档种群未满，则将前 popsize-cnt 个F(i)>=1的个体加入到下一代存档种群中
        if (cnt < GlobData.popsize) {
            globpop.globFitSort();      // 连同个体编号一同排序，排序后获得前 popsize-cnt 个 F(i) >= 1 的个体
            for (int i = 0; i < gpopsize; i++) {
                if (cnt < GlobData.popsize && globpop.fitness[1][i] >= 1) {
                    cnt++;
                    arcind = (int) globpop.fitness[0][i];
                    archive.flag[arcind] = 1;
                }
                if (cnt >= GlobData.popsize)
                    break;
            }
        } else if (cnt > GlobData.popsize) {    // 第三步：若下一代存档种群溢出，则采用截断操作，迭代地从archive中删除个体
            // 首先找到最小的min，位于 (row, col) 的位置
            for (int k = 0; k < cnt - GlobData.popsize; k++) {  // 需要移除cnt-popsize个个体
                min = Integer.MAX_VALUE;
                for (int i = 0; i < gpopsize; i++) {
                    for (int j = 0; j < gpopsize; j++) {
                        // 要保证两个个体 i 和 j 都在归档种群中
                        if (archive.flag[i] == 1 && archive.flag[j] == 1 && GlobData.archdis[i][j] < min) {
                            min = GlobData.dis[i][j];
                            row = i;        // // 组合种群中最小距离的两个个体(row, col)
                            col = j;
                        }
                    }
                }
                // 将个体 row 和 col 的距离复制到辅助数组，连同下标一起复制
                for (int i = 0; i < gpopsize; i++) {
                    // 第 row 个个体
                    GlobData.disrow[0][i] = i;                          // 第一行存放下标
                    GlobData.disrow[1][i] = GlobData.archdis[row][i];   // 第二行存放下标对应个体的距离

                    // 第 col 个个体
                    GlobData.discol[0][i] = i;
                    GlobData.discol[1][i] = GlobData.archdis[col][i];
                }

                uniSort();      // 对两个辅助数组连同下标进行升序排序

                // 依次对比disrow 和 discol 上的第二行，将最先出现更小值的个体用minind 保存
                minind = row;
                for (int i = 0; i < gpopsize; i++) {
                    ind1 = (int) GlobData.disrow[0][i];
                    ind2 = (int) GlobData.discol[0][i];

                    if (archive.flag[ind1] == 1 && archive.flag[ind2] == 1 && GlobData.disrow[1][i] != GlobData.discol[1][i]) {
                        if (GlobData.disrow[1][i] < GlobData.discol[1][i])
                            minind = row;
                        else
                            minind = col;
                        break;
                    }
                }
                archive.flag[minind] = 0;       // 标志个体minind 移出存档
            }
        }

        // 到此就完成了环境选择操作，将标记为 1 的个体加入到matepop中
        cnt = 0;
        for (int i = 0; i < gpopsize; i++) {
            if (archive.flag[i] == 1) {
                matepop.ind[cnt++] = globpop.ind[i];
            }
        }
    }

    // 对两个辅助数组 disrow，discol 连同下标一起升序排序
    private void uniSort() {
        boolean flag;
        double tmp1, tmp2;

        for (int i = 1; i < 2 * GlobData.popsize; i++) {
            flag = true;
            for (int j = 0; j < 2 * GlobData.popsize - i; j++) {
                if (GlobData.disrow[1][j + 1] < GlobData.disrow[1][j]) {
                    tmp1 = GlobData.disrow[1][j + 1];
                    GlobData.disrow[1][j + 1] = GlobData.disrow[1][j];
                    GlobData.disrow[1][j] = tmp1;

                    tmp2 = GlobData.disrow[0][j + 1];
                    GlobData.disrow[0][j + 1] = GlobData.disrow[0][j];
                    GlobData.disrow[0][j] = tmp2;

                    flag = false;
                }
            }
            if (flag)
                break;
        }
    }

    public double calIGD(double[][] pf){
        double sumD = 0;
        int samplePointNum = pf.length;
        for (int i = 0 ; i < samplePointNum; i ++){
            double minDis = Double.MAX_VALUE;
            for(int j = 0 ; j < GlobData.popsize; j ++){
                double dis = 0;
                for (int k = 0; k < GlobData.nfunc; k++){
                    dis += Math.pow(matepop.ind[j].fvalue[k] - pf[i][k],2);
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
