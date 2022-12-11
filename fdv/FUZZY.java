package fdv;

// 注意事项：多目标问题的测试含糊是在Fitness文件中，对于该算法的参数设置外面可以在main函数以及Initial函数中进行修改
// 模糊算法是在NSGA-II算法上实现的，但是模糊算法多了一个模糊操作，即在生成个体时使用了模糊函数（功能：类似于对一个数进行四舍五入）
// 输入：多目标问题的测试函数
// 输出：优化后的多目标问题的目标值
// 效果：得到的结果跟pf面接近


import javax.swing.*;
import java.io.*;
import java.util.Random;

import mytools.MyPaint;

public class FUZZY {
    int Gens;      // 迭代次数
    int Pop_Size;  // 种群大小
    int Dim_Y;     //目标维度
    int Dim_X;     // 决策变量维度
    double pc ;    // 交叉概率 pc 设置为 1.0
    double pm ;    // 变异概率 pm 设置为 1/D
    int nc ;       // 交叉的分布指数设置为 nc=20；
    int nm ;       // 突变的分布指数设置为 nm=20；
    double Rate ;  // 模糊进化比率
    double iter;   // 进化比率
    double[] MIN_X;// 存放决策变量的下界
    double[] MAX_X;// 存放决策变量的上界
    String testFun;// 测试函数名称

    FUZZY.Particle Population[];
    FUZZY.Particle Offspring[];
    FUZZY.Particle SumPop[];

    Random rand = new Random(System.currentTimeMillis());
    Fitness ft = new Fitness();
    Nondoninatedsort ND;

    class Particle {
        double[] position;  // 决策变量
        double[] fitness;   // 目标值
        public Particle(int dimX, int dimY) {
            position = new double[dimX];
            fitness = new double[dimY];
        }

        public void setPosition(double[] position) {
            this.position = position;
        }

        public void setFitness(double[] fitness) {
            this.fitness = fitness;
        }

        public double[] getFitness() {
            return this.fitness;
        }

        public double[] getPosition() {
            return this.position;
        }
    }

    double CrowdDistance[];
    int FrontRank[];
    double Zmin[];
    double Zmax[];

    public FUZZY(String ttN, int Populaton_Size, int Generation_Num, int DimDecision_X, int DimObjective_Y, double[] MinValue_X, double[] MaxValue_X) {
        Initial(ttN, Populaton_Size, Generation_Num, DimDecision_X, DimObjective_Y, MinValue_X, MaxValue_X);
        Initialize();
    }

    //设置初始参数
    private void Initial(String tN, int Pop_Size, int Gens, int DimDecision_X, int DimObjective_Y, double[] MinValue_X, double[] MaxValue_X) {
        this.testFun = tN;
        this.Pop_Size = Pop_Size;
        this.Gens = Gens;
        this.Dim_X = DimDecision_X;
        this.Dim_Y = DimObjective_Y;
        this.pc    = 1.0;
        this.pm    = 1.0/Dim_X;
        this.nc    = 20;
        this.nm    = 20;
        this.Rate  = 0.8;
        this.MIN_X = MinValue_X;
        this.MAX_X = MaxValue_X;
    }

    // Initialize函数：进行种群的初始化
    private void Initialize() {
        // 声明父代种群、子代种群、父子代合并的种群
        SumPop = new FUZZY.Particle[Pop_Size * 2];
        Population = new FUZZY.Particle[Pop_Size];
        Offspring = new FUZZY.Particle[Pop_Size];
        CrowdDistance = new double[Pop_Size];
        FrontRank = new int[Pop_Size];
        // 随机生成初始种群
        for (int i = 0; i < Pop_Size; i++) {
            Population[i] = new FUZZY.Particle(Dim_X, Dim_Y);
            Offspring[i] = new FUZZY.Particle(Dim_X, Dim_Y);
            for (int j = 0; j < Dim_X; j++) {
                Population[i].position[j] = rand.nextDouble() * (MAX_X[j] - MIN_X[j]) + MIN_X[j];
            }
        }
        for (int i = 0; i < Pop_Size * 2; i++) {
            SumPop[i] = new FUZZY.Particle(Dim_X, Dim_Y);
        }
        // 非支配排序的外部类，初始化备用
        ND = new Nondoninatedsort(Pop_Size * 2);
        //初始化Zmin和Zmax，用来归一化
        Zmin = new double[Dim_Y];
        Zmax = new double[Dim_Y];
    }

    // Init函数：计算初始种群的适应值，非支配等级，拥挤度距离
    // 输入：种群的个体
    // 输出：适应度值，非支配等级，个体的拥挤度距离
    private void Init() {
        Population = CalFitness(Population, Pop_Size);
        int[] FrontNo = ND.NDsort(Population, Pop_Size, Dim_Y);
        double[] CrowdDis = CrowdingDistance(Population, Pop_Size, FrontNo);
        for (int i = 0; i < Pop_Size; i++) {
            FrontRank[i] = FrontNo[i];
            CrowdDistance[i] = CrowdDis[i];
        }
    }

    // MatingPool函数：二元锦标赛选择，选择更优的个体进行遗传算子的操作
    // 输入：父代个体
    // 输出：经过锦标赛选择后的子代个体
    private void MatingPool() {
        for (int i = 0; i < Pop_Size; i++) {
            // 随机选择两个个体
            int parent1_index = new Random().nextInt(Pop_Size);
            int parent2_index = new Random().nextInt(Pop_Size);
            // 确保个体不同
            while (Is_equal(Population, parent1_index, parent2_index)) {
                parent2_index = new Random().nextInt(Pop_Size);
            }
            // 首先比较选择个体的PF等级
            if (FrontRank[parent1_index] < FrontRank[parent2_index]) {
                for (int j = 0; j < Dim_X; j++) {
                    Offspring[i].position[j] = Population[parent1_index].position[j];
                }
            }
            else if (FrontRank[parent1_index] > FrontRank[parent2_index]) {
                for (int j = 0; j < Dim_X; j++) {
                    Offspring[i].position[j] = Population[parent2_index].position[j];
                }
            }
            // 如果等级相同，随机
            else {
                if (CrowdDistance[parent1_index] < CrowdDistance[parent2_index]) {
                    for (int j = 0; j < Dim_X; j++) {
                        Offspring[i].position[j] = Population[parent2_index].position[j];
                    }
                } else {
                    for (int j = 0; j < Dim_X; j++) {
                        Offspring[i].position[j] = Population[parent1_index].position[j];
                    }
                }
            }
        }
    }

    // Is_equal函数：判断在二元锦标赛选择的个体是否是同一个个体
    // 输入：种群，选择的两个个体
    // 输出：逻辑值
    private boolean Is_equal(FUZZY.Particle[] Pop, int a, int b) {
        boolean logical = true;
        for (int i = 0; i < Dim_X; i++) {
            if (Population[a].position[i] != Population[b].position[i]) {
                logical = false;
                break;
            }
        }
        return logical;
    }

    // OperationGA函数：进行遗传算子操作，即交叉变异操作
    // 输入：父代个体
    // 输出：子代个体以及父子代合并的个体
    private void OperationGA() {
        // 模拟二元交叉
        for (int i = 0; i < Pop_Size; i += 2) {
            double rdc = rand.nextDouble();
            double beta1;
            if (rdc > 0.5)
                beta1= Math.pow(1 / (2 - rdc * 2), 1.0 / (nc + 1));
            else
                beta1= Math.pow(rdc * 2, 1.0 / (nc + 1));
            for (int j = 0; j < Dim_X; j++) {
                if (rand.nextDouble() < pc) {
                    double X1 = Offspring[i].position[j];
                    double X2 = Offspring[i + 1].position[j];
                    Offspring[i].position[j] = 0.5 * (1 + beta1) * X1 + (1 - beta1) * X2;
                    Offspring[i + 1].position[j] = 0.5 * (1 - beta1) * X1 + (1 + beta1) * X2;
                }
            }
        }
        // 多项式变异
        for (int i = 0; i < Pop_Size; i++) {
            for (int j = 0; j < Dim_X; j++) {
                if (rand.nextDouble() < pm) {
                    double rdm = rand.nextDouble();
                    double beta2;
                    if (rdm > 0.5)
                        beta2 = 1 - Math.pow(2 * (1 - rdm), 1.0 / (nm + 1));
                    else
                        beta2= Math.pow(rdm * 2, 1.0 / (nm + 1)) - 1;
                    Offspring[i].position[j] = Offspring[i].position[j] + beta2 * (MAX_X[j] - MIN_X[j]);
                }
                // 边界控制，防止前面的操作使得X超出范围
                if (Offspring[i].position[j] > MAX_X[j])
                    Offspring[i].position[j] = MAX_X[j];
                else if (Offspring[i].position[j] < MIN_X[j])
                    Offspring[i].position[j] = MIN_X[j];
            }
        }
        // 合并父子代个体
        for (int i = 0; i < Pop_Size; i++) {
            for (int j = 0; j < Dim_X; j++) {
                SumPop[i + Pop_Size].position[j] = Offspring[i].position[j];
            }
        }
        for (int i = 0; i < Pop_Size; i++) {
            for (int j = 0; j < Dim_X; j++) {
                SumPop[i].position[j] = Population[i].position[j];
            }
        }
    }

    // Fdv函数：模糊化操作
    // 输入：子代个体
    // 输出：模糊化后的子代个体及父子代合并的种群
    private void Fdv(double Iter) {              //论文算法核心思想
        int Total = 1;
        double  Acc = 0.4;
        int S;
        double a = Math.floor(Math.sqrt(2 * Rate * Total / Acc));
        S = (int) a;
        double[][] Step = new double[1][S + 2];
        for (int j = 0; j < S + 2; j++) {
            Step[0][j] = 0;
        }
        for (int i = 1; i <= S; i++) {
            Step[0][i] = (S * i * 1.0 - 1.0 * Math.pow(i, 2) / 2.0) * Acc;
        }
        Step[0][S + 1] = Rate * Total;
        // 模糊操作
        double R = 1.0;
        double[] gamma_a = new double[Dim_X];
        double[] gamma_b = new double[Dim_X];
        double[] miu1 = new double[Dim_X];
        double[] miu2 = new double[Dim_X];
        // 对每个个体的所有决策变量进行模糊化操作
        for (int i = 0; i < S + 1; i++) {
            if (Iter > Step[0][i] && Iter <= Step[0][i + 1]) {
                for (int j = 0; j < Pop_Size; j++) {
                    for (int k = 0; k < Dim_X; k++) {
                        gamma_a[k] = R * Math.pow(10, -(i + 1)) * Math.floor(Math.pow(10, i + 1) * Math.pow(R, -1) * (Offspring[j].position[k] - MIN_X[0])) + MIN_X[0];
                        gamma_b[k] = R * Math.pow(10, -(i + 1)) * Math.ceil(Math.pow(10, i + 1) * Math.pow(R, -1) * (Offspring[j].position[k] - MIN_X[0])) + MIN_X[0];
                        miu1[k] = 1.0 / (Offspring[j].position[k] - gamma_a[k]);
                        miu2[k] = 1.0 / (gamma_b[k] - Offspring[j].position[k]);
                        if (miu1[k] > miu2[k]) {
                            Offspring[j].position[k] = gamma_a[k];
                        }
                        else {
                            if (miu1[k] == miu2[k]) {
                                Random r = new Random();
                                if (gamma_a[k] < gamma_b[k]) {
                                    Offspring[j].position[k] = (gamma_b[k] - gamma_a[k]) * r.nextDouble() + gamma_a[k];
                                }
                                else if (gamma_a[k] == gamma_b[k]) {
                                    Offspring[j].position[k] = gamma_a[k];
                                }
                                else {
                                    Offspring[j].position[k] = (gamma_a[k] - gamma_b[k]) * r.nextDouble() + gamma_b[k];
                                }
                            }
                            else {
                                Offspring[j].position[k] = gamma_b[k];
                            }
                        }

                    }
                }
            }
        }
        // 修正决策变量
        for (int i = 0; i < Pop_Size; i++) {
            for (int j = 0; j < Dim_X; j++) {
                if (Offspring[i].position[j] > MAX_X[j])
                    Offspring[i].position[j] = MAX_X[j];
                else if (Offspring[i].position[j] < MIN_X[j])
                    Offspring[i].position[j] = MIN_X[j];
            }
        }
        // 合并父子代种群
        for (int i = 0; i < Pop_Size; i++) {
            for (int j = 0; j < Dim_X; j++) {
                SumPop[i + Pop_Size].position[j] = Offspring[i].position[j];

            }
        }
        for (int i = 0; i < Pop_Size; i++) {
            for (int j = 0; j < Dim_X; j++) {
                SumPop[i].position[j] = Population[i].position[j];
            }
        }

    }

    // NSGAII函数：对个体进行环境选择
    // 输入：父子代合并后的种群
    // 输出：经过非支配排序、拥挤度距离选择后的新的父代种群
    private void NSGAII() {
        // 使用NSGAII算法进行环境选择
        int PopSize = Pop_Size * 2;
        SumPop = CalFitness(SumPop, PopSize);
        int[] FrontNo = ND.NDsort(SumPop, PopSize, Dim_Y);
        double[] CrowdDis = CrowdingDistance(SumPop, PopSize, FrontNo);
        // 选择个体进入下一代
        int[] UniqLoc = Unique(PopSize);
        int[] Next = new int[PopSize];
        int[] Loc = new int[PopSize];
        int MaxFNo;
        int P1 = 0;
        int P2 = 0;
        MaxFNo = FrontNo[PopSize + 1];
        for (int i = 0; i < PopSize; i++) {
            if (FrontNo[i] < MaxFNo) {
                Next[P1++] = i;
            } else if (FrontNo[i] == MaxFNo) {
                Loc[P2++] = i;
            }
        }
        // 已经找到了P1个下一代个体，还需要 Pop_Size - P1个
        for (int i = P1; i < Pop_Size; i++) {
            int P = 0;
            double Min = 0.0;
            for (int j = 0; j < P2; j++) {
                if (Loc[j] == 0) continue;
                if (CrowdDis[Loc[j]] > Min) {
                    Min = CrowdDis[Loc[j]];
                    P = j;
                }
            }
            Next[i] = Loc[P];
            Loc[P] = 0;
        }
        // 生成新的父代种群，计算适应度值、拥挤度距离、非支配等级
        for (int i = 0; i < Pop_Size; i++) {
            for (int j = 0; j < Dim_X; j++) {
                Population[i].position[j] = SumPop[Next[i]].position[j];
            }
            for (int j = 0; j < Dim_Y; j++) {
                Population[i].fitness[j] = SumPop[Next[i]].fitness[j];
            }
            CrowdDistance[i] = CrowdDis[Next[i]];
            FrontRank[i] = FrontNo[Next[i]];
        }
    }

    // Unique函数：确保合并种群的个体没有重复的
    // 输入：种群大小
    // 输出：不重复个体的位置
    private int[] Unique(int PopSize) {
        int[] Location = new int[Pop_Size * 2 + 1];
        int p = 1;
        // 第一个解必不可能重复
        Location[0] = 0;
        for (int i = 1; i < PopSize; i++) {
            boolean flag = true;
            for (int j = i - 1; j >= 0; j--) {
                flag = true;
                for (int k = 0; k < Dim_X; k++) {
                    if (SumPop[i].position[k] != SumPop[j].position[k]) {
                        flag = false;
                        break;
                    }
                }
                // 如果在比较了一个解之后发现完全一致，则确定这个解被去掉，可以退出循环
                if (flag) break;
            }
            // 如果是独立解，则加入Location里面，不是则跳过
            if (!flag) {
                Location[p] = i;
                p++;
            }
        }
        Location[Pop_Size * 2] = p - 1;
        return Location;
    }

    // CalFitness函数：计算种群适应度值
    // 输入：种群，种群数量
    // 输出：含有适应度值的种群
    private FUZZY.Particle[] CalFitness(FUZZY.Particle[] Pop, int Num) {

        for (int i = 0; i < Num; i++) {
            switch (this.testFun){
                case "ZDT1":
                    double[] fitnessZDT1 = ft.CalFitnessZDT1(Pop[i].position);
                    for (int j = 0; j < Dim_Y; j++) {
                        Pop[i].fitness[j] = fitnessZDT1[j];
                    }
                    break;
                case "ZDT2":
                    double[] fitnessZDT2 = ft.CalFitnessZDT2(Pop[i].position);
                    for (int j = 0; j < Dim_Y; j++) {
                        Pop[i].fitness[j] = fitnessZDT2[j];
                    }
                    break;
                case "ZDT3":
                    double[] fitnessZDT3 = ft.CalFitnessZDT3(Pop[i].position);
                    for (int j = 0; j < Dim_Y; j++) {
                        Pop[i].fitness[j] = fitnessZDT3[j];
                    }
                    break;
                case "DTLZ1":
                    double[] fitnessDTLZ1 = ft.CalFitnessDTLZ1(Pop[i].position, Dim_Y);
                    for (int j = 0; j < Dim_Y; j++) {
                        Pop[i].fitness[j] = fitnessDTLZ1[j];
                    }
                    break;
                case "DTLZ3":
                    double[] fitnessDTLZ3 = ft.CalFitnessDTLZ3(Pop[i].position, Dim_Y);
                    for (int j = 0; j < Dim_Y; j++) {
                        Pop[i].fitness[j] = fitnessDTLZ3[j];
                    }
                    break;
                case "DTLZ7":
                    double[] fitnessDTLZ7 = ft.CalFitnessDTLZ7(Pop[i].position, Dim_Y);
                    for (int j = 0; j < Dim_Y; j++) {
                        Pop[i].fitness[j] = fitnessDTLZ7[j];
                    }
                    break;
            }
        }
        return Pop;
    }

    // CrowdingDistance函数：计算个体的拥挤度距离
    // 输入：种群，种群大小，排序好的个体
    // 输出：个体的拥挤度距离
    private double[] CrowdingDistance(FUZZY.Particle[] Pop, int PopSize, int[] FrontNo) {
        double[] CrowdDis = new double[PopSize];
        double INF = 10000.0;
        int MaxFNo = FrontNo[PopSize];
        for (int i = 1; i <= MaxFNo; i++) {
            int[] Front = find(FrontNo, i, PopSize);
            int Front_num = Front[0];
            if (Front_num == 0 )
                continue;
            for (int j = 0; j < Dim_Y; j++) {
                int[] Rank = sortFront(Pop, Front, j);
                // 对排在第一个和最后有个个体的拥挤度设为无穷大，保证多样性
                CrowdDis[Front[Rank[0]]] = INF;
                CrowdDis[Front[Rank[Front_num - 1]]] = INF;
                if (Front_num == 1 || Front_num == 2) continue;
                // 计算拥挤度
                for (int k = 1; k < Front_num - 1; k++) {
                    CrowdDis[ Front[ Rank[k] ] ] = CrowdDis[ Front[ Rank[k] ] ]
                            + (Pop[ Front[ Rank[k + 1] ] ].fitness[j] - Pop[ Front[ Rank[k - 1] ] ].fitness[j] )
                            /(Zmax[j] - Zmin[j]);
                }
            }
        }
        return CrowdDis;
    }

    // 找到该非支配排序等级的解的位置
    private int[] find(int[] FrontNo, int point, int PopSize) {
        int num = 0;//数组的第一位记录个数
        int P = 1;
        int[] Front = new int[PopSize];
        for (int i = 0; i < PopSize; i++) {
            if (FrontNo[i] == point) {
                Front[P] = i;
                num++;
                P++;
            }
        }
        Front[0] = num;
        return Front;
    }

    // sortFront函数：对种群中的每个个体按照第一维的适应度进行排序（从小到大）
    // 输入：种群
    // 输出：每个个体按照第一维的适应度排好序的二维数组
    private int[] sortFront(FUZZY.Particle[] Pop, int[] Front, int Dim) {
        int Front_num = Front[0];//首先取出该层的解的这一维的坐标
        int[] Rank = new int[Front_num];
        double[] Obj = new double[Front_num];
        double INF = 100000.0;
        for (int i = 0; i < Front_num; i++) {
            int P = Front[i + 1];
            Obj[i] = Pop[P].fitness[Dim];
            Rank[i] = 0;
        }
        for (int i = 0; i < Front_num; i++) {
            double min = INF;
            int j;
            for (j = 0; j < Front_num; j++) {
                if (Obj[j] < min) {
                    Rank[i] = j;
                    min = Obj[j];
                }
            }
            Obj[Rank[i]] = INF;
        }
        for (int i = 0; i < Front_num; i++) {
            Rank[i] += 1;
        }
        Zmin[Dim] = Pop[Front[Rank[0]]].fitness[Dim];
        Zmax[Dim] = Pop[Front[Rank[Front_num - 1]]].fitness[Dim];
        return Rank;
    }

    // A_Fuzzy_Decision_Variables_Framework_for_Large-scale_Multiobjective_Optimization的框架
    public void run() {
        Init();                        // 对种群进行初始化
        for (int i = 1; i <= Gens; i++) {
            MatingPool();             // 交配池（二元锦标赛选择）
            OperationGA();           // 遗传算子操作（交叉变异）
            iter = i * 1.0 / Gens;   // 模糊进化比率
            if (iter <= Rate) {
                Fdv(iter);           // 模糊操作
            }
            NSGAII();  // NSGA-II算法（环境选择，选择更优的个体生成新的父代种群）
        }
        // 输出结果
        String dataStr = "";
        for (int i = 0; i < Pop_Size; i++) {
            int lengthRec = dataStr.length();
            if(FrontRank[i] > 1){
                continue;
            }
            for (int j = 0; j < Dim_Y; j++) {
                System.out.print(Population[i].fitness[j] + "\t");
                dataStr = dataStr + Population[i].fitness[j] + "\t";

            }
            dataStr +="\n";
            System.out.println();
        }

        try {
            File file = new File("E:\\algorithmTest\\res\\fdvPF_"+ testFun +"obj"+Dim_Y + "var"+Dim_X +".txt");
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


    public double calIGD(double[][] pf){

        double sumD = 0;
        int samplePointNum = pf.length;
        for (int i = 0 ; i < samplePointNum; i ++){
            double minDis = Double.MAX_VALUE;
            for(int j = 0 ; j < Pop_Size; j ++){
                double dis = 0;
                for (int k = 0; k < Dim_Y; k++){
                    dis += Math.pow(Population[j].fitness[k] - pf[i][k],2);
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

    public void paint2Obj(){
        JFrame f = new JFrame();
        double[] myData = new double[this.Pop_Size * 2];
        for (int i = 0; i < Pop_Size; i = i + 2){
            myData[i] = Population[i].fitness[0];
            myData[i + 1] = Population[i].fitness[1];
        }
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.add(new MyPaint(myData));
        f.setSize(300,300);
        f.setLocation(200,200);
        f.setVisible(true);
    }

}

