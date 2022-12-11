package nsgaii;

import java.io.*;
import java.util.Random;
import fdv.Fitness;

public class NSGAII {

    int gens;
    int Popsize;
    double pc ;
    double pm ;
    int dis_c;
    int dis_m;
    int dim_x;
    int dim_y;
    double[] MIN_X;
    double[] MAX_X;
    String testFun;// 测试函数名称
	
    class Particle {
        double[] position;
        double[] fitness;
        public Particle(int dimX,int dimY) {
            position = new double[dimX];
            fitness  = new double[dimY];
        }
    }
	
	Random rand = new Random(System.currentTimeMillis());
    Fitness ft = new Fitness();
	
    Particle[] Population;
    Particle[] Offspring;
    Particle[] Sumpop;
    double[] CrowdDistance;
    int[] FrontRank;
    int[] Loc;
    double[] Zmin;
    double[] Zmax;

    public NSGAII(String testFun, int populaton_size,int generation_num,int dim_decision,int dim_objective,double[] Min_Xvalue,double[] Max_Xvalue){
        this.testFun = testFun;
        set_mainpara(populaton_size,generation_num,dim_decision,dim_objective,Min_Xvalue,Max_Xvalue);
        initialise();
    }

    //设置算法的参数
    private void set_mainpara(int Popsize, int gens, int dim_decision, int dim_objective, double[] Min_Xvalue, double[] Max_Xvalue){
        //设置活动参数
        this.Popsize = Popsize;
        this.gens = gens;
        this.dim_x = dim_decision;
        this.dim_y = dim_objective;
        this.pc    = 0.9;
        this.pm    = 0.05;
        this.dis_c = 20;
        this.dis_m = 20;

        this.MIN_X = Min_Xvalue;
        this.MAX_X = Max_Xvalue;
    }

    //initialise函数；初始化需要参数设置的内容
    private void initialise() {
        //将三个种群申明备用
        //主种群的特殊属性，PF等级以及拥挤度距离，用来进行个体选择
        Sumpop = new Particle[Popsize * 2];
        Population    = new Particle[Popsize];
        Offspring     = new Particle[Popsize];
        CrowdDistance = new double[Popsize];
        FrontRank     = new int[Popsize];
        //非支配排序的记录数组
        Loc = new int[Popsize * 2];
        //种群的正式初始化
        for (int i = 0; i < Popsize; i++) {
            Population[i] = new Particle(dim_x, dim_y);
            Offspring[i] = new Particle(dim_x, dim_y);
            for (int j = 0; j < dim_x; j++) {
                Population[i].position[j] = rand.nextDouble() * (MAX_X[j] - MIN_X[j]) + MIN_X[j];
            }
        }
        for (int i = 0; i < Popsize * 2; i++) {
            Sumpop[i] = new Particle(dim_x, dim_y);
        }
        //初始化Zmin和Zmax，用来归一化
        Zmin = new double[dim_y];
        Zmax = new double[dim_y];
    }
	
	//run函数：调用该函数，进行算法的迭代
    public void run(){
        //NSGAII的总体框架，选择，遗传算子，环境选择
        init();
        for (int i = 0; i < gens; i++) {
            mating_selection();
            operator_GA();
            enviroment_selection();
        }
        //结果的输出
        String dataStr = "";
        for (int i = 0; i < Popsize; i++) {
            for (int j = 0; j < dim_y; j++) {
                System.out.print(Population[i].fitness[j] + "\t");
                dataStr = dataStr + Population[i].fitness[j] + "\t";
            }
            dataStr +="\n";
            System.out.println();
        }

        try {
            File file = new File("E:\\algorithmTest\\res\\nsgaiiPF_"+ testFun +"obj"+dim_y + "var"+dim_x +".txt");
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
	
	//init函数：计算初始种群的适应值，非支配等级，拥挤度距离
    private void init(){
        Population = cal_fitness(Population, Popsize);
        int[] FrontNo = nondoninated_sort(Population, Popsize, dim_y);
        double[] CrowdDis = crowd_distance(Population, Popsize,FrontNo);
        //更新对应位置个体的数值
        for (int i = 0; i < Popsize; i++) {
            FrontRank[i] = FrontNo[i];
            CrowdDistance[i] = CrowdDis[i];
        }
    }
	
	//mating_selection函数：在种群中选择个体来产生下一代种群
    private void mating_selection(){
        //二元竞标赛，挑两个个体，选择更收敛的个体
        for (int i = 0; i < Popsize; i++){
            //随机选择两个个体
            int parent1_index = (int)(Math.random() * Popsize);
            int parent2_index = new Random().nextInt(Popsize);
            //确保个体不同
            while (isequal(Population,parent1_index,parent2_index)){
                parent2_index = new Random().nextInt(Popsize);
            }
            //首先比较选择个体的PF等级
            if (FrontRank[parent1_index] < FrontRank[parent2_index]) {
                for (int j = 0; j < dim_x; j++) {
                    Offspring[i].position[j] = Population[parent1_index].position[j];
                }
            }
            else if(FrontRank[parent1_index] > FrontRank[parent2_index]){
                for (int j = 0; j < dim_x; j++) {
                    Offspring[i].position[j] = Population[parent2_index].position[j];
                }
            }
            //如果等级相同，选择拥挤度距离大的个体
            else {
                if(CrowdDistance[parent1_index] < CrowdDistance[parent2_index]){
                    for (int j = 0; j < dim_x; j++) {
                        Offspring[i].position[j] = Population[parent2_index].position[j];
                    }
                }
                else{
                    for (int j = 0; j < dim_x; j++) {
                        Offspring[i].position[j] = Population[parent1_index].position[j];
                    }
                }
            }
        }
    }

    //isequal函数：判断两个个体是否完全一致
    //输入；种群，个体a的位置，个体b的位置
    //输出：是否相同
    private boolean isequal(Particle[] Pop, int a,int b){
        boolean tag = true;
        for (int i = 0; i < dim_x; i++) {
            if(Population[a].position[i] != Population[b].position[i]){
                tag = false;
                break;
            }
        }
        return tag;
    }

    //operator_GA函数：使用遗传算子产生下一代个体
    private void operator_GA(){
        //模拟二进制交叉
        for (int i = 0; i < Popsize; i+=2){
            double rdc = rand.nextDouble();
            double control_cross;
            if(rdc > 0.5)
                control_cross = Math.pow( 1/( 2 - rdc * 2 ) , 1.0 / ( dis_c + 1 ) );
            else
                control_cross = Math.pow(  rdc * 2 , 1.0 / ( dis_c + 1 ) );
            for(int j = 0; j < dim_x; j++){
                if(rand.nextDouble() < pc){
                    double X1 = Offspring[i].position[j];
                    double X2 = Offspring[i+1].position[j];
                    Offspring[i].position[j] = 0.5*( (1 + control_cross) * X1 + (1 - control_cross) * X2);
                    Offspring[i+1].position[j] = 0.5*( (1 - control_cross) * X1 + (1 + control_cross) * X2);
                }
            }
        }
        //多项式变异
        for (int i = 0; i < Popsize; i++){
            for(int j = 0; j < dim_x; j++){
                if(rand.nextDouble() < pm){
                    double rdm = rand.nextDouble();
                    double control_mutate;
                    if(rdm > 0.5)
                        control_mutate = 1 - Math.pow(  2*( 1 - rdm ) , 1.0 / ( dis_m + 1 ) );
                    else
                        control_mutate = Math.pow(  rdm * 2 , 1.0 / ( dis_m + 1 ) ) - 1;
                    Offspring[i].position[j] = Offspring[i].position[j] + control_mutate * (MAX_X[j] - MIN_X[j]) ;
                }
                //边界控制，防止前面的操作使得X超出范围
                if(Offspring[i].position[j] > MAX_X[j])
                    Offspring[i].position[j] = MAX_X[j];
                else if (Offspring[i].position[j] < MIN_X[j])
                    Offspring[i].position[j] = MIN_X[j];
            }

        }
        //将两个种群合并
        for (int i = 0; i < Popsize; i++) {
            for (int j = 0; j < dim_x; j++) {
                Sumpop[i + Popsize].position[j] = Offspring[i].position[j];
            }
        }
        for (int i = 0; i < Popsize; i++) {
            for (int j = 0; j < dim_x; j++) {
                Sumpop[i ].position[j] = Population[i].position[j];
            }
        }
    }

    //enviroment_selection函数：选择好的个体进入下一代种群
    private void enviroment_selection(){
        //计算种群的等级以及适应度
        int pop_size = Popsize * 2;
        Sumpop = cal_fitness(Sumpop, pop_size);
        int[] FrontNo = nondoninated_sort(Sumpop, pop_size, dim_y);
        double[] CrowdDis = crowd_distance(Sumpop, pop_size,FrontNo);
        //选择个体进入下一代
        int[] UniqLoc = unique(pop_size);
        int[] Next = new int[pop_size];
        int[] Loc = new int[pop_size];
        int max_fNo;
        int p1 = 0;
        int p2 = 0;
        //FrontNo[各个个体的非支配等级;最大非支配等级;需要的最大非支配等级]
        //填入比最大非支配等级小的个体，因为计算过没有超出，所有直接加入下一代种群
        max_fNo = FrontNo[pop_size + 1];
        for(int i = 0; i < pop_size; i++){
            if(FrontNo[i] < max_fNo){
                Next[p1++] = i;
            }
            else if (FrontNo[i] == max_fNo){
                Loc[p2++]  = i;
            }
        }
        //已经找到了p1个下一代个体，还需要 pop_size - P1个
        for (int i = p1; i < Popsize; i++) {
            int P = 0;
            double Min = 0.0;
            for (int j = 0; j < p2; j++) {
                if(Loc[j]== 0) continue;
                if(CrowdDis[Loc[j]] > Min){
                    Min = CrowdDis[Loc[j]];
                    P = j;
                }
            }
            Next[i] = Loc[P];
            Loc[P]=0;
        }
        //前面记录的都是坐标，现在正式将选好的个体赋值给主种群,并记录好选中个体的非支配等级和拥挤度距离
        for (int i = 0; i < Popsize; i++) {
            for (int j = 0; j < dim_x; j++) {
                Population[i].position[j] = Sumpop[Next[i]].position[j];
            }
            for (int j = 0; j < dim_y; j++) {
                Population[i].fitness[j] = Sumpop[Next[i]].fitness[j];
            }
            CrowdDistance[i] = CrowdDis[Next[i]];
            FrontRank[i] = FrontNo[Next[i]];
        }
    }

    //unique函数：确保合并种群的个体都是独一无二的
    //输入：种群大小
    //输出：不重复个体的位置
    private int[] unique(int pop_size){
        int[] Location = new int[Popsize * 2 + 1];
        int p = 1;
        //第一个解必不可能重复
        Location[0] = 0;
        for (int i = 1; i < pop_size; i++) {
            boolean flag = true;
            for (int j = i - 1; j>= 0 ; j--) {
                flag = true;
                for (int k = 0; k < dim_x; k++) {
                    if(Sumpop[i].position[k] != Sumpop[j].position[k]){
                        flag =false;
                        break;
                    }
                }
                //如果在比较了一个解之后发现完全一致，则确定这个解被去掉，可以退出循环
                if(flag) break;
            }
            //如果是独立解，则加入Location里面
            //不是则跳过
            if (!flag){
                Location[p] = i;
                p++;
            }
        }
        Location[Popsize * 2] = p - 1;
        return Location;
    }

    //cal_fitness函数：计算种群个体的适应值
    //输入：种群，个体数量
    //输出：拥有正确适应值的种群
    private Particle[] cal_fitness(Particle[] Pop, int pop_size){
        for (int i = 0; i < pop_size; i ++){
            switch (this.testFun){
                case "ZDT1":
                    double[] fitnessZDT1 = ft.CalFitnessZDT1(Pop[i].position);
                    for (int j = 0; j < dim_y; j++) {
                        Pop[i].fitness[j] = fitnessZDT1[j];
                    }
                    break;
                case "ZDT2":
                    double[] fitnessZDT2 = ft.CalFitnessZDT2(Pop[i].position);
                    for (int j = 0; j < dim_y; j++) {
                        Pop[i].fitness[j] = fitnessZDT2[j];
                    }
                    break;
                case "ZDT3":
                    double[] fitnessZDT3 = ft.CalFitnessZDT3(Pop[i].position);
                    for (int j = 0; j < dim_y; j++) {
                        Pop[i].fitness[j] = fitnessZDT3[j];
                    }
                    break;
                case "DTLZ1":
                    double[] fitnessDTLZ1 = ft.CalFitnessDTLZ1(Pop[i].position, dim_y);
                    for (int j = 0; j < dim_y; j++) {
                        Pop[i].fitness[j] = fitnessDTLZ1[j];
                    }
                    break;
                case "DTLZ3":
                    double[] fitnessDTLZ3 = ft.CalFitnessDTLZ3(Pop[i].position, dim_y);
                    for (int j = 0; j < dim_y; j++) {
                        Pop[i].fitness[j] = fitnessDTLZ3[j];
                    }
                    break;
                case "DTLZ7":
                    double[] fitnessDTLZ7 = ft.CalFitnessDTLZ7(Pop[i].position, dim_y);
                    for (int j = 0; j < dim_y; j++) {
                        Pop[i].fitness[j] = fitnessDTLZ7[j];
                    }
                    break;
            }
        }
        return Pop;
    }

    //nondoninated_sort函数：对种群进行非支配排序，非支配指的是没有个体能够在所有的维度不劣于该个体，且至少有一维优于
    //输入：种群，个体数量，目标的维度
    //输出：[不同的个体的层级，最大的层级]
    private int[] nondoninated_sort(Particle[] Pop,int num,int dim){
        int[] FrontNo = new int[num];
        int[] FrontP  = new int[num + 2];
        double[][] Obj = sort_onedim(Pop,num,dim);
        int Ranked_count = 0;
        int MaxFNo = 0;
        int M = dim;
        while(Ranked_count < num){
            MaxFNo = MaxFNo + 1;
            for (int i = 0; i < num; i ++){
                if(FrontNo[i] == 0){
                    boolean Dominated = false;
                    //因为是从小到大排，所以我们只需要考虑，在他们前面的是不是支配他的，后面的肯定不会支配他们
                    for (int j = i - 1 ; j >= 0; j --){
                        if(FrontNo[j] == MaxFNo ){
                            //从第二维开始比较，计算支配关系
                            int m = 1;
                            while (m <= M-1 && Obj[i][m] >= Obj[j][m])
                                m ++;
                            //比较支配的维度
                            Dominated = m > M - 1;
                            //如果他被支配，则这次标记不用管他
                            if(Dominated)
                                break;
                        }
                    }
                    //如果是非支配的，则分配支配等级
                    if(!Dominated){
                        Ranked_count++;
                        FrontNo[i] = MaxFNo;
                        if(Ranked_count == num/2)
                            FrontP[num + 1] = MaxFNo;
                    }

                }
            }
        }
        //将排序之后的支配等级赋到原来的位置
        for (int i = 0; i < num; i ++){
            FrontP[ Loc[ i ] ] = FrontNo[ i ];
        }
        FrontP[num] = MaxFNo;
        return FrontP;
    }

    //sort_onedim函数：根据第一维的目标值对于种群进行排序，方便后续的非支配排序进行操作
    //输入：种群，个体数目，目标维度
    //输出：按第一维从小到大排的目标值，同时将其原来的位置记录在Loc数组里
    private double[][] sort_onedim(Particle[] Pop, int num, int dim){
        //保存个体的适应值，方便操作
        double[][] Ans = new double[num][dim];
        double[] Obj2 = new double[num];
        double INF = 10000.0;
        for (int i = 0; i < num; i++) {
            Obj2[i] = Pop[i].fitness[0];
        }
        //排序 O（N^2）
        for (int i = 0; i < num; i ++){
            double min = INF;
            int P ;
            P = 0;
            for (int j = 0; j < num; j ++){
                if(Obj2[j] < min){
                    min = Obj2[j];
                    P = j;
                }
            }
            Loc[i] = P;
            Obj2[P] = INF;
        }
        for (int i = 0; i < num; i++) {
            for (int j = 0; j < dim; j++) {
                Ans[i][j] = Pop[ Loc[i] ].fitness[j];
            }
        }
        return Ans;
    }

    //crowd_distance函数：根据非支配等级计算种群个体的拥挤度距离
    //随着代数的增加，非支配级数会慢慢变成一，所有只有前期才会有很多拥挤度距离很大的个体
    //输入：待计算的种群，种群大小，每个个体的非支配等级
    //输出：对应个体的拥挤度距离
    private double[] crowd_distance(Particle[] Pop, int pop_size, int[] FrontNo){
        // Rank =            [1 ~ Front_num]
        // Front = [Front_num;对应的个体在种群的位置]
        double[] CrowdDis = new double[pop_size];
        double INF = 10000.0;
        int max_fNo = FrontNo[pop_size];
        //一层一层非支配等级地计算拥挤度距离
        for(int i = 1; i <= max_fNo; i ++){
            int[] Front = find(FrontNo,i, pop_size);
            int front_num = Front[0];
            if(front_num ==0)
                continue;
            for(int j = 0; j < dim_y; j ++){
                int[] Rank = sort_front(Pop,Front,j);
                //每一层最边上的个体我们需要保留，这是有一个目标优化的很好的个体
                CrowdDis[ Front[ Rank[0] ] ] = INF;
                CrowdDis[ Front[ Rank[front_num - 1] ] ] = INF;
                //如果这一层个体不超过两个，那么他们都是最边上的个体，直接进入下一代
                if(front_num == 1 || front_num ==2 ) continue;
                //计算中间个体的拥挤度距离
                for(int k = 1; k < front_num - 1; k ++){
                    CrowdDis[ Front[ Rank[k] ] ] = CrowdDis[ Front[ Rank[k] ] ]
                            + (Pop[ Front[ Rank[k + 1] ] ].fitness[j] - Pop[ Front[ Rank[k - 1] ] ].fitness[j] )
                            /(Zmax[j] - Zmin[j]);
                }

            }
        }
        return CrowdDis;
    }

    //find函数：找到输入的非支配排序等级的解的位置
    //输入：所有的非支配的数组，需要的非支配等级，种群大小
    //输出：[解的数量，解的位置]
    private int[] find(int[] FrontNo,int point,int pop_size){
        //数组的第一位记录个数
        int num = 0;
        int p = 1;
        int[] Front = new int[pop_size];
        for(int i = 0; i < pop_size; i++){
            if(FrontNo[i]==point){
                Front[p] = i;
                num ++;
                p++;
            }
        }
        Front[0] = num;
        return Front;
    }

    //sort_front函数：根据输入非支配等级在Dim这一维进行排序
    //输入：种群，待排序个体的位置，待排序的维度
    //输出：排序后的坐标
    private int[] sort_front(Particle[] Pop, int[] Front, int dim){
        //首先取出该层的解的这一维的坐标
        int front_num = Front[0];
        int[] Rank = new int[front_num];
        double[] Obj = new double[front_num];
        double INF = 100000.0;
        for (int i = 0; i < front_num; i ++){
            int P = Front[i + 1];
            Obj[i] = Pop[P].fitness[dim];
        }
        //排序得到从小到大的点的坐标对于的解位置
        //在每次该算法过程中找到该非支层在这一维的最大值最小值，进行之后的归一化操作
        for (int i = 0; i < front_num; i ++){
            double min = INF;
            int j;
            for (j = 0; j < front_num; j ++){
                if(Obj[j] < min){
                    Rank[i] = j;
                    min = Obj[j];
                }
            }
            Obj[ Rank[ i ] ] = INF;
        }
        for (int i = 0; i < front_num; i++) {
            Rank[i] += 1;
        }
        //记录这一层的最大值和最小值，用来进行归一化
        //每次调用，记录的都是这一层的最大最小值
        Zmin[dim] = Pop[ Front[Rank[0]] ].fitness[dim];
        Zmax[dim] = Pop[ Front[Rank[front_num -1]] ].fitness[dim];
        return Rank;
    }

    public double calIGD(double[][] pf){

        double sumD = 0;
        int samplePointNum = pf.length;
        for (int i = 0 ; i < samplePointNum; i ++){
            double minDis = Double.MAX_VALUE;
            for(int j = 0 ; j < Popsize; j ++){
                double dis = 0;
                for (int k = 0; k < dim_y; k++){
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
}
