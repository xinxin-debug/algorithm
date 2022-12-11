package testTADF;

import java.util.*;


public class populations {
    private int N = 100;//种群的大小
    private int maxFE = 600;//迭代次数
    private int dimentionx = 2;//个体决策变量数
    private double pc = 1;//发生交叉行为的概率
    private double pm = 1;//发生变异的概率
    private static individual[] pop;//种群的个体
    private static UniformPoint V;//权重向量对象
    private static double[][] R;//权重向量
    private double Z[];//参考点
    private static individual[] CA;//收敛性归档集
    private int[] c_count=new int[N];//记录每个权重向量上已经加入CA归档集的个体数
    private static individual[] DA;//多样性归档集
    private int NCA = N;//收敛性归档集的大小
    private int NDA = N;//多样性归档集的大小
    private static individual[] Cpop;//子种群
    private static individual[] Ppop;//交配种群

    public populations() { //构造函数
        V = new UniformPoint();
        R = V.UniformPoint(N, 2);//初始化生成N个2维的均匀权重向量 R[N][2]
    }

    //**********运行*******
    //**迭代进行优化，首先初始化种群，之后迭代进行CA与DA操作（收敛性、分布性归档集）
    //测试问题：MMF1，位于individual
    //测试结果：输出dec和obj空间的数据
    public void run() {
        init();

        CA=UpdateCA(CA, pop, NCA, R, Z);
        DA=UpdateDA(DA,pop,CA,NDA,R,Z);

        int itr=0;
        double tim;
        individual[] Parents;
        while (itr<maxFE){
            System.out.print("\r");
            tim=itr/(double)maxFE*100;
            System.out.print((int)tim + "% --------");

            Parents=MatingSelection(CA,DA);
            individual[] popi;
            popi=OperatorGA(Parents);
            relese(popi);

            CA=UpdateCA(CA, popi, NCA, R, Z);
            DA=UpdateDA(DA,popi,CA,NDA,R,Z);

            pop=DA;
            itr++;
           // data_out.dataout(pop,N,itr);//导出数据
        }
        System.out.print("\r");
        System.out.println((100 + "% --------"));
        System.out.println("dec:");
        for (int i=0;i<N;i++){      //输出最终每个个体的目标空间结果
            System.out.println(pop[i].obj[0]+"\t"+pop[i].obj[1]);
        }
    }



    //*************************************
    //********* 更新pop ********************
    //对新产生的种群进行属性更新(包括目标值）
    //输入：种群
    //输出：更新后的种群
    public void relese(individual[] pop){
        for (int i=0;i<pop.length;i++){
                pop[i].fitness = 0;
                pop[i].distance = 0;
                pop[i].ri = -1;
                pop[i].rj = -1;
                pop[i].S = 0;
                pop[i].n_var = new double[100];//暂定为有100个决策变量
                pop[i].number = 0;
            pop[i].calculate_MMF1();//计算该个体的目标空间值
        }
    }

    //*********************************************************
    //********* 初始化每个个体,及他们的决策空间、目标空间，和参考点Z值**
    //输出：初始化的新种群
    private void init() {
        pop = new individual[N];//初始化N个个体individual
        pop[0] = new individual(0);//记录下第一个个体目标空间值，方便后面比较出最小的
        double a = pop[0].obj[0];
        double b = pop[0].obj[1];

        for (int i = 1; i < N; i++) {
            pop[i] = new individual(i);//调用构造函数，对种群中的每一个个体进行初始化操作 产生一个随机点,并计算该点对应目标空间的值
            //找到目标空间的一个最小向量作为参考点Z
            if (pop[i].obj[0] < a) {
                a = pop[i].obj[0];
            }
            if (pop[i].obj[1] < b) {
                b = pop[i].obj[1];
            }
        }
        this.Z = new double[2];
        Z[0] = a;
        Z[1] = b;
    }

    //********************************************************************
    //********************* 选择 ******************************************
    //从CA、DA归档集挑选3/4、1/4的个体到下一代，二元锦标赛选择
    //输入：种群CA、DA
    //输出：选择好的种群
    public individual[] MatingSelection(individual[] CAa,individual[] DAa){
        individual[] CAx=new individual[CAa.length];
        individual[] DAx=new individual[DAa.length];

        //////
        int []a=new int[CAa.length];
        int []b=new int[DAa.length];
        for(int i=0;i< CAa.length;i++){
            a[i]=i;
            b[i]=i;
        }
        Random r=new Random();
        int c1,c2,ctem;
        for(int i=0;i< CAa.length;i++){
            c1= r.nextInt(NCA);
            c2= r.nextInt(NCA);
            while(c1==i)
                c1= r.nextInt(NCA);
            ctem=a[c1];
            a[c1]=a[i];
            a[i]=ctem;

            while(c2==i)
                c2= r.nextInt(NCA);
            ctem=b[c2];
            b[c2]=b[i];
            b[i]=ctem;
        }

        for (int i=0;i<CAa.length;i++){//拷贝
            CAx[i]=new individual(CAa[i]);
        }
        for (int i=0;i<DAa.length;i++){//拷贝
            DAx[i]=new individual(DAa[i]);
        }

        individual[] Parent_M=new individual[3*NCA/4+NDA/4];
        for (int i=0;i<3*NCA/4;i++){
            if(CAx[a[i]].fitness<CAx[b[i]].fitness){
                Parent_M[i]=new individual(CAx[a[i]]);
            }else {
                Parent_M[i]=new individual(CAx[b[i]]);
            }
        }
        Random r3=new Random();
        int p3;
        int []sav = new int[NCA/4]; ///
        for (int i=0;i<NCA/4;i++){
            p3=r3.nextInt(NDA);
            int sam=1;
            while(sam==1) {
                int j=0;
                for (; j < i; j++) {
                    if (p3 == sav[j]) {
                        p3=r3.nextInt(NDA);break;
                    }
                }
                if(j==i||j==0)
                    sam=0;
            }
            sav[i]=p3;              ///
            Parent_M[3*NCA/4+i]=new individual(DAx[p3]);
        }
        return Parent_M;
    }

    //**************************************************
    //*************** 遗传操作 **************************
    //输入：父代种群
    //输出：子代种群
    public individual[] OperatorGA(individual[] Parent){
        individual[] Parent1=new individual[Parent.length/2];
        individual[] Parent2=new individual[Parent.length/2];
        double[][] offspring_dec1=new double[Parent.length][2];
        double[][] offspring_dec2=new double[Parent.length/2][2];
        //将Parent分开
        for (int i=0;i<Parent.length/2;i++){
            Parent1[i]=Parent[i];
        }

        for (int i=Parent.length/2;i<Parent.length;i++){
            Parent2[i-Parent.length/2]=Parent[i];
        }

        for (int i=0;i<Parent.length/2;i++){
            double[] beta=new double[2];//beta[0] beta[1]

            for (int j=0;j<2;j++){
                double mu=Math.random();
                if(mu<=0.5){
                    beta[j]=Math.pow(2*mu,0.0476);
                }else {
                    beta[j]=Math.pow((2-2*mu),-0.0476);
                }
                beta[j]=beta[j]*Math.pow(-1,((int)(10*Math.random()))%2);
                if(Math.random()<=0.5){//0.5
                    beta[j]=1;
                }
                offspring_dec1[i][j]=(Parent1[i].dec[j]+Parent2[i].dec[j])/2 + beta[j]*(Parent1[i].dec[j]-Parent2[i].dec[j])/2;
                offspring_dec2[i][j]=(Parent1[i].dec[j]+Parent2[i].dec[j])/2 - (1-beta[j])*(Parent1[i].dec[j]-Parent2[i].dec[j])/2;
            }
        }
        //合并子代
        for (int i=Parent.length/2;i<Parent.length;i++){
            for (int j=0;j<2;j++){
                offspring_dec1[i][j]=offspring_dec2[i-Parent.length/2][j];
            }
        }

        double[] Lower=new double[2];
        double[] Upper=new double[2];
        Lower[0]=1; Lower[1]=-1;
        Upper[0]=3; Upper[1]=1;
        for (int i=0;i<Parent.length;i++){
            for (int j=0;j<2;j++){
                if(offspring_dec1[i][j]<Lower[j])//保证范围在Lower~Upper之间
                    offspring_dec1[i][j]=Lower[j];
                if(offspring_dec1[i][j]>Upper[j])
                    offspring_dec1[i][j]=Upper[j];

                boolean[] Site=new boolean[2];
                double[] mu=new double[2];
                boolean[] temp=new boolean[2];

                if(Math.random()<=0.5)
                    Site[j]=true;
                else
                    Site[j]=false;
                mu[j]=Math.random();
                temp[j]=Site[j] && mu[j]<=0.5;
                if(temp[j]){
                    offspring_dec1[i][j]=offspring_dec1[i][j]+(Upper[j]-Lower[j]) * (Math.pow(2*mu[j]+(1-2*mu[j])*Math.pow(1-(offspring_dec1[i][j]-Lower[j])/(Upper[j]-Lower[j]),20+1),0.0476)-1);
                }
                temp[j]=Site[j] && mu[j]>0.5;
                if(temp[j]){
                    offspring_dec1[i][j]=offspring_dec1[i][j]+(Upper[j]-Lower[j]) * (1-Math.pow(2*(1-mu[j])+2*(mu[j]-0.5)*Math.pow(1-(Upper[j]-offspring_dec1[i][j])/(Upper[j]-Lower[j]),20+1),0.0476));
                }
                //复制回去
                Parent[i].dec[j]=offspring_dec1[i][j];

            }
        }
        return Parent;
    }

    //**********************************************************
    //***************** 更新CA归档集 *****************************
    //对权重向量上绑定的个体根据决策空间距离和适应度进行选择
    //输入：CA归档集、子代种群pop，群体大小、权重向量
    //输出：新的CA
    private individual[] UpdateCA(individual[] CAx, individual[] popx, int NCAx, double[][] Rx, double[] Zx) {
        int Nx = N;//记录当前Pop里的个体数
        int L = 3;//测试出来是3
        double[][] distancR;//个体到权重向量的距离
        double[][] d;//记录个体决策空间之间的欧式距离

        if(CAx!=null){
            Nx=2*N;
        }
        individual[] deepcopy = new individual[Nx]; //种群和CA合并为Pop
        for (int i = 0; i < N; i++) {
            deepcopy[i] =new individual(popx[i]);
        }
        if (CAx != null) { //CA归档集有个体才合并到Pop里
            for (int i = N; i < 2 * N; i++) {
                deepcopy[i] = new individual(CAx[i-N]);
            }
        }
        individual[] Pop=new individual[Nx];

        for (int i=0;i<Nx;i++){//拷贝
            Pop[i]=new individual(deepcopy[i]);
        }

        //绑定向量
        distancR = new double[Nx][Rx.length];//每个个体分别到每个权重向量的距离  到这里才放new是因为在前面放的话Nx值还没算出来
        double[] sum=new double[Nx];
        for (int i = 0; i < Nx; i++) {
            sum[i] = Math.sqrt(Pop[i].obj[0]*Pop[i].obj[0]+Pop[i].obj[1]*Pop[i].obj[1]);
        }
        for (int i = 0; i < Nx; i++) {
            double mindis = 100000;//初始化为最大
            for (int j = 0; j < Rx.length; j++) {
                double x1 = Pop[i].obj[0], x2 = R[j][0];
                double y1 = Pop[i].obj[1], y2 = R[j][1];
                double cos = 0;
                cos = (x1 * x2 + y1 * y2) / (Math.sqrt(x1 * x1 + y1 * y1) * Math.sqrt(x2 * x2 + y2 * y2));
                distancR[i][j] = sum[i]* Math.sqrt(1 - cos * cos);
                if (distancR[i][j] < mindis) {
                    mindis = distancR[i][j];
                    Pop[i].ri = j;//记录每个个体绑定的最近权重向量序号
                }
            }
        }
        //计算每个个体决策向量之间的 欧式距离
        d = new double[Nx][Nx];
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Nx; j++) {
                if(i!=j){
                    d[i][j] = pdist2(Pop[i].dec, Pop[j].dec);
                }
            }
        }
        //计算个体适应度
        indator_C(Pop);//更新到了Pop[].fitness里面
        //更新
        int[][] W_I=new int[Rx.length][Nx];//记录每个权重向量上所有按适应度从小到大排序好的个体序号
        int[][] C=new int[Rx.length][Nx];//记录在每个权重向量上加入收敛归档集CA的个体序号
        for (int i=0;i<Rx.length;i++) {
            for (int j = 0; j < Nx; j++) {
                C[i][j]=-1;
            }
        }
        int[] w_count=new int[Rx.length];//记录每个权重向量上目前的个体数
        c_count=new int[Rx.length];
        int[] Index = new int[Pop.length];//记录按个体决策空间欧式距离排序好后对应的索引
        int sum_count=0;//统计所有收敛归档集个体的总个数

        MyComparator myComparator1 = new MyComparator();//new一个自定义比较的对象
        myComparator1.setIndex(2);//设置按适应度排序个体
        Arrays.sort(Pop, myComparator1);
        for (int i=0;i<Rx.length;i++){
            int s=0;
            for (int j=0;j<Nx;j++){//依次遍历排序好的个体
                if(Pop[j].ri==i){
                    W_I[i][s++]=j;//第i个权重向量上依次绑定了满足条件的个体序号
                    w_count[i]++;
                }
            }
        }
        for (int i=0;i<Rx.length;i++){
            if(w_count[i]!=0){
                C[i][0]=W_I[i][0];//依次从有个体的权重向量上选择一个个体加入收敛归档集，用C来记录，相应也要从W_I里删去这个个体
                for (int j=0;j<w_count[i];j++){//更新W_I
                    W_I[i][j]=W_I[i][j+1];
                }
                w_count[i]--;
                c_count[i]++;
                sum_count++;
            }
        }

        int jj=0;//记录权重向量序号
        while (sum_count<NCA){
            while (w_count[jj]==0){//当第jj个向量上的个体都加入归档集后，就换到下一条向量    确保当前向量上是有个体的
                jj=(jj+1)%Rx.length;
            }

            int kk=0;
            int Temp=0;
            while (kk<w_count[jj]&&Temp==0){
                int tempx=0;
                Index=ArrayHelper.Arraysort(d[W_I[jj][kk]]);//W_I[jj][kk]:第jj个向量上的第kk个个体序号
                for (int k=0;k<L;k++){
                    for (int i=0;i<Rx.length;i++){
                        for (int j=0;j<c_count[i];j++){
                            if(Index[k]==C[i][j]){//如果靠近W_I[jj][kk]个体的前L个个体有在收敛归档集中
                                tempx++;
                            }
                        }
                    }
                }
                if(tempx==0){//如果靠近W_I[jj][kk]个体的前L个个体没有一个在收敛归档集中
                    C[jj][c_count[jj]]=W_I[jj][kk];//将当前W_I[jj][kk]并入归档集中
                    for (int j=kk;j<w_count[jj];j++){//更新W_I
                        W_I[jj][j]=W_I[jj][j+1];
                    }
                    w_count[jj]--;
                    c_count[jj]++;
                    sum_count++;
                    Temp=1;
                }
                else {
                    kk++;
                }
            }
            if(Temp==0){//如果在当前第jj个向量上 找的所有kk个体的前L个个体都有在收敛归档集中的
                C[jj][c_count[jj]]=W_I[jj][0];//就将排序好的第一个并入归档集
                for (int j=0;j<w_count[jj];j++){//更新W_I
                    W_I[jj][j]=W_I[jj][j+1];
                }
                w_count[jj]--;
                c_count[jj]++;
                Temp=1;
                sum_count++;
            }
            jj=(jj+1)%Rx.length;
        }

        int p=0;
        CA=new individual[NCA];
        for (int i=0;i<Rx.length;i++){
            for (int j=0;j<c_count[i];j++){
                CA[p++]=new individual(Pop[C[i][j]]);
            }
        }
        return CA;
    }

    //**********************************************************
    //********** 计算欧式距离 ************************************
    //输入：两个向量
    //输出：两个向量距离
    public double pdist2(double[] vector1, double[] vector2) {
        double distance = 0;
        if (vector1.length == vector2.length) {
            for (int i = 0; i < vector1.length; i++) {
                double temp = Math.pow((vector1[i] - vector2[i]), 2);
                distance += temp;
            }
            distance = Math.sqrt(distance);
        }
        return distance;
    }

    //****************************************
    //******** 计算个体适应度 ******************
    //输入：种群
    //输出：计算种群的个体属性：适应度
    public void indator_C(individual[] Popx) {
        int[][] Dominate = new int[Popx.length][Popx.length];//记录个体间的支配关系
        for (int i = 0; i < Popx.length; i++) {
            for (int j = 0; j < Popx.length; j++) {
                if (j != i) {
                    if ((Popx[i].obj[0] <= Popx[j].obj[0]) && (Popx[i].obj[1] <= Popx[j].obj[1])) {
                        Dominate[i][j] = 1;//i支配j  这里只记录了纵轴支配横轴
                        Popx[i].S += 1;//计算SPEA2算法中的S值
                    }
                }
            }
        }
        //计算SPEA2算法中的R值
        int[] Rx = new int[Popx.length];
        Arrays.fill(Rx, 0);
        for (int j = 0; j < Popx.length; j++) {
            for (int i = 0; i < Popx.length; i++) {
                if (j != i) {
                    if (Dominate[i][j] == 1) {//如果j被i支配
                        Rx[j] += Popx[i].S;//将支配j的个体的S值都加起来
                    }
                }
            }
        }
        //计算SPEA2算法中的拥挤度D值
        double[] Dx=new double[Popx.length];
        crowding a=new crowding();
        Dx=a.crowding(Popx);//计算D值

        for (int i=0;i<Popx.length;i++){//计算适应度
            Popx[i].fitness=Rx[i]+Dx[i];
            Popx[i].S=0;//归零
        }
    }


    //**********************************************************
    //************** 更新多样性归档集DA ***************************
    //首先进行全部个体的非支配排序，在最大需要的支配层内进行选择：
    //目标空间距离大的全部选择，对于不足DA大小：对每个权重考虑：根据CA中权重向量绑定的个体的情况，个体越少该权重向量优先选择；
    // 选择的权重向量上的多个个体选择：找到支配层最少的一个个体，在目标空间中距离该个体最近的个体中选择拥挤度最大的
    //输入：DA、子代种群、CA的绑定权重信息、种群大小、权重向量
    //输出：更新的DA
    public individual[] UpdateDA(individual[] DAx, individual[] popx,individual[] CAx,int NDAx, double[][] Rx, double[] Zx){
        int Nx = N;//记录当前Pop里的个体数
        double[][] distancR;//个体到权重向量的距离
        double[][] d;//记录个体目标空间之间的欧式距离
        if(DAx!=null){
            Nx=2*N;
        }
        individual[] deepcopy = new individual[Nx]; //种群和DA合并为Pop
        for (int i = 0; i < N; i++) {
            deepcopy[i] = new individual(popx[i]);
        }
        if (DAx != null) { //DA归档集有个体才合并到Pop里
            for (int i = N; i < 2 * N; i++) {
                deepcopy[i] = new individual(DAx[i-N]);
            }
        }
        individual[] Pop=new individual[Nx];
        for (int i=0;i<Nx;i++){
            Pop[i]=new individual(deepcopy[i]);
        }

        //计算每个个体目标空间的欧式距离Obj
        d = new double[Nx][Nx];
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Nx; j++) {
                if(i!=j){
                    d[i][j] = pdist2(Pop[i].obj, Pop[j].obj);
                }
            }
        }
        //计算Dec拥挤度
        crowding b=new crowding();
        b.crowding(Pop);
        double avg_crowd_dist=b.avg_crowd_dist;//平均拥挤距离，个体的拥挤距离在个体对象的属性里有
        fast_Nondominated_sort(Pop);//快速非支配排序
        // ==============================================================================
        int k=0;//记录目前能完整按层加入归档集的个体总数
        int MaxFNo=0;//记录前几层能完全加入归档集的总层数
        while (k<NDAx){
            MaxFNo++;
            int sumi=0;
            for (int i=0;i<Nx;i++){
                if(Pop[i].rank==MaxFNo){
                    sumi++;//记录当前MaxFNo层共有多少个体
                }
            }
            k=k+sumi;
        }

        int countx=1;//记录当前层数
        while (countx<MaxFNo){//第MaxFNo层暂时先不管
            for (int i=0;i<Nx;i++){
                if(Pop[i].rank==countx){//找到第countx层的所有个体
                    if(Pop[i].distance>0.5*avg_crowd_dist){// 满足条件则并入多样性归档集，不满足条件的进入LastSelection
                        Pop[i].rank=0;//把该个体处于的层数记为0层，即将该个体从第countx层删除  最后若有rank不为0的则表示剩下没被选择的个体
                    }
                }
            }
            countx++;
        }

        //绑定向量
        distancR = new double[Nx][Rx.length];//每个个体分别到每个权重向量的距离  到这里才放new是因为在前面放的话Nx值还没算出来
        int[] F_count=new int[Rx.length];//计算目前每个向量上有多少个体
        double[] sum=new double[Nx];
        for (int i = 0; i < Nx; i++) {
            sum[i] = Math.sqrt(Pop[i].obj[0]*Pop[i].obj[0]+Pop[i].obj[1]*Pop[i].obj[1]);
        }
        for (int i = 0; i < Nx; i++) {
            double mindis = 100000;//初始化为最大
            for (int j = 0; j < Rx.length; j++) {
                double x1 = Pop[i].obj[0], x2 = R[j][0];
                double y1 = Pop[i].obj[1], y2 = R[j][1];
                double cos = 0;
                cos = (x1 * x2 + y1 * y2) / (Math.sqrt(x1 * x1 + y1 * y1) * Math.sqrt(x2 * x2 + y2 * y2));
                distancR[i][j] = sum[i]* Math.sqrt(1 - cos * cos);
                if (distancR[i][j] < mindis) {
                    mindis = distancR[i][j];
                    Pop[i].rj = j;//记录每个个体绑定的最近权重向量序号
                }
            }
        }
        for (int i=0;i<Nx;i++){//计算目前每个向量上有多少个体
            for (int j=0;j<Rx.length;j++){
                if(Pop[i].rj==j&&Pop[i].rank<=MaxFNo){
                    F_count[j]++;
                }
            }
        }
        int sumx=0;//记算已经归入多样性归档集DA的个体总个数
        for (int i=0;i<Nx;i++){
            if(Pop[i].rank==0)
                sumx++;
        }
        int NQ=NDAx-sumx;//记录剩余需要加入DA的个体数
        //LastSelection    对DA中临界层和前面层的未选择个体进行挑选，先看CA中哪些向量绑定个体少，逐次增加；对这些向量看
                            //DA中选定的个体是不是也少，再在这个向量绑定的个体中挑选；先选出支配层最小的，找它附近拥挤度大的
        //这步绑定权重向量可以忽略，前面每个个体绑定的信息都算出来了
        int[] count=new int[Rx.length];//统计每个权重向量上加入了收敛档案CA的个体数
        count=c_count;//复制一下CAsel
        int[] D=new int[Nx];//存储加入多样性归档集DA里个体的序号   D[i]若为1就表示个体i加入了归档集
        int itr=1;
        int[] itr_temp=new int[Rx.length];//记录count(在CA里绑定的)里面<itr的向量序号
        int[] itr_count=new int[Rx.length];//记录满足最终条件的向量序号
        int sumnq=0;//记录从剩下的里面共选了多少个体

        while (sumnq<NQ){
            int s=0,p=0;
            for (int i=0;i<Rx.length;i++){
                if(count[i]<itr){
                    itr_temp[s++]=i;
                }
            }
            for (int i=0;i<s;i++){
                if(F_count[itr_temp[i]]<itr){
                    itr_count[p++]=itr_temp[i];//记录满足最终条件的向量序号
                }
            }

            for (int i=0;i<p;i++){//依次对找到的向量进行操作
                int[] temp_s=new int[Nx];
                int q=0;
                int empty=1;
                for (int j=0;j<Nx;j++){
                    if(Pop[j].rank!=0&&Pop[j].rank<=MaxFNo){//找到剩下个体
                        if(itr_count[i]==Pop[j].rj){//若itr_count(i)中有上面非支配排序剩余个体所绑定的向量序号
                            temp_s[q++]=j;//记录该个体的序号
                            empty=0;
                        }
                    }
                }
                if(empty==1)//如果temp_s个体集为空则继续for循环
                    continue;

                if(q==1){//如果只有一个个体满足条件，直接加入归档集
                    D[temp_s[0]]=1;
                    Pop[temp_s[0]].rank=0;
                }else {
                    int min=10000;
                    int temp_A=0;
                    for (int j=0;j<q;j++){
                        if(Pop[temp_s[j]].rank<min){
                            min=Pop[temp_s[j]].rank;
                            temp_A=temp_s[j];//找到temp_s个体集里面层数最小的一个个体序号
                        }
                    }
                    double max=-1;
                    int maxi=0;
                    for (int j=0;j<q;j++){//找到temp_s个体集里面离temp_A个体拥挤度最大的个体序号
                        if(d[temp_A][temp_s[j]]>max){
                            max=d[temp_A][temp_s[j]];
                            maxi=temp_s[j];
                        }
                    }
                    D[maxi]=1;//加入归档集
                    Pop[maxi].rank=0;
                }
                sumnq=0;
                for (int j=0;j<Nx;j++){
                    if(D[j]==1)
                        sumnq++;
                }
                if(sumnq==NQ)
                    break;
            }
            itr++;
        }
        int n=0;
        DA=new individual[NDAx];
        for (int i=0;i<Nx;i++){
            if(Pop[i].rank==0){
                DA[n++]=new individual(Pop[i]);
            }
        }
        return DA;
    }
    //*****************************************************
    //************** 快速非支配排序 *************************
    //输入：种群
    //输出：种群的pareto支配等级
    private void fast_Nondominated_sort(individual[] popx){

        ArrayList<ArrayList<individual>> F=new ArrayList<>();//帕累托等级列表
        //快速非支配排序，
        ArrayList<individual> fx=new ArrayList<>();//x等级的列表
        for(int i=0;i<popx.length;i++){
            popx[i].Sp=new ArrayList<>();
            popx[i].nq=0;
            for (int j=0;j<popx.length;j++){
                if (j!=i){
                    if (popx[i].obj[0]<=popx[j].obj[0] && popx[i].obj[1]<=popx[j].obj[1]){//两个函数均小于
                        //支配关系 i支配j
                        if (popx[i].obj[0]==popx[j].obj[0] && popx[i].obj[1]==popx[j].obj[1]){
                            //对于相等的来说，没有任何操作才是对的
                        }else{
                            popx[i].Sp.add(popx[j]);//就把j个体加入i个体的个体支配个体的集合
                        }
                    }
                    else if(popx[i].obj[0]>=popx[j].obj[0] && popx[i].obj[1]>=popx[j].obj[1]){//两个函数均大于
                        //被支配关系 i被j支配
                        if (popx[i].obj[0]==popx[j].obj[0] && popx[i].obj[1]==popx[j].obj[1]){
                            //对于相等的来说，没有任何操作才是对的
                        }else{
                            popx[i].nq++;//就把i个体被支配的次数+1
                        }
                    }
                }
            }
            if (popx[i].nq==0){//如果没有人可支配i
                popx[i].rank=1;
                fx.add(popx[i]);//非支配个体
            }
        }
        F.add(fx);//添加第一层支配层的个体列表
        int i=0;
        while(F.get(i).size() !=0){ //第i层有个体，找下一层
            //开始一层一层的去掉帕累托层进行分级
            fx=new ArrayList<>();   //存放下一层个体
            for (int j=0;j<F.get(i).size();j++){    //第i层的每一个个体
                //实际在整个种群中的位置参数
                for (int m=0;m<F.get(i).get(j).Sp.size();m++){      //get列表里面第i行j列
                    individual h=F.get(i).get(j).Sp.get(m);
                    h.nq=h.nq-1;//相当于第一个被删去了。
                    if (h.nq==0){
                        h.rank=i+2;
                        fx.add(h);//新的一轮开始了
                    }
                }
            }
            F.add(fx);
            i++;
        }
        for (int j=0;j<popx.length;j++){
            popx[j].Sp=new ArrayList<>();//归零
        }
    }



    public static void main(String[] args) {
        populations p = new populations();
        p.run();
    }
}
