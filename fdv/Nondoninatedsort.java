package fdv;

public class Nondoninatedsort {
    int[] Loc;
    public Nondoninatedsort(int Num){
        this.Loc = new int[Num];
    }
    public int[] NDsort(FUZZY.Particle[] Pop, int Num, int Dim){
        int[] FrontNo ;
        if(Dim < 3 || Num < 500){
            FrontNo = ENS_SS(Pop,Num,Dim);
        }
        else{
            FrontNo = T_ENS(Pop,Num);
        }
        return FrontNo;
    }

    // ENS_SS函数：对种群进行非支配排序，非支配指的是没有个体能够在所有的维度不劣于该个体，且至少有一维优于
    // 输入：种群，个体数量，目标的维度
    // 输出：[不同的个体的层级，最大的层级]
    private int[] ENS_SS(FUZZY.Particle[] Pop, int Num, int Dim){
        int[] FrontNo = new int[Num];
        int[] FrontP  = new int[Num + 2];
        double[][] Obj = SortOneDim(Pop,Num,Dim);

        int Ranked_count = 0;
        int MaxFNo = 0;
        int M = Dim;
        while(Ranked_count < Num){
            MaxFNo = MaxFNo + 1;
            for (int i = 0; i < Num; i ++){
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
                        if(Ranked_count == Num/2)
                            FrontP[Num + 1] = MaxFNo;
                    }

                }
            }
        }
        //将排序之后的支配等级赋到原来的位置
        for (int i = 0; i < Num; i ++){
            FrontP[ Loc[ i ] ] = FrontNo[ i ];
        }

        FrontP[Num] = MaxFNo;

        return FrontP;
    }
    // SortOneDim函数：根据第一维的目标值对于种群进行排序，方便后续的非支配排序进行操作
    // 输入：种群，个体数目，目标维度
    // 输出：按第一维从小到大排的目标值，同时将其原来的位置记录在Loc数组里
    private double[][] SortOneDim(FUZZY.Particle[] Pop, int Num, int Dim){
        // 保存个体的适应值，方便操作
        double[][] ans = new double[Num][Dim];
        double[] Obj2 = new double[ Num ];
        double INF = 10000.0;
        for (int i = 0; i < Num; i++) {
            Obj2[i] = Pop[i].fitness[0];
        }
        // 对每个个体进行非支配排序
        for (int i = 0; i < Num; i ++){
            double min = INF;
            int P ;
            P = 0;
            for (int j = 0; j < Num; j ++){
                if(Obj2[j] < min){
                    min = Obj2[j];
                    P = j;
                }
            }
            Loc[i] = P;
            Obj2[P] = INF;
        }
        for (int i = 0; i < Num; i++) {
            for (int j = 0; j < Dim; j++) {
                ans[i][j] = Pop[ Loc[i] ].fitness[j];
            }
        }
        return ans;
    }

    private int[] T_ENS(FUZZY.Particle[] Pop, int Num){
        int[] FrontNo = new int[Num];
        return FrontNo;
    }
}

