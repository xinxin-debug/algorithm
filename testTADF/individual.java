package testTADF;

import java.util.*;
import java.util.ArrayList;
import java.util.Random;


public class individual{
    static int dimension=2;//维度
    static int M=2;//目标数
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

    int number;
    double high=1.0;
    double low=0.0;

    public individual(int i){
        this.init(i);
    }


    //*******个体初始化*******
    public void init(int i){
        this.fitness=0;
        this.distance=0;
        this.ri=-1;
        this.rj=-1;
        this.S=0;
        this.n_var=new double[100];//暂定为有100个决策变量
        this.dec=new double[dimension];
        this.dec[0]=2*Math.random()+1;//x1属于1~3
        this.dec[1]=2*Math.random()-1;//x2属于-1~1
        this.number=0;
        this.calculate_MMF1();//计算该个体的目标空间值
    }


    //深拷贝构造函数
    public individual(individual that){
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
        this.dec=new double[2];
        this.obj=new double[2];
        for(int i=0;i<2;i++) {
            this.dec[i] = that.dec[i];
            this.obj[i] = that.obj[i];
        }
    }


    //********** MMF1 *******
    //计算函数目标值
    public void calculate_MMF1(){ //根据决策空间计算MMF1问题的目标空间值,这里还不是适应度
        //修复Dec里无效的解(即保证决策空间的值在指定范围内)
        if(this.dec[0]<1){
            this.dec[0]=1;
        }
        if(this.dec[0]>3){
            this.dec[0]=3;
        }
        if(this.dec[1]<-1){
            this.dec[1]=-1;
        }
        if(this.dec[1]>1){
            this.dec[1]=1;
        }
        //计算目标空间值
        this.obj=new double[M];
        this.obj[0]=Math.abs(this.dec[0]-2);
        this.obj[1]=1 - Math.sqrt(this.obj[0]) + 2*(this.dec[1]-Math.sin(6*Math.PI*this.obj[0]+Math.PI))*(this.dec[1]-Math.sin(6*Math.PI*this.obj[0]+Math.PI));
    }
}
