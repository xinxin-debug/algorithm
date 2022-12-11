package testTADF;

import java.util.Arrays;

//计算决策空间拥挤度和SPEA2算法里的D值
public class crowding {
    double avg_crowd_dist;//平均拥挤度
    private int dimentionx = 2;//个体决策变量数



    public double[] crowding(individual[] Pop) {
        individual[] Popx=new individual[Pop.length];
        for (int i=0;i<Pop.length;i++){
            Popx[i]=new individual(Pop[i]);//拷贝
            Popx[i].number=i;//后面会进行排序，所以先记录一下原来个体的序号
        }
        MyComparator myComparator = new MyComparator();//new一个自定义比较的对象
        for (int j = 0; j < dimentionx; j++) {
            myComparator.setIndex(j);//按决策变量第j维从小到大排序
            Arrays.sort(Popx, myComparator);
            double f_max, f_min;
            f_min = Popx[0].dec[j];
            f_max = Popx[Popx.length - 1].dec[j];
            if(Popx.length==1){//如果个体只有一个，就把排序信息置为1
                Popx[0].n_var[j]=1;
            }else {//边界解第j维决策向量排序的信息
                Popx[Popx.length-1].n_var[j]=2*(Popx[Popx.length-1].dec[j]-Popx[Popx.length-2].dec[j])/(f_max-f_min);
                Popx[0].n_var[j]=2*(Popx[1].dec[j]-Popx[0].dec[j])/(f_max-f_min);
            }
            for (int i = 1; i < Popx.length-1; i++) {//非边界解第j维决策向量排序的信息``````
                double next_obj=Popx[i+1].dec[j];
                double previous_obj=Popx[i-1].dec[j];
                if(f_max-f_min==0)
                    Popx[i].n_var[j]=1;
                else
                    Popx[i].n_var[j]=(next_obj-previous_obj)/(f_max-f_min);
            }
        }
        //计算个体在决策空间的拥挤度
        double distance_sum=0;
        for (int i=0;i<Popx.length;i++){
            for (int j=0;j<dimentionx;j++){
                Popx[i].distance=Popx[i].distance+Popx[i].n_var[j];
            }
            Popx[i].distance=Popx[i].distance/dimentionx;
            Pop[Popx[i].number].distance=Popx[i].distance;//同步回去

            distance_sum+=Popx[i].distance;//求所有个体拥挤度之和
        }
        avg_crowd_dist=distance_sum/Popx.length;//计算平均拥挤度
        //计算D值
        double[] D=new double[Popx.length];
        for (int i=0;i<Popx.length;i++){
            D[i]=avg_crowd_dist/Pop[i].distance;
        }
        return D;
    }
}
