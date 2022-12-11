package tadf;

import java.util.Comparator;



public class MyComparator implements Comparator<individual> {
    private int index;
    private int model;

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    @Override
    public int compare(individual o1, individual o2) {
        if(index==0||index==1){//如果是根据决策空间来排序
            if (o1.dec[index]>o2.dec[index]){
                return 1;
            } else if (o1.dec[index]==o2.dec[index]){
                return 0;
            }
        }else if(index==2){//如果根据个体适应度排序
            if (o1.fitness>o2.fitness){
                return 1;
            } else if (o1.fitness==o2.fitness){
                return 0;
            }
        }
        return -1;
    }
}
