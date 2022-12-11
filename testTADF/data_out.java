package testTADF;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class data_out {

    public static void dataout(individual[] pop, int N, int itr){
        double [][]data_obj = new double[N][2];
        double [][]data_dec = new double[N][2];
        for (int i=0;i<N;i++){//每个个体的结果
            //System.out.println(pop[i].obj[0]+" "+pop[i].obj[1]);
            data_obj[i][0]=pop[i].obj[0];
            data_obj[i][1]=pop[i].obj[1];
            data_dec[i][0]=pop[i].dec[0];
            data_dec[i][1]=pop[i].dec[1];
        }
        new File(".\\data").mkdir();

        //目标空间数据
        try {
            File writeName = new File(".\\data\\data_obj_"+itr+".txt"); // 相对路径，如果没有则要建立一个新的output.txt文件
            if(!writeName.exists()) {
                writeName.createNewFile(); // 创建新文件,有同名的文件的话直接覆盖
            }
            FileWriter writer = new FileWriter(writeName);
            BufferedWriter out = new BufferedWriter(writer);

            //二维数组按行存入到文件中
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < 2; j++) {
                    //将每个元素转换为字符串
                    String content = String.valueOf(data_obj[i][j]) + "";
                    out.write(content + "\t");
                }
                out.write("\r\n");
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        //jc空间数据
        try {
            File writeName = new File(".\\data\\data_dec_"+itr+".txt"); // 相对路径，如果没有则要建立一个新的output.txt文件
            if(!writeName.exists()) {
                writeName.createNewFile(); // 创建新文件,有同名的文件的话直接覆盖
            }
            FileWriter writer = new FileWriter(writeName);
            BufferedWriter out = new BufferedWriter(writer);

            //二维数组按行存入到文件中
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < 2; j++) {
                    //将每个元素转换为字符串
                    String content = String.valueOf(data_dec[i][j]) + "";
                    out.write(content + "\t");
                }
                out.write("\r\n");
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
