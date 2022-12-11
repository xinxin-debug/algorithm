import function.MyFunction;

import java.io.IOException;
import java.util.Arrays;

public class testMain {
    public static void main(String[] args) throws IOException {
        //
        int runTime = 31;
        int populationSize = 100;
        int generationNum = 200;
        int objNum = 3;
        int constraintNum = 3;
        int[] varNumList = {10,30,50,100};
        String testFunName ="DTLZ1";     // ZDT1 ZDT2 ZDT3 DTLZ1 DTLZ3 DTLZ7
        String testAlgoName = "FDV";  //FDV  NSGAII  SPEA2 MOEAD TADF
        double[][] resAll = new double[runTime][2];
        double[] meanIGD = new double[varNumList.length];
        double[] meanRunTime = new double[varNumList.length];

        //
        MyFunction func = new MyFunction() {
            @Override
            public double[] objectives(double[] X) {

                // 这里写用户自定义的目标函数 （以ZDT1为例）
                double sum, g;
                double[] f = new double[moead.GlobData.nfunc];

                f[0] = X[0];
                sum = 0.0;

                for (int j = 1; j < moead.GlobData.nvar; j++) {
                    sum += X[j];
                }
                g = 1 + 9 * (sum / (moead.GlobData.nvar - 1));
                f[1] = g * (1 - (Math.pow(f[0] / g, 0.5)));

                return f;

            }

            @Override
            public double[] constraints(double[] X) {
                return new double[0];
            }
        };


        for (int j = 0; j < varNumList.length;j++){
            int varNum = varNumList[j];
            double sumIgd = 0;
            double sumRunTime = 0;
            for (int i = 0; i < runTime ; i++){
                Test tt = new Test(populationSize, generationNum, objNum, constraintNum, varNum, testFunName, testAlgoName);
                resAll[i] = tt.runTest(func);
                sumIgd += resAll[i][0];
                sumRunTime += resAll[i][1];
                System.out.println("varNum: " + varNum + ", runTime: " + i + "\n");
            }
            meanIGD[j] = sumIgd/runTime;
            meanRunTime[j] = sumRunTime/runTime;
        }
        System.out.println(Arrays.toString(meanIGD));
        System.out.println(Arrays.toString(meanRunTime));




    }
}
