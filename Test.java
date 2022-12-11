import fdv.*;
import function.MyFunction;
import mytools.TrueParetoFront;
import nsgaii.*;
import spea2.MyAlgorithm;
import tadf.*;

import java.io.*;

public class Test {
    int genNum;
    int popSize;
    int objNum;
    int constraintNum;
    int varNum;
    String testFunName;
    String testAlgoName;

    public Test(int pS, int gN, int oN, int cN, int vN, String tFN, String tAN){
        this.popSize = pS;
        this.genNum = gN;
        this.objNum = oN;
        this.constraintNum = cN;
        this.varNum = vN;
        this.testFunName = tFN;
        this.testAlgoName = tAN;
    }

    public double[] runTest(MyFunction func) throws IOException {
        double MinValue_X[] = new double[this.varNum];  // 决策变量的下界
        double MaxValue_X[] = new double[this.varNum];  // 决策变量的上界
        // 设置决策变量上下界
        for (int i = 0; i < this.varNum; i++) {
            MinValue_X[i] = 0.0;
            MaxValue_X[i] = 1.0;
        }

        TrueParetoFront TPF = new TrueParetoFront();
        int samplePointNum = varNum;
        double[][] opts = TPF.getTruePF(this.testFunName, samplePointNum, objNum);

        switch (this.testAlgoName){
            case "FDV":
                long startTime = System.currentTimeMillis();
                FUZZY algoFuzzy = new FUZZY(this.testFunName, this.popSize, this.genNum, this.varNum, this.objNum, MinValue_X, MaxValue_X);
                algoFuzzy.run();
                long endTime = System.currentTimeMillis();
                long runTime = endTime - startTime;
                double igd = algoFuzzy.calIGD(opts);
                try {
                    File file = new File("E:\\algorithmTest\\res\\fdvIGDrunTime_"+ testFunName +"obj"+ objNum + "var"+ varNum +".txt");
                    if(!file.exists()) {
                        file.createNewFile();
                    }
                    FileOutputStream fos = new FileOutputStream(file,true);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);
                    BufferedWriter bw = new BufferedWriter(osw);
                    bw.write( igd+ "\t"+ runTime);
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
                System.out.println("IGD:" + igd);
                double[] res = {igd,runTime};
                return res;
            case "NSGAII":
                long startTimeNSGAII = System.currentTimeMillis();
                NSGAII GA = new NSGAII(this.testFunName, this.popSize, this.genNum, this.varNum, this.objNum, MinValue_X, MaxValue_X);
                GA.run();
                long endTimeNSGAII = System.currentTimeMillis();
                long runTimeNSGAII = endTimeNSGAII - startTimeNSGAII;
                double igdNSGAII = GA.calIGD(opts);
                try {
                    File file = new File("E:\\algorithmTest\\res\\nsgaiiIGDrunTime_"+ testFunName +"obj"+ objNum + "var"+ varNum +".txt");
                    if(!file.exists()) {
                        file.createNewFile();
                    }
                    FileOutputStream fos = new FileOutputStream(file,true);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);
                    BufferedWriter bw = new BufferedWriter(osw);
                    bw.write( igdNSGAII+ "\t"+ runTimeNSGAII);
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
                System.out.println("IGD:" + igdNSGAII);
                double[] resNSGAII = {igdNSGAII,runTimeNSGAII};
                return resNSGAII;
            case "SPEA2":
                int funInput = 0;
                if(testFunName == "ZDT1"){
                    funInput = 1;
                }else if(testFunName == "ZDT2"){
                    funInput = 2;
                }else if(testFunName == "ZDT3"){
                    funInput = 3;
                }
                else if(testFunName == "DTLZ1"){
                    funInput = 4;
                }
                else if(testFunName == "DTLZ3"){
                    funInput = 5;
                }else if(testFunName == "DTLZ7"){
                    funInput = 6;
                }
                long startTimeSPEA2 = System.currentTimeMillis();
                MyAlgorithm ea = new spea2.MyAlgorithm(funInput,popSize,genNum,varNum,objNum);
                long endTimeSPEA2 = System.currentTimeMillis();
                long runTimeSPEA2 = endTimeSPEA2 - startTimeSPEA2;
                double igdSPEA2 = ea.calIGD(opts);
                try {
                    File file = new File("E:\\algorithmTest\\res\\spea2IGDrunTime_"+ testFunName +"obj"+ objNum + "var"+ varNum +".txt");
                    if(!file.exists()) {
                        file.createNewFile();
                    }
                    FileOutputStream fos = new FileOutputStream(file,true);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);
                    BufferedWriter bw = new BufferedWriter(osw);
                    bw.write( igdSPEA2+ "\t"+ runTimeSPEA2);
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
                System.out.println("IGD:" + igdSPEA2);
                double[] resSPEA2 = {igdSPEA2, runTimeSPEA2};
                return resSPEA2;
            case "MOEAD":
                int funInputMOEAD = 0;
                if(testFunName == "ZDT1"){
                    funInputMOEAD = 1;
                }else if(testFunName == "ZDT2"){
                    funInputMOEAD = 2;
                }else if(testFunName == "ZDT3"){
                    funInputMOEAD = 3;
                }
                else if(testFunName == "DTLZ1"){
                    funInputMOEAD = 4;
                }
                else if(testFunName == "DTLZ3"){
                    funInputMOEAD = 5;
                }else if(testFunName == "DTLZ7"){
                    funInputMOEAD = 6;
                }
                long startTimeMOEAD = System.currentTimeMillis();
                moead.MyAlgorithm moeadIns = new moead.MyAlgorithm(func,funInputMOEAD,popSize,genNum,varNum,objNum);
                long endTimeMOEAD = System.currentTimeMillis();
                long runTimeMOEAD = endTimeMOEAD - startTimeMOEAD;
                double igdMOEAD = moeadIns.calIGD(opts);
                try {
                    File file = new File("E:\\algorithmTest\\res\\moeadIGDrunTime_"+ testFunName +"obj"+ objNum + "var"+ varNum +".txt");
                    if(!file.exists()) {
                        file.createNewFile();
                    }
                    FileOutputStream fos = new FileOutputStream(file,true);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);
                    BufferedWriter bw = new BufferedWriter(osw);
                    bw.write( igdMOEAD+ "\t"+ runTimeMOEAD);
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
                System.out.println("IGD:" + igdMOEAD);
                double[] resMOEAD = {igdMOEAD, runTimeMOEAD};
                return resMOEAD;
            case "TADF":
                long st = System.currentTimeMillis();
                populations p = new populations(testFunName,objNum,varNum);
                p.run();
                long et = System.currentTimeMillis();
                long rt = et - st;
                double igdTADF = p.calIGD(opts);
                try {
                    File file = new File("E:\\algorithmTest\\res\\tadfIGDrunTime_"+ testFunName +"obj"+ objNum + "var"+ varNum +".txt");
                    if(!file.exists()) {
                        file.createNewFile();
                    }
                    FileOutputStream fos = new FileOutputStream(file,true);
                    OutputStreamWriter osw = new OutputStreamWriter(fos);
                    BufferedWriter bw = new BufferedWriter(osw);
                    bw.write( igdTADF+ "\t"+ rt);
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
                System.out.println("IGD:" + igdTADF);
                double[] resTADF = {igdTADF, rt};
                return resTADF;
        }
        double[] res = {0 , 0};
        return res;
    }
}
