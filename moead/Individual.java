package moead;

// 个体类
public class Individual {
    double[] xreal;     //  决策变量
    double[] fitness;    // 函数值适应度

    // 构造方法
    public Individual() {
        xreal = new double[GlobData.nvar];     // 决策变量维度
        fitness = new double[GlobData.nfunc];  // 目标函数个数
    }

    // 计算个体的函数值
    /*
    public void indFunc(int fun_input) {
        double sum, g;
        double[] x = new double[GlobData.nvar];
        double[] f = new double[GlobData.nfunc];

        System.arraycopy(xreal, 0, x, 0, GlobData.nvar);

        // 计算函数值
        if (fun_input == 1) {       // ZDT1
            f[0] = x[0];
            sum = 0.0;

            for (int j = 1; j < moead.GlobData.nvar; j++) {
                sum += x[j];
            }
            g = 1 + 9 * (sum / (moead.GlobData.nvar - 1));
            f[1] = g * (1 - (Math.pow(f[0] / g, 0.5)));
        } else if(fun_input == 2){ //ZDT2
            f[0] = x[0];
            sum = 0.0;

            for (int j = 1; j < moead.GlobData.nvar; j++) {
                sum += x[j];
            }
            g = 1 + 9 * (sum / (moead.GlobData.nvar - 1));
            f[1] = g * (1 - (Math.pow(f[0] / g, 2)));
        }else if (fun_input == 3){     //ZDT3
            f[0] = x[0];
            sum = 0.0;

            for (int j = 1; j < moead.GlobData.nvar; j++) {
                sum += x[j];
            }
            g = 1 + 9 * (sum / (moead.GlobData.nvar - 1));
            f[1] = g *(1 -Math.pow((f[0]/g),0.5)-(f[0]/g)*Math.sin((10*Math.PI*f[0])));
        }else if (fun_input == 4) {    // DTLZ1
            sum = 0.0;
            for (int j = moead.GlobData.nfunc; j < moead.GlobData.nvar; j++) {
                sum += (Math.pow(x[j] - 0.5, 2) - Math.cos(20 * Math.PI * (x[j] - 0.5)));
            }

            double g1 = 100.0 * (moead.GlobData.nvar - moead.GlobData.nfunc) + 100 * sum;
            if(moead.GlobData.nfunc <= 2){
                f[0] = 0.5 * (1 + g1) * x[0];
                f[1] = 0.5 * (1 + g1) * (1 - x[0]);
            }else{
                f[0] = 1;
                for(int i1 = 0; i1 < moead.GlobData.nfunc - 1; i1++){
                    f[0] = f[0] * x[i1];
                }
                f[0] = 0.5 * (1 + g1) * f[0];
                for(int i1 = 1; i1 < moead.GlobData.nfunc; i1++){
                    f[i1] = 1;
                    for(int j = 0; j < moead.GlobData.nfunc - i1 - 1; j++) {
                        f[i1] = f[i1] * x[j];
                    }
                    f[i1] = f[i1] * (1 - x[moead.GlobData.nfunc - i1 - 1]) * 0.5 * (1 + g1);
                }
                // Y[objNum-1] = (1 + g) * (1 - X[0]);
            }
        } else if (fun_input == 5) { //DTLZ3
            double sum2 = 0.0;
            for (int j = moead.GlobData.nfunc; j < moead.GlobData.nvar; j++) {
                sum2 += (Math.pow(x[j] - 0.5, 2) - Math.cos(20 * Math.PI * (x[j] - 0.5)));
            }

            g = 100.0 * (moead.GlobData.nvar - moead.GlobData.nfunc) + 100 * sum2;
            if(moead.GlobData.nfunc <= 2){
                f[0] = (1 + g) * Math.cos( x[0] * Math.PI/2);
                f[1] = (1 + g) * Math.sin( x[0] * Math.PI/2);
            }else{
                f[0] = 1;
                for(int i1 = 0; i1 < moead.GlobData.nfunc - 1; i1++){
                    f[0] = f[0] * Math.cos( x[i1] * Math.PI/2);
                }
                f[0] = (1 + g) * f[0];
                for(int i2 = 1; i2 < moead.GlobData.nfunc; i2++){
                    f[i2] = 1;
                    for(int j = 0; j < moead.GlobData.nfunc - i2 - 1; j++) {
                        f[i2] = f[i2] * Math.cos( x[j] * Math.PI/2);
                    }
                    f[i2] = f[i2] * (1 + g) * Math.sin( x[moead.GlobData.nfunc - i2 - 1] * Math.PI/2);
                }
                // Y[objNum-1] = (1 + g) * (1 - X[0]);
            }
        }else if (fun_input == 6){   //DTLZ7
            double sum7 = 0;
            for (int j = moead.GlobData.nfunc; j < moead.GlobData.nvar; j++) {
                sum7 += x[j];
            }
            double g7 = 1 + 9 / (moead.GlobData.nvar - moead.GlobData.nfunc) * sum7;
            double sumH = 0;
            for(int i3 = 0; i3 < moead.GlobData.nfunc - 1; i3++){
                f[i3] = x[i3];
                sumH += f[i3]/(1+g7)*(1+Math.sin(3*Math.PI*f[i3]));
            }
            f[moead.GlobData.nfunc - 1] = (1 + g7) * (moead.GlobData.nfunc-sumH);
        }

        // 更新个体的函数值
        System.arraycopy(f, 0, fitness, 0, moead.GlobData.nfunc);
    }
*/

}
