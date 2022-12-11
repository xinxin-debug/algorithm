package fdv;

public class Fitness {
    public double[] CalFitnessZDT1(double[] X){
        double[] Y = new double[2];
        int len = X.length;
        Y[0] = X[0];
        double Sum = 0;
        for (int i = 1; i < len ; i++) {
            Sum += X[i];
        }
        double g = 1.0 + 9.0 * (Sum/(len - 1));
        Y[1] = g * (1 -Math.pow((Y[0]/g),0.5));      //ZDT1;
        //Y[1] = g * (1 -Math.pow((Y[0]/g),2));        //ZDT2;
        //Y[1] = g *(1 -Math.pow((Y[0]/g),0.5)-(Y[0]/g)*Math.sin((10*Math.PI*Y[0])));  //ZDT3
        return Y;
    }

    public double[] CalFitnessZDT2(double[] X){
        double[] Y = new double[2];
        int len = X.length;
        Y[0] = X[0];
        double Sum = 0;
        for (int i = 1; i < len ; i++) {
            Sum += X[i];
        }
        double g = 1.0 + 9.0 * (Sum/(len - 1));
        //Y[1] = g * (1 -Math.pow((Y[0]/g),0.5));      //ZDT1;
        Y[1] = g * (1 -Math.pow((Y[0]/g),2));        //ZDT2;
        //Y[1] = g *(1 -Math.pow((Y[0]/g),0.5)-(Y[0]/g)*Math.sin((10*Math.PI*Y[0])));  //ZDT3
        return Y;
    }

    public double[] CalFitnessZDT3(double[] X){
        double[] Y = new double[2];
        int len = X.length;
        Y[0] = X[0];
        double Sum = 0;
        for (int i = 1; i < len ; i++) {
            Sum += X[i];
        }
        double g = 1.0 + 9.0 * (Sum/(len - 1));
        //Y[1] = g * (1 -Math.pow((Y[0]/g),0.5));      //ZDT1;
        //Y[1] = g * (1 -Math.pow((Y[0]/g),2));        //ZDT2;
        Y[1] = g *(1 -Math.pow((Y[0]/g),0.5)-(Y[0]/g)*Math.sin((10*Math.PI*Y[0])));  //ZDT3
        return Y;
    }

    public double[] CalFitnessDTLZ1(double[] X, int objNum){
        double[] Y = new double[objNum];
        int len = X.length;
        double sum = 0.0;
        for (int j = objNum; j < len; j++) {
            sum += (Math.pow(X[j] - 0.5, 2) - Math.cos(20 * Math.PI * (X[j] - 0.5)));
        }

        double g = 100.0 * (len - objNum) + 100 * sum;
        if(objNum <= 2){
            Y[0] = 0.5 * (1 + g) * X[0];
            Y[1] = 0.5 * (1 + g) * (1 - X[0]);
            return Y;
        }else{
            Y[0] = 1;
            for(int i = 0 ; i < objNum - 1; i++){
                Y[0] = Y[0] * X[i];
            }
            Y[0] = 0.5 * (1 + g) * Y[0];
            for(int i = 1; i < objNum; i++){
                Y[i] = 1;
                for(int j = 0; j < objNum - i - 1;j++) {
                    Y[i] = Y[i] * X[j];
                }
                Y[i] = Y[i] * (1 - X[objNum - i - 1]) * 0.5 * (1 + g);
            }
            // Y[objNum-1] = (1 + g) * (1 - X[0]);
            return Y;
        }
    }

    public double[] CalFitnessDTLZ3(double[] X, int objNum){
        double[] Y = new double[objNum];
        int len = X.length;
        double sum = 0.0;
        for (int j = objNum; j < len; j++) {
            sum += (Math.pow(X[j] - 0.5, 2) - Math.cos(20 * Math.PI * (X[j] - 0.5)));
        }

        double g = 100.0 * (len - objNum) + 100 * sum;
        if(objNum <= 2){
            Y[0] = (1 + g) * Math.cos( X[0] * Math.PI/2);
            Y[1] = (1 + g) * Math.sin( X[0] * Math.PI/2);
        }else{
            Y[0] = 1;
            for(int i = 0 ; i < objNum - 1; i++){
                Y[0] = Y[0] * Math.cos( X[i] * Math.PI/2);
            }
            Y[0] = (1 + g) * Y[0];
            for(int i = 1; i < objNum; i++){
                Y[i] = 1;
                for(int j = 0; j < objNum - i - 1;j++) {
                    Y[i] = Y[i] * Math.cos( X[j] * Math.PI/2);
                }
                Y[i] = Y[i] * (1 + g) * Math.sin( X[objNum - i - 1] * Math.PI/2);
            }
            // Y[objNum-1] = (1 + g) * (1 - X[0]);
        }
        return Y;
    }

    public double[] CalFitnessDTLZ7(double[] X, int objNum){
        double[] Y = new double[objNum];
        int len = X.length;
        double sum = 0.0;
        for (int j = objNum; j < len; j++) {
            sum += X[j];
        }
        double g = 1 + 9 / (len - objNum) * sum;
        double sumH = 0;
        for(int i = 0 ; i < objNum - 1; i++){
            Y[i] = X[i];
            sumH += Y[i]/(1+g)*(1+Math.sin(3*Math.PI*Y[i]));
        }
        Y[objNum - 1] = (1 + g) * (objNum-sumH);
        return Y;
    }
}


