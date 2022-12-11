package mytools;

public class TrueParetoFront {
    public double[][] getTruePF(String funName, int N, int M){
        UniformPoint UP = new UniformPoint();
        double[][] P = UP.UniformPoint(N,M);
        switch (funName){
            case "ZDT1":
                for(int i = 0; i < N; i++){
                    P[i][0] = (double) i/(N-1);
                    P[i][1] = 1 - Math.sqrt(P[i][0]);
                }
                break;
            case "ZDT2":
                for(int i = 0; i < N; i++){
                    P[i][0] = (double) 1/(N-1)*i;
                    P[i][1] = 1 - Math.pow(P[i][0],2);
                }
                break;

            case "ZDT3":
                for(int i = 0; i < N; i++){
                    P[i][0] = (double)1/(N-1)*i;
                    P[i][1] = 1 - Math.sqrt(P[i][0]) - P[i][0] * Math.sin(10 * Math.PI * P[i][0]);
                }
                break;
            case "DTLZ1":
                int samplePointNumDTLZ1 = P.length;
                for(int i = 0; i < samplePointNumDTLZ1; i++){
                    for (int j = 0; j< M; j++){
                        P[i][j] *= 0.5;
                    }
                }
                return P;
            case "DTLZ3":
                int samplePointNum = P.length;
                double[] sum = new double[samplePointNum];
                for (int i = 0; i < samplePointNum; i++){
                    sum[i] = 0;
                    for (int j = 0; j < M; j++){
                        sum[i] = sum[i] + P[i][j] * P[i][j];
                    }
                }
                for (int i = 0 ; i < samplePointNum; i++){
                    for (int j = 0; j < M; j++){
                        P[i][j] = P[i][j] / Math.sqrt(sum[i]);
                    }
                }
                break;
            case "DTLZ7":
                double[] interval = new double[]{0,0.251412,0.631627,0.859401};
                double median = (double)(interval[1] - interval[0]) / (interval[3] - interval[2] + interval[1] - interval[0]);
                double[][] X = new double[N][M-1];
                if(M == 2){
                    for(int k = 0; k < N; k++){
                        X[k][0] = (double) k / N - 1;
                    }
                }else
                    X = UP.UniformPoint(N, M-1);
                for (int i = 0; i < N; i++){
                    double sumX = 0;
                    for (int j = 0; j < M-1; j++){
                        if(X[i][j] <= median){
                            X[i][j] = X[i][j] * (interval[1] - interval[0])/median + interval[0];
                        }else {
                            X[i][j] = X[i][j] * (interval[3] - interval[2])/(1 - median) + interval[3];
                        }
                        sumX += X[i][j] /2*(1+Math.sin(3*Math.PI*X[i][j]));
                        P[i][j] = X[i][j];
                    }
                    P[i][M-1] = 2 * (M - sumX);
                    break;
                }
        }
        return P;
    }
}
