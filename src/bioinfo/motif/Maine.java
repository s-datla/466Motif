package bioinfo.motif;

import java.io.File;

public class Maine {
    private int defaultIcpc = 2;
    private int defaultMl = 8;
    private int defaultSl = 500;
    private int defaultSC = 10;

    public static void main(String[] args) {
        new Maine().run();
    }

    void run() {
        String directory = "benchmark";
        Main benchmark = new Main();

        double[] icpc = {1, 1.5};
        int[] ml = {6, 7};
        int[] sc = {5, 20};

        String dataSet = directory + "/defaultParameter";
        benchmark.createDataSet(defaultIcpc, defaultMl, defaultSl, defaultSC, dataSet);

        for (int i = 0; i < icpc.length; i++) {
            for (int j = 0; j < 10; j++) {
                dataSet = directory + "/icpc=" + icpc[i] + "&num=" + j;
                File dir = new File(dataSet);
                dir.mkdirs();
                benchmark.createDataSet(icpc[i], defaultMl, defaultSl, defaultSC, dataSet);
            }
        }

        for (int i = 0; i < ml.length; i++) {
            for (int j = 0; j < 10; j++) {
                dataSet = directory + "/ml=" + ml[i] + "&num=" + j;
                File dir = new File(dataSet);
                dir.mkdirs();
                benchmark.createDataSet(defaultIcpc, ml[i], defaultSl, defaultSC, dataSet);
            }
        }

        for (int i = 0; i < sc.length; i++) {
            for (int j = 0; j < 10; j++) {
                dataSet = directory + "/sc=" + sc[i] + "&num=" + j;
                File dir = new File(dataSet);
                dir.mkdirs();
                benchmark.createDataSet(defaultIcpc, defaultMl, defaultSl, sc[i], dataSet);
            }
        }


    }

}
