package bioinfo.motif;

import java.io.File;

public class Main {
    private int defaultIcpc = 2;
    private int defaultMl = 8;
    private int defaultSl = 500;
    private int defaultSC = 10;

    public static void main(String[] args) {
        if (args.length != 0 && args[0].equals("benchmark"))
            new Main().run();
        else {
            String[] datasetnames = {"defaultParameter", "icpc=1.0", "icpc=1.5", "ml=6", "ml=7", "sc=5", "sc=20"};
            for (int dsnum = 0; dsnum < datasetnames.length; dsnum++) {
                for (int num = 7; num < 10; num++) {
                    String dir = "src/benchmark/" + datasetnames[dsnum] + "&num=" + num;
                    new Evaluation().doEvalutation(dir);
                }
            }
        }
    }

    void run() {
        String directory = "benchmark";
        Benchmark benchmark = new Benchmark();

        double[] icpc = {1, 1.5};
        int[] ml = {6, 7};
        int[] sc = {5, 20};

        String dataSet;
        File dir;

        for (int j = 0; j < 10; j++) {
            dataSet = directory + "/defaultParameter" + "&num=" + j;
            dir = new File(dataSet);
            dir.mkdirs();
            benchmark.createDataSet(defaultIcpc, defaultMl, defaultSl, defaultSC, dataSet);
        }

        for (int i = 0; i < icpc.length; i++) {
            for (int j = 0; j < 10; j++) {
                dataSet = directory + "/icpc=" + icpc[i] + "&num=" + j;
                dir = new File(dataSet);
                dir.mkdirs();
                benchmark.createDataSet(icpc[i], defaultMl, defaultSl, defaultSC, dataSet);
            }
        }

        for (int i = 0; i < ml.length; i++) {
            for (int j = 0; j < 10; j++) {
                dataSet = directory + "/ml=" + ml[i] + "&num=" + j;
                dir = new File(dataSet);
                dir.mkdirs();
                benchmark.createDataSet(defaultIcpc, ml[i], defaultSl, defaultSC, dataSet);
            }
        }

        for (int i = 0; i < sc.length; i++) {
            for (int j = 0; j < 10; j++) {
                dataSet = directory + "/sc=" + sc[i] + "&num=" + j;
                dir = new File(dataSet);
                dir.mkdirs();
                benchmark.createDataSet(defaultIcpc, defaultMl, defaultSl, sc[i], dataSet);
            }
        }


    }

}
