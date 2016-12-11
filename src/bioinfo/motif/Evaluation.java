package bioinfo.motif;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class Evaluation {
    private int ML;
    private double[][] motif;
    private double[][] predictedMotif;
    private List<Integer> sites = new ArrayList<Integer>();
    private List<Integer> predictedSites = new ArrayList<Integer>();


    public void doEvalutation(String dir) {
        try {
            readMotifLength(dir);
            readSites(dir);
            readMotif(dir);
            readPredictedSites(dir);
            readPredictedMotif(dir);

        } catch (IOException e) {
            e.printStackTrace();
        }


        double relativeEntropy =  relativeEntropy(motif, predictedMotif);
        double overlappingPositions = overlappingPositions(sites, predictedSites);
        System.out.println("dir: " + dir);
        System.out.println("relative entropy:" +relativeEntropy);
        System.out.println("ratio overlapping positions: " + overlappingPositions);

        printEvaluation(dir, relativeEntropy, overlappingPositions);
    }

    //more relevant as a ratio
    private double overlappingPositions(List<Integer> sites, List<Integer> predictedSites){
        double nbOverlappingPositions = 0;
        for (int i=0; i<sites.size(); i++){
            int diff = Math.abs(sites.get(i)- predictedSites.get(i));
            if (diff < ML){
                nbOverlappingPositions += ML - diff;
            }
        }
        return nbOverlappingPositions/(ML*sites.size());
    }



    private double relativeEntropy(double[][] motif, double[][] predictedMotif ){
        double entropy = 0;
        for (int i=0; i<ML; i++){
            //log2??
            entropy += calculateEntropy(motif[i], predictedMotif[i]);
        }
        return entropy;
    }

    private double calculateEntropy(double[] nucleotideDistribution, double[] predictedNucleotideDistribution){
        double entropyPosition = 0;
        for (int i=0; i<4; i++){
            if (nucleotideDistribution[i] == 0 && predictedNucleotideDistribution[i] !=0){
                System.out.println("ERROR we can't predict the relative entropy if the probability of a nucleotid is 0 for the motif and different of 0 for the predicted motif!!!!");
                entropyPosition = 1000000000;
            }
            else{
                if (predictedNucleotideDistribution[i] != 0) {
                    entropyPosition += predictedNucleotideDistribution[i]*Math.log(predictedNucleotideDistribution[i]/nucleotideDistribution[i]);
                }
            }
        }
        return entropyPosition;
    }

    private void readMotifLength(String directory) throws IOException {
        BufferedReader inputStream = null;
        try {
            inputStream = new BufferedReader(new FileReader(directory + "/motiflength.txt"));
            String ML = inputStream.readLine();
            this.ML = Integer.valueOf(ML);
        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
        }
    }

    private void readMotif(String directory) throws IOException {
        BufferedReader inputStream = null;
        String[] nucleotides;
        motif = new double[ML][4];
        try {
            inputStream = new BufferedReader(new FileReader(directory + "/motif.txt"));
            inputStream.readLine();
            for (int i=0; i<ML; i++){
                nucleotides = inputStream.readLine().split(" ");
                for (int j=0; j<4; j++) {
                    motif[i][j] = Double.valueOf(nucleotides[j]);
                }
            }
        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
        }
    }

    private void readSites(String directory) throws IOException {
        BufferedReader inputStream = null;
        try {
            inputStream = new BufferedReader(new FileReader(directory + "/sites.txt"));
            String site;
            while ((site = inputStream.readLine())!=null) {
                sites.add(Integer.valueOf(site));

            }
        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
        }
    }

    private void readPredictedMotif(String directory) throws IOException {
        BufferedReader inputStream = null;
        String[] nucleotides;
        predictedMotif = new double[ML][4];
        try {
            inputStream = new BufferedReader(new FileReader(directory + "/predictedmotif.txt"));
            inputStream.readLine();
            for (int i=0; i<ML; i++){
                nucleotides = inputStream.readLine().split(" ");
                for (int j=0; j<4; j++) {
                    predictedMotif[i][j] = Double.valueOf(nucleotides[j]);
                }
            }
        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
        }
    }

    private void readPredictedSites(String directory) throws IOException {
        BufferedReader inputStream = null;
        try {
            inputStream = new BufferedReader(new FileReader(directory + "/predictedsites.txt"));
            String site;
            while ((site = inputStream.readLine())!=null) {
                predictedSites.add(Integer.valueOf(site));
            }
        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
        }
    }

    void printEvaluation(String dir, double relativeEntropy, double overlaps){
        try {
            File f = new File("src/bioinfo/motif/evaluation.csv");

            if (!(f.exists() && !f.isDirectory())) {
                f.createNewFile();
            }
            FileWriter fw = new FileWriter(f, true);
            fw.write(dir+";" + relativeEntropy + ";" + overlaps);
            fw.write("\n");
            fw.flush();
            fw.close();
        } catch (Exception e){
            e.printStackTrace();
        }

    }
}
