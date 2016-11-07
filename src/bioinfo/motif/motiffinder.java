package bioinfo.motif;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MotifFinder {
    private List<String> sequences;
    private int ML;

    public int startPos(double[][] probDist, String sequence){
        double[] probabilities = new double[sequence.length()-probDist.length];
        int ret = 0;
        for (int i = 0; i<=sequence.length()-probDist.length; i++){
            probabilities[i] = 1;
            for (int j = i; j<i+sequence.length(); j++){if (sequence.charAt(i)=='A'){
                if (sequence.charAt(i)=='A'){
                    probabilities[i]=probabilities[i]*probDist[j][0];
                }
                else if (sequence.charAt(i)=='C'){
                    probabilities[i]=probabilities[i]*probDist[j][1];
                }
                else if (sequence.charAt(i)=='G'){
                    probabilities[i]=probabilities[i]*probDist[j][2];
                }
                else if (sequence.charAt(i)=='T'){
                    probabilities[i]=probabilities[i]*probDist[j][3];
                }
            }

            }
        }
        ret = getRandomProbability(probabilities);
        return ret;
    }

    private int getRandomProbability(double[] probabilities) {
        double sum = 0;
        for (int i = 0; i<probabilities.length; i++){
            sum=sum+probabilities[i];
        }
        int ret = 0;
        double random = Math.random()*sum;
        sum=0;
        for (int i = 0; i<probabilities.length; i++){
            sum=sum+probabilities[i];
            if (sum>=random){
                ret = i;
                i=probabilities.length;
            }
        }
        return ret;
    }


    private void readSequences(String directory) throws IOException {
        BufferedReader inputStream = null;
        String sequence = null;
        List<String> sequences = new ArrayList<String>();
        try {
            inputStream.readLine();
            inputStream = new BufferedReader(new FileReader(directory + "/sequences.fa"));
            while ((sequence = inputStream.readLine()) != null) {
                sequences.add(sequence);
            }

        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
        }
        this.sequences = sequences;
    }

    private void readMotifLength(String directory) throws IOException{
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

    public void run(){
        String dir = "benchmark/defaultParameter&num=0";
        try {
            readSequences(dir);
            readMotifLength(dir);
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("hi");
    }

    public MotifFinder() {
    }
}
