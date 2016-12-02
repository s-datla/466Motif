package miniProject;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class motifFinder {
    private List<String> sequences;
    private int ML;

    public int startPos(double[][] probDist, String sequence){
        double[] probabilities = new double[sequence.length()-probDist.length];
        int ret = 0;
        for (int i = 0; i<sequence.length()-probDist.length; i++){
        	probabilities[i] = 1;
            for (int j = 0; j<probDist.length; j++){
                if (sequence.charAt(i+j)=='A'){
                    probabilities[i]=probabilities[i]*probDist[j][0];
                }
                else if (sequence.charAt(i+j)=='C'){
                    probabilities[i]=probabilities[i]*probDist[j][1];
                }
                else if (sequence.charAt(i+j)=='G'){
                    probabilities[i]=probabilities[i]*probDist[j][2];
                }
                else if (sequence.charAt(i+j)=='T'){
                    probabilities[i]=probabilities[i]*probDist[j][3];
                }
            }
        }
        ret = 0;
        for (int i = 0; i<probabilities.length; i++){
        	if (probabilities[i]>probabilities[ret]){
        		ret = i;
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
            inputStream = new BufferedReader(new FileReader(directory + "/sequences.fa"));
            //inputStream.readLine();
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

        Random randomizer = new Random();
        int r = 0;
        String seq;
        int pos;
        List<String> otherSequences = new ArrayList<String>();
        double[][] PWM = null;
        double[] probs;
        List<Integer> positions = new ArrayList<Integer>();
        int lengthSeq = sequences.get(0).length();
        //we initialize randomly each position
        for (int i = 0; i<sequences.size() ; i++) {
            positions.add(randomizer.nextInt(lengthSeq - ML));
        }

        //for now the condition is a number of iterations, but we can change it for a duration when our algorithm works
        int K = 89000000;
        for (int k=0; k<K; k++) {
            otherSequences = sequences;
            r = randomizer.nextInt(sequences.size());
            seq = sequences.get(r);
            //I'm not sure that we need the previous position of the sequence z
            pos = positions.get(r);
            PWM = calculatePWM(sequences, positions);
            //probs = calculateP(sequences);
            pos=startPos(PWM,seq);
            positions.set(r, pos);
            if (k%500000==0){
            	PWM = calculatePWM(sequences,positions);
            	System.out.println(k);
                for  (int i = 0; i<PWM.length; i++){
                	System.out.println(""+PWM[i][0]+" "+PWM[i][1]+" "+PWM[i][2]+" "+PWM[i][3]);
                	
                }
                for (int i = 0; i<10; i++){
                	System.out.println(sequences.get(i).substring(positions.get(i), positions.get(i)+8));
                }
                System.out.println("");
                System.out.println("");
            }
            /* HERE WE CAN CALCULATE THE SCORE AND UPDATE THE POSITION FOR THE PARTICULAR SEQUENCE*/


        }
        PWM = calculatePWM(sequences,positions);
        for  (int i = 0; i<PWM.length; i++){
        	
        	System.out.println(""+PWM[i][0]+" "+PWM[i][1]+" "+PWM[i][2]+" "+PWM[i][3]);
        	System.out.println(positions.get(i));
        	System.out.println(sequences.get(i));
        	System.out.println(sequences.get(i).substring(positions.get(i), positions.get(i)+8));
        }
        
    }

    //calculate the probability of generating x according to the background model (=~1/4 for each nucleotide)
    private double[] calculateP (List<String> sequences) {
        double[] probabilities = new double[4];
        int N = 0;
        for (String seq : sequences) {
            for (int i = 0; i<seq.length(); i++) {
                if (seq.charAt(i)=='A'){
                    probabilities[0] ++;
                }
                else if (seq.charAt(i)=='C'){
                    probabilities[1] ++;
                }
                else if (seq.charAt(i)=='G'){
                    probabilities[2] ++;
                }
                else if (seq.charAt(i)=='T'){
                    probabilities[3] ++;
                }
                N++;
            }
        }
        for (int i = 0; i<4; i++) {
            probabilities[i] /= N;
        }
        return probabilities;
    }

    //calculate the position weight matrix
    private double[][] calculatePWM(List<String> sequences, List<Integer> positions) {
        double[][] PWM = new double[ML][4];
        for (int i = 0; i<sequences.size(); i++) {
            String seq = sequences.get(i).substring(positions.get(i), positions.get(i)+ML);
            for (int j = 0; j<ML; j++){
                int pos = j;
                if (seq.charAt(pos)=='A'){
                    PWM[j][0] ++;
                }
                else if (seq.charAt(pos)=='C'){
                    PWM[j][1] ++;
                }
                else if (seq.charAt(pos)=='G'){
                    PWM[j][2] ++;
                }
                else if (seq.charAt(pos)=='T'){
                    PWM[j][3] ++;
                }
            }
        }
        //pseudo-count of 0.25 + transformation in probability 0<PWM<1
        for (int j = 0; j<ML; j++){
            for (int k = 0; k<4; k++) {
                PWM[j][k] = (PWM[j][k] + 1);// / (sequences.size()+0.4) ;
            }
        }
        return PWM;
    }

    public motifFinder() {
    }


}
