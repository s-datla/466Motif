package miniProject;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class motifFinder {
	private static double pseudocount = 0.25;
    private List<String> sequences;
    private int ML;
    private List<Integer> bestPositions = new ArrayList<Integer>();
    private double bestscore=0;
    int bestK;
    private double[] background;
    private double getScore(List<String> sequences, double[][] pwm, List<Integer> positions){
    	double score=0;
		double Aback = background[0];
		double Cback = background[1];
		double Gback = background[2];
		double Tback = background[3];
    	for (int i = 0; i<ML; i++){
    		double Aprob = pwm[i][0];
    		double Cprob = pwm[i][1];
    		double Gprob = pwm[i][2];
    		double Tprob = pwm[i][3];
    		for (int j = 0; j<sequences.size(); j++){
    			char c = sequences.get(j).charAt(positions.get(j)+i);
    			if (c=='A'){
    				score = score+(Aprob*Math.log(Aprob/Aback));
    			}
    			if (c=='C'){
    				score = score+(Cprob*Math.log(Cprob/Cback));
    			}
    			if (c=='G'){
    				score = score+(Gprob*Math.log(Gprob/Gback));
    			}
    			if (c=='T'){
    				score = score+(Tprob*Math.log(Tprob/Tback));
    			}
    		}
    	}
    	return score;
    }
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
        String dir = "benchmark/icpc=1.0&num=0";
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
        List<Integer> positions = new ArrayList<Integer>();
        int lengthSeq = sequences.get(0).length();
        //we initialize randomly each position
        for (int i = 0; i<sequences.size() ; i++) {
            positions.add(randomizer.nextInt(lengthSeq - ML));
        }

        //for now the condition is a number of iterations, but we can change it for a duration when our algorithm works
        int K = 1000000;
        background = this.calculateP(sequences);
        for (int k=0; k<K; k++) {
            otherSequences = sequences;
            r = randomizer.nextInt(sequences.size());
            seq = sequences.get(r);
            //I'm not sure that we need the previous position of the sequence z
            pos = positions.get(r);
            PWM = calculatePWM(sequences, positions,r);
            //probs = calculateP(sequences);
            pos=startPos(PWM,seq);
            positions.set(r, pos);
            double[][] tempPWM = calculatePWM(sequences,positions,-1);
            double score = this.getScore(sequences, tempPWM, positions);
            if (score>bestscore){
            	bestscore = score;
            	bestPositions.clear();
            	for (int i = 0; i<positions.size(); i++){
            		bestPositions.add(positions.get(i));
            		bestK = k;
            	}
            	
            }
            if (k%100000==0){
            	
            	System.out.println(k);
                for  (int i = 0; i<tempPWM.length; i++){
                	System.out.println(""+tempPWM[i][0]+" "+tempPWM[i][1]+" "+tempPWM[i][2]+" "+tempPWM[i][3]);
                	
                }
                for (int i = 0; i<sequences.size(); i++){
                	System.out.println(sequences.get(i).substring(positions.get(i), positions.get(i)+8));
                }
                
                System.out.println(this.getScore(sequences, PWM, positions));
                positions = phaseShiftCheck(positions,sequences,PWM);
                System.out.println(score);
                System.out.println("");
                double[][] tempPWM2 = calculatePWM(sequences,bestPositions,-1);
                double score2 = this.getScore(sequences, tempPWM2, bestPositions);
                if (score2>bestscore){
                	bestscore = score;
                	bestPositions.clear();
                	for (int i = 0; i<positions.size(); i++){
                		bestPositions.add(positions.get(i));
                		bestK = k;
                	}
                }
                System.out.println("bscore " +bestscore);
                System.out.println("bestk "+bestK);
                System.out.println("");
                if (checkMaxScore(sequences)){
                	k=K;
                }
            }
            /* HERE WE CAN CALCULATE THE SCORE AND UPDATE THE POSITION FOR THE PARTICULAR SEQUENCE*/


        }
        //PRINT EVALUATION
        PWM = calculatePWM(sequences,bestPositions,-1);
        System.out.println("alg stop best score");
    	System.out.println(bestK);
    	System.out.println(bestscore);
    	for (int i = 0; i<sequences.size(); i++){
    		System.out.println(sequences.get(i).substring(bestPositions.get(i), bestPositions.get(i)+8));
    	}
        //for  (int i = 0; i<PWM.length; i++){
        	
        	//System.out.println(""+PWM[i][0]+" "+PWM[i][1]+" "+PWM[i][2]+" "+PWM[i][3]);
        	//System.out.println(positions.get(i));
        	//System.out.println(sequences.get(i));
        	//System.out.println(sequences.get(i).substring(positions.get(i), positions.get(i)+8));
        //}
        
    }
    private boolean checkMaxScore(List<String> sequences){
    	for (int i = 0; i<sequences.size()-1; i++){
    		int p1 = bestPositions.get(i);
    		int p2 = bestPositions.get(i+1);
    		String s1 = sequences.get(i).substring(p1,p1+8);
    		String s2 = sequences.get(i+1).substring(p2, p2+8);
    		if (!s1.equals(s2)){
    			return false;
    		}
    	}
    	return true;
    }
    private List<Integer> phaseShiftCheck(List<Integer> positions, List<String> sequences, double[][] pwm){
    	double score = this.getScore(sequences, pwm, positions);
    	double leftscore = 0;
    	double rightscore = 0;
    	boolean canMoveRight = true;
    	boolean canMoveLeft = true;
    	for (int i = 0; i<positions.size(); i++){
    		if (positions.get(i)-1<0){
    			canMoveLeft = false;
    		}
    		if (positions.get(i)+1>=sequences.get(0).length()){
    			canMoveRight = false;
    		}
    	}
    	if (canMoveLeft){
    		for (int i = 0; i<positions.size(); i++){
    			positions.set(i, positions.get(i)-1);
    		}
    		pwm = calculatePWM(sequences,positions,-1);
    		leftscore = getScore(sequences,pwm,positions);
    		for (int i = 0; i<positions.size(); i++){
    			positions.set(i, positions.get(i)+1);
    		}
    	}
    	if (canMoveRight){
    		for (int i = 0; i<positions.size(); i++){
    			positions.set(i, positions.get(i)+1);
    		}
       		pwm = calculatePWM(sequences,positions,-1);
    		rightscore = getScore(sequences,pwm,positions);
    		for (int i = 0; i<positions.size(); i++){
    			positions.set(i, positions.get(i)-1);
    		}
    	}
    	if (leftscore>score && leftscore>rightscore){
    		for (int i = 0; i<positions.size(); i++){
    			positions.set(i, positions.get(i)-1);
    		}
    		return positions;
    	}
    	if (rightscore>score){
    		for (int i = 0; i<positions.size(); i++){
    			positions.set(i, positions.get(i)+1);
    		}
    	}
    	return positions;
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
    private double[][] calculatePWM(List<String> sequences, List<Integer> positions, int sequence) {
        double[][] PWM = new double[ML][4];
        for (int i = 0; i<sequences.size(); i++) {
        	if (i==sequence){
        		continue;
        	}
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
            	int size = sequences.size();
            	if (sequence>=0){
            		size--;
            	}
                PWM[j][k] = (PWM[j][k] + pseudocount)/(size+pseudocount*4);
            }
        }
        return PWM;
    }

    public motifFinder() {
    }


}
