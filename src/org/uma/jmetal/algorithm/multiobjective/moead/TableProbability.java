/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.uma.jmetal.algorithm.multiobjective.moead;

/**
 *
 * @author UZIEL
 */
public class TableProbability {

    private final int methodsQuality;
    private int[] tableProbability;
    private double[][] solutionBestAndBadByMethod;
    private int[] porcetageByMethod;
    private int porcentage = 100;

    public TableProbability() {
        this(4);
    }

    public TableProbability(int tam) {
        methodsQuality = tam;
        porcetageByMethod = new int[tam];
        tableProbability = new int[porcentage];
    }

    public void initializeProbability() {
        int porcentageIndividual = (int) Math.ceil(100 / methodsQuality);
        setPorcentageForEachMethod(porcentageIndividual);
        repairPorcentageIsViolatedLimite();
        fillTableProbability();
    }

    private void setPorcentageForEachMethod(int porcetage) {
        for (int i = 0; i < methodsQuality; i++) {
            porcetageByMethod[i] = porcetage;
        }
    }

    private void repairPorcentageIsViolatedLimite() {
        int sumPorcentage = getSumPorcentageByMethod();

        if (sumPorcentage < porcentage || sumPorcentage > porcentage ) {
            if(sumPorcentage < porcentage){
                repairBottomPorcentage(sumPorcentage);     
            }else{
                repairTopPorcentage(sumPorcentage);    
            }
        } 
    }
    
    private void repairTopPorcentage( int sumPorcentage ){                
        int rest = sumPorcentage - porcentage;
        int[] serie  = generateSerie();
        
        //corregir probabilidades        
        for (int i = 0; i < serie.length; i++) {
            if (rest > 0) {
                if (porcetageByMethod[serie[i]] != 1) {
                    porcetageByMethod[serie[i]]--; 
                    rest--;
                }
            } else {
                break;
            }
        }                       
    }  
    
    private void repairBottomPorcentage( int sumPorcentage ){
        int cant = porcentage - sumPorcentage;
        int[] serie  = generateSerie();
        
         //corregir probabilidades        
        for (int i = 0; i < serie.length; i++) {
            if (cant > 0) {
                if (porcetageByMethod[serie[i]] != 1) {
                    porcetageByMethod[serie[i]]++; 
                    cant--;
                }
            } else {
                break;
            }
        }                        
    }
    
 
    public int getSumPorcentageByMethod() {
        int sum = 0;
        for (int i = 0; i < methodsQuality; i++) {
            sum += porcetageByMethod[i];
        }

        return sum;
    }
    
    private int[] generateSerie(){
        int[] serie = new int[methodsQuality * methodsQuality];
        int value = 0;
        
        for (int i = 0; i < serie.length; i++) {            
            if(value % methodsQuality == 0 ){
               value = 0;
            }                
            value++;
            serie[i] = value - 1;  
        }
        return serie;        
    }
    
    //***********

    public void setporcetageByMethod(int[] arr) {
        this.porcetageByMethod = arr;
    }

 

    public void fillTableProbability() {
        int start = 0, end = 0;

        for (int i = 1; i <= porcetageByMethod.length; i++) {
            end += porcetageByMethod[i - 1];
            for (int j = start; j < end; j++) {
                this.tableProbability[j] = i;
            }
            start = end;
        }
    }

    public void calculateProbabilityForEachMethod() {
        double e = 0.01, probabilityMethod = 0, sum;

        double sumTotal = getDoubleWithTwoDecimals(getTotalSumatorias());

        for (int i = 0; i < this.methodsQuality; i++) {
            sum = solutionBestAndBadByMethod[i][2];
            probabilityMethod = sum / sumTotal + e;

            solutionBestAndBadByMethod[i][3] = getDoubleWithTwoDecimals(probabilityMethod);
        }
    }

    public void setporcetageByMethod() {
        int sum = 0;
        for (int i = 0; i < methodsQuality; i++) {
            double num = solutionBestAndBadByMethod[i][3] * 100;
            porcetageByMethod[i] = (int) (solutionBestAndBadByMethod[i][3] * 100);
            sum += porcetageByMethod[i];
        }

        int tope = (sum - 100) * 2;

        //corregir probabilidades
        //while (tope > 0) {
        for (int i = 0; tope > i; i++) {
            if (porcetageByMethod[i] != 1) {
                porcetageByMethod[i]--;
                tope--;
            }
        }
        //}        
    }

    private double getTotalSumatorias() {
        double sum = 0;
        for (int i = 0; i < this.methodsQuality; i++) {
            sum += solutionBestAndBadByMethod[i][2];
        }
        return sum;
    }

    public void calculateSumatoria() {
        double e = 0.01;
        double SUM = 0, bestSolutions = 0, badSolutions = 0;

        for (int i = 0; i < methodsQuality; i++) {
            bestSolutions = solutionBestAndBadByMethod[i][0];
            badSolutions = solutionBestAndBadByMethod[i][1];

            SUM = bestSolutions / (bestSolutions + badSolutions + e);
            solutionBestAndBadByMethod[i][2] = getDoubleWithTwoDecimals(SUM);
        }

    }

    private double getDoubleWithTwoDecimals(double number) {
        number = Math.round(number * 100);
        number /= 100;
        return number;
    }

    public void setSolutionBestAndBadByMethod(double[][] table) {
        this.solutionBestAndBadByMethod = table;
    }

    public int[] getTableProbability() {
        return this.tableProbability;
    }
}
