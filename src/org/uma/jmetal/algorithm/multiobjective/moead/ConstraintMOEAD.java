package org.uma.jmetal.algorithm.multiobjective.moead;

import org.uma.jmetal.algorithm.multiobjective.moead.util.MOEADUtils;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.comparator.impl.ViolationThresholdComparator;

import java.util.List;

/**
 * This class implements a constrained version of the MOEAD algorithm based on
 * the one presented in the paper: "An adaptive constraint handling approach
 * embedded MOEA/D". DOI: 10.1109/CEC.2012.6252868
 *
 * @author Antonio J. Nebro
 * @author Juan J. Durillo
 * @version 1.0
 */
@SuppressWarnings("serial")
public class ConstraintMOEAD extends AbstractMOEAD<DoubleSolution> {

    private DifferentialEvolutionCrossover differentialEvolutionCrossover;
    private ViolationThresholdComparator<DoubleSolution> violationThresholdComparator;

    private double[][] solutionBestAndBadByMethod;
    //private Probability probability;
    private int[] tablePropability;

    public ConstraintMOEAD(Problem<DoubleSolution> problem,
            int populationSize,
            int resultPopulationSize,
            int maxEvaluations,
            MutationOperator<DoubleSolution> mutation,
            CrossoverOperator<DoubleSolution> crossover,
            FunctionType functionType,
            String dataDirectory,
            double neighborhoodSelectionProbability,
            int maximumNumberOfReplacedSolutions,
            int neighborSize) {
        super(problem, populationSize, resultPopulationSize, maxEvaluations, crossover, mutation, functionType,
                dataDirectory, neighborhoodSelectionProbability, maximumNumberOfReplacedSolutions,
                neighborSize);

        differentialEvolutionCrossover = (DifferentialEvolutionCrossover) crossoverOperator;
        violationThresholdComparator = new ViolationThresholdComparator<DoubleSolution>();
        solutionBestAndBadByMethod = new double[4][4];
        
       
        tablePropability = new int[100];
    }

    @Override
    public void run() {
        initializeUniformWeight();
        initializeNeighborhood();
        initializePopulation();
        idealPoint.update(population);

        violationThresholdComparator.updateThreshold(population);

        evaluations = populationSize;

        tablePropability = _initializeProbability();
                

        int LP = 0, CLP = 0;
        LP = (int) Math.round(0.5 * population.get(0).getNumberOfVariables()) + 2;

        do {
            
            
            differentialEvolutionCrossover._setProbability(tablePropability);
            
            differentialEvolutionCrossover._setSFSandSIS(population);

            int[] permutation = new int[populationSize];
            MOEADUtils.randomPermutation(permutation, populationSize);

            for (int i = 0; i < populationSize; i++) {
               
                int subProblemId = permutation[i];

                NeighborType neighborType = chooseNeighborType();
                List<DoubleSolution> parents = parentSelection(subProblemId, neighborType);

                //***
                differentialEvolutionCrossover._setPopulation(population);

                differentialEvolutionCrossover.setCurrentSolution(population.get(subProblemId));
                List<DoubleSolution> children = differentialEvolutionCrossover.execute(parents);

                DoubleSolution child = children.get(0);
                mutationOperator.execute(child);
                problem.evaluate(child);

                if (_isParticleWasRepair()) {                                                         
                    int flag = _evaluateMethod(subProblemId, child);
                    _updateSolutionsBestAndBadByMethod(flag);                                        
                }

                evaluations++;

                idealPoint.update(child.getObjectives());
                updateNeighborhood(child, subProblemId, neighborType);
            }

            CLP++;

            if (CLP % LP == 0) {
                _calculateSumatoria();
                _getProbabilityForEachMethod();
                int[] viabilityForMethod = _getViabilityByMethod();
                
                if(viabilityForMethod[0] == 0){
                    tablePropability = _initializeProbability();
                }else{
                    tablePropability = _fillTableProbability(viabilityForMethod);                       
                }
                                
                _clearSolutionBestAndBadByMethod();
            }

            violationThresholdComparator.updateThreshold(population);

        } while (evaluations < maxEvaluations);
    }

    private int[] _initializeProbability() {

        int[] viabilityByMethod = {25, 25, 25, 25};
        int[] probability = _fillTableProbability(viabilityByMethod);

        return probability;
    }

    private int[] _getViabilityByMethod() {
        int[] probability = new int[solutionBestAndBadByMethod.length];
        int sum = 0;
        
        for (int i = 0; i < solutionBestAndBadByMethod.length; i++) {            
            probability[i] = (int) (solutionBestAndBadByMethod[i][3] * 100);
            sum += probability[i];
        }

        probability = repairProbability(probability);
        
        return probability;
    }

    private int[] repairProbability(int[] arr) {
        int[] opc = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
        int[] probability = new int[arr.length];        
        System.arraycopy(arr, 0, probability, 0, arr.length);
        
        int tope = getResto(arr);
        
        //corregir probabilidades        
        for (int i = 0; i < opc.length; i++) {
            if (tope > 0) {
                if (probability[opc[i]] != 1) {
                    probability[opc[i]]--;
                    tope--;
                }
            } else {
                break;
            }
        }
       
        return probability;
    }
    
    private int getResto(int[] arr){
        int sum = 0;
        for (int i = 0; i < arr.length; i++) {
            sum += arr[i];
        }
        
        return sum % 100;
    }

    private int[] _fillTableProbability(int[] viabilityByMethod) {
        int[] probability = new int[100];
        int start = 0, end = 0;

        for (int i = 1; i <= viabilityByMethod.length; i++) {
            end += viabilityByMethod[i - 1];
            for (int j = start; j < end; j++) {
                probability[j] = i;
            }
            start = end;
        }
        return probability;
    }

    private boolean _isParticleWasRepair() {
        return differentialEvolutionCrossover._isParticleRepaired();
    }

    private void _updateSolutionsBestAndBadByMethod(int flag) {
        final int BEST_METHOD = 0, BAD_METHOD = 1;
        int appliedMethod = differentialEvolutionCrossover._getMethodUsed() - 1;
        
       
        if (flag == BEST_METHOD) {
            setBestMethod(appliedMethod);
        } else if (flag == BAD_METHOD) {
            setBadtMethod(appliedMethod);
        } else {
            setBadtMethod(appliedMethod);
        }
    }

    private void setBestMethod(int appliedMethod) {
        solutionBestAndBadByMethod[appliedMethod][0]++;
    }

    private void setBadtMethod(int appliedMethod) {
        solutionBestAndBadByMethod[appliedMethod][1]++;
    }

    private int _evaluateSolutionFeasible(int subProblemId, DoubleSolution child) {
        final int A_DOMINANCE_B = 1, B_DOMINANCE_A = -1, NO_DOMINANCE = 0;
        int resp = 0;

        int dominance = checkDominance(population.get(subProblemId), child);

        switch (dominance) {
            case A_DOMINANCE_B:
                resp = 1;
                break;
            case B_DOMINANCE_A:
                resp = 0;
                break;
            case NO_DOMINANCE: // mejoro la solucion
                resp = 0;
                break;
            default:
                resp = 0;
                break;
        }

        return resp;
    }

    private int _evaluateMethod(int subProblemId, DoubleSolution child) {
        final int EQUAL_SOLUTIONS = 0, BAD_SOLUTION = -1, IMPROVE_SOLUTION = 1;

        int method = compareIfMethodSolutionImprove(subProblemId, child);
        int resp = 0;

        switch (method) {
            case BAD_SOLUTION:
                resp = 1;
                break;
            case IMPROVE_SOLUTION:
                resp = 0;
                break;
            case EQUAL_SOLUTIONS:
                resp = _evaluateSolutionFeasible(subProblemId, child);
                break;
        }

        return resp;
    }

    private int compareIfMethodSolutionImprove(int subProblemId, DoubleSolution child) {
        return violationThresholdComparator.compare(population.get(subProblemId), child);
    }

    private void _clearSolutionBestAndBadByMethod() {
        solutionBestAndBadByMethod = new double[4][4];
    }

    private void _calculateSumatoria() {
        double e = 0.01;
        double SUM = 0, bestSolutions = 0, badSolutions = 0;

        for (int i = 0; i < solutionBestAndBadByMethod.length; i++) {
            bestSolutions = solutionBestAndBadByMethod[i][0];
            badSolutions = solutionBestAndBadByMethod[i][1];

            SUM = bestSolutions / (bestSolutions + badSolutions + e);
            solutionBestAndBadByMethod[i][2] = _getDoubleWithTwoDecimals(SUM);
        }

    }

    private double _getDoubleWithTwoDecimals(double number) {
        number = Math.round(number * 100);
        number /= 100;
        return number;
    }

    private void _getProbabilityForEachMethod() {
        double e = 0.01, probabilityMethod = 0, sum;

        double sumTotal = _getDoubleWithTwoDecimals(_getTotalSumatorias());

        for (int i = 0; i < solutionBestAndBadByMethod.length; i++) {
            sum = solutionBestAndBadByMethod[i][2];
            probabilityMethod = sum / sumTotal + e;

            solutionBestAndBadByMethod[i][3] = _getDoubleWithTwoDecimals(probabilityMethod);
        }
    }

    private double _getTotalSumatorias() {
        double sum = 0;
        for (int i = 0; i < solutionBestAndBadByMethod.length; i++) {
            sum += solutionBestAndBadByMethod[i][2];
        }
        return sum;
    }

    /**
     * check the dominance relationship between a and b: 1 -> a dominates b, -1
     * -> b dominates a 0 -> non-dominated with each other
     *
     * @param solutionB
     */
    public int checkDominance(DoubleSolution solutionA, DoubleSolution solutionB) {
        int flag1 = 0;
        int flag2 = 0;

        for (int i = 0; i < problem.getNumberOfObjectives(); i++) {
            if (solutionA_dominate_B(solutionA.getObjective(i), solutionB.getObjective(i))) {
                flag1 = 1;
            } else {
                if (solutionB_dominate_A(solutionA.getObjective(i), solutionB.getObjective(i))) {
                    flag2 = 1;
                }
            }
        }
        if (flag1 == 1 && flag2 == 0) {
            return 1;
        } else {
            if (flag1 == 0 && flag2 == 1) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    private boolean solutionA_dominate_B(double objetiveA, double objetiveB) {
        return objetiveA < objetiveB;
    }

    private boolean solutionB_dominate_A(double objetiveA, double objetiveB) {
        return objetiveA > objetiveB;
    }

    public void initializePopulation() {
        for (int i = 0; i < populationSize; i++) {
            DoubleSolution newSolution = (DoubleSolution) problem.createSolution();

            problem.evaluate(newSolution);
            population.add(newSolution);
        }
    }

    @Override
    protected void updateNeighborhood(DoubleSolution individual, int subproblemId, NeighborType neighborType) {
        int size;
        int time;

        time = 0;

        if (neighborType == NeighborType.NEIGHBOR) {
            size = neighborhood[subproblemId].length;
        } else {
            size = population.size();
        }
        int[] perm = new int[size];

        MOEADUtils.randomPermutation(perm, size);

        for (int i = 0; i < size; i++) {
            int k;
            if (neighborType == NeighborType.NEIGHBOR) {
                k = neighborhood[subproblemId][perm[i]];
            } else {
                k = perm[i];
            }
            double f1, f2;

            f1 = fitnessFunction(population.get(k), lambda[k]);
            f2 = fitnessFunction(individual, lambda[k]);

            if (violationThresholdComparator.needToCompare(population.get(k), individual)) {
                int flag = violationThresholdComparator.compare(population.get(k), individual);
                if (flag == 1) {
                    population.set(k, (DoubleSolution) individual.copy());
                } else if (flag == 0) {
                    if (f2 < f1) {
                        population.set(k, (DoubleSolution) individual.copy());
                        time++;
                    }
                }
            } else {
                if (f2 < f1) {
                    population.set(k, (DoubleSolution) individual.copy());
                    time++;
                }
            }

            if (time >= maximumNumberOfReplacedSolutions) {
                return;
            }
        }
    }

    @Override
    public String getName() {
        return "cMOEAD";
    }

    @Override
    public String getDescription() {
        return "Multi-Objective Evolutionary Algorithm based on Decomposition with constraints support";
    }

}
/*probability.setSolutionBestAndBadByMethod(solutionBestAndBadByMethod);
                probability.calculateSumatoria();                
                probability.calculateProbabilityForEachMethod();
                probability.setViabilityByMethod();
                probability.fillTableProbability();
                
                _clearSolutionBestAndBadByMethod();
 */


 /* private int[] _getProbability(int[][] values) {
      int[] newValues = _getValues(values);
      int[] probability = _getProbability(newValues);
      return probability;
  }
  
  private int[] _getValues(int[][] values){
      int[] arr = new int[values.length + 1];
      int total = 0;
     
      int sum = 0;
      
      for (int i = 0; i < values.length; i++) {
          sum = values[i][0] + values[i][1];
          if(sum != 0)
            arr[i] = Math.round(values[i][0] * 100 / sum );          
              
          total += arr[i];
      }
      arr[values.length] = total;
      
      return arr;
  }
  
  private int[] _getProbability(int[] values){
      int total = values[values.length - 1];
       if(total == 0){
          System.out.println("dfgd");
      }
      int resto = 0, sum = 0;
      int[] probability = new int[values.length - 1];
      
      for (int i = 0; i < probability.length; i++) {
        probability[i] = Math.round( (values[i] * 100) / total );
        sum += probability[i];
      }
      resto = 100 - sum;
      int mayor = probability[0];
      int idx = 0;
      for (int i = 1; i < probability.length; i++) {
          if(probability[i] > mayor){
              mayor = probability[i];
              idx = i;
          }        
      }      
      
      probability[idx] += resto;
      
      probability = _updateProbability(probability);
      
      return probability;
  }*/
